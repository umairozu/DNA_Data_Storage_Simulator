

import argparse
import os
import random
import re

import Levenshtein
import numpy as np
import json
from GC_content import gc_error_probability
from Error_module import Error_simulation

err_rates = {
             "1": {"raw_rate": 0.000043, "substitution": 0.99, "deletion": 0.01, "insertion": 0}, #Taq
             "2": {"raw_rate": 0.0000024, "substitution": 1, "deletion": 0, "insertion": 0}, #Pwo
             "3": {"raw_rate": 0.0000028, "substitution": 1, "deletion": 0, "insertion": 0}, #Pfu
             "4": {"raw_rate": 0.0000026, "substitution": 0.84, "deletion": 0.08, "insertion": 0.08}, #Phusion
             "5": {"deletion": 0.3333, "substitution": 0.33340000000000003, "insertion": 0.3333, "raw_rate": 0.0} #None
             }

mutation_attributes = {
                        #Taq
                        "1": {"deletion": {"position": {"homopolymer": 0, "random": 1},
                                           "pattern": {"G": 0.25, "C": 0.25, "A": 0.25, "T": 0.25}},
                              "insertion": {"position": {"homopolymer": 0, "random": 1},
                                            "pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},
                              "substitution": {"pattern": {"A": {"G": 0.97, "T": 0.01, "C": 0.02},
                                                       "T": {"C": 0.97, "A": 0.01, "G": 0.02},
                                                       "G": {"A": 1, "T": 0, "C": 0},
                                                       "C": {"T": 1, "G": 0, "A": 0}}}},

                        #Pwo
                        "2": {"deletion": {"position": {"homopolymer": 0, "random": 1},
                                           "pattern": {"G": 0.25, "C": 0.25, "A": 0.25, "T": 0.25}},
                              "insertion": {"position": {"homopolymer": 0, "random": 1},
                                            "pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},
                              "substitution": {"pattern": {"A": {"G": 1, "T": 0, "C": 0},
                                                       "T": {"C": 0.67, "A": 0.33, "G": 0},
                                                       "G": {"A": 0.57, "T": 0, "C": 0.43},
                                                       "C": {"T": 1, "G": 0, "A": 0}}}},

                        #Pfu
                        "3": {"deletion": {"position": {"homopolymer": 0, "random": 1},
                                           "pattern": {"G": 0.25, "C": 0.25, "A": 0.25, "T": 0.25}},
                              "insertion": {"position": {"homopolymer": 0, "random": 1},
                                            "pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},
                              "substitution": {"pattern": {"A": {"G": 0.75, "T": 0.25, "C": 0},
                                                       "T": {"C": 0.75, "A": 0.25, "G": 0},
                                                       "G": {"A": 1, "T": 0, "C": 0},
                                                       "C": {"T": 1, "G": 0, "A": 0}}}},

                        #Phusion
                        "4": {"deletion": {"position": {"homopolymer": 0, "random": 1},
                                           "pattern": {"G": 0.25, "C": 0.25, "A": 0.25, "T": 0.25}},
                              "insertion": {"position": {"homopolymer": 0, "random": 1},
                                            "pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},
                              "substitution": {"pattern": {"A": {"G": 1, "T": 0, "C": 0},
                                                       "T": {"C": 1, "A": 0, "G": 0},
                                                       "G": {"A": 1, "T": 0, "C": 0},
                                                       "C": {"T": 1, "G": 0, "A": 0}}}},

                        #None
                        "5": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                           "position": {"homopolymer": 0.5, "random": 0.5}},
                              "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                            "position": {"homopolymer": 0.5, "random": 0.5}},
                              "substitution": {"pattern": {}}}
                        }

BASE_DIR = fr'{os.getcwd()}'
os.makedirs(fr'{BASE_DIR}/pcr', exist_ok=True)
STORAGE_DIR = fr'{os.getcwd()}/storage'
PCR_DIR = fr'{os.getcwd()}/pcr'

metadata_path = fr'{BASE_DIR}/metadata.json'
with open(metadata_path, "r") as f:
    metadata = json.load(f)

primer_F = metadata["primer_F"]
primer_R = metadata["primer_R"]
orig_length = metadata["orig_length"]
pF_length = metadata["pF_length"]
pR_length = metadata["pR_length"]



_p = argparse.ArgumentParser(description="pcr.py run")
_p.add_argument("--s", type = lambda x: float(x) if 0.0 < float(x) <= 100.0 else _p.error("keep 0 < Sampling fraction <= 100"),
                                        default = "1.0" , help = "sampling fraction > 0.0")
_p.add_argument("--n", type = lambda x: int(x) if int(x) > 0 else _p.error("Pcr cycle should be > 0 [better to keep it 30]"),
                                        default = "30" , help = "number of pcr cycles > 0")
_p.add_argument("--mut", default = "2" , help = "Mutation intensity (0-3)")
_p.add_argument("--c", type = lambda x: int(x) if  0 < int(x) <= 100 else _p.error("keep 0 < custom VALVE <= 100"),
                                        help = "Optional custom VALVE (0-100), it is basically a mutation knob" )
_p.add_argument("--in_file", required = True, help="Input file name (e.g  storage_in)")
_p.add_argument("--out_file", default = "pcr_output.txt", help = "Output file name [default is 'pcr_output'] ")



args = _p.parse_args()
sampling_frac = float(args.s) / 100
num_cycles = int(args.n)
in_file = fr'{args.in_file}'
out_file = fr'{PCR_DIR}/{args.out_file}'


def amp_factor(eff_i):
    return 1 + eff_i

def seq_check(eff_i, sequence):
    seq_pF = sequence[:pF_length]
    seq_pR = sequence[-pR_length:]

    # dropping efficiency if GC deviates from 45% - 55% range
    gc_error = gc_error_probability(sequence)
    if gc_error <= 0.0:
        pass
    elif gc_error <= 0.10:
        eff_i *= 0.90
    elif gc_error <= 0.25:
        eff_i *= 0.75
    elif gc_error <= 0.40:
        eff_i *= 0.55
    elif gc_error <= 0.60:
        eff_i *= 0.35
    else:
        eff_i *= 0.20
        #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
        """OPTIONAL: Local GC Content affect"""

    sim_pF = Levenshtein.ratio(primer_F, seq_pF)
    sim_pR = Levenshtein.ratio(primer_R, seq_pR)

    # dropping efficiency if mismatches at the primer sites
    if (sim_pF and sim_pR) < 0.90:
        eff_i *= 0.7
    elif (sim_pF and sim_pR) < 0.80:
        eff_i *= 0.6
    elif (sim_pF and sim_pR) < 0.70:
        eff_i *= 0.5
    else:
        eff_i = eff_i
    return eff_i


# PCR pre-filtering --> removing non_specific amplicons as a cleanup for sequencing
with open(fr'{in_file}') as f:
    next(f)
    data = []
    for line in f:
        if "," in line:
            count, seq, _ = line.split(",")
            count_int = int(count.strip())

            if count_int > 0:  # Only append if count is positive
                data.append({'count': count_int, 'seq': seq.strip()})
    """
    #200823,AATGGTTTACCCATA
    #count = [200823], seq [AATGGTTTACCCATA]
    pairs = [line.strip().split(",") for line in f if line.strip()]
    count, seq = zip(*pairs)
    """

dropouts = []
filtered_lines = []

for index in data:
    count = index['count']
    seq = index['seq']

    pF = seq[:pF_length]
    pR = seq[-pR_length:]
    primer_F = primer_F
    primer_R = primer_R

    similarity_pF = Levenshtein.ratio(primer_F, pF)
    similarity_pR = Levenshtein.ratio(primer_R, pR)
    #print(fr"Similarity  -->  pF: {similarity_pF}, pR: {similarity_pR}")

    if (similarity_pF < 0.75 and similarity_pR < 0.75) or (abs(len(seq) - orig_length) > 20):
        dropouts.append((count, seq))
    else:
        filtered_lines.append((count, seq))


with open(fr'{PCR_DIR}/pcr_dropouts.txt', "w") as f:
    f.write("count, sequence, length\n")
    for c, s in dropouts:
        l = len(s)
        f.write(f"{c},{s},{l}\n")

if len(filtered_lines) == 0:
    raise ValueError("PCR FAILED!!!")


with open(fr'{PCR_DIR}/pcr_filtered.txt', "w") as f:
    f.write("count, sequence, length\n")
    for c, s in filtered_lines:
        l = len(s)
        f.write(f"{c},{s},{l}\n")


# PCR sampling
#sampling_frac = float(input("Enter PCR sampling fraction (in %age): ")) / 100

with open(fr'{PCR_DIR}/pcr_filtered.txt') as f:
    next(f)
    data = []
    for line in f:
        if "," in line.strip():
            count, seq, length = line.split(",")
            data.append({ 'count': int(count.strip()), 'seq': seq.strip(), 'length': length.strip() })

for index in data:
    count = index['count']
    seq = index['seq']
    length = index['length']
    bnd = np.random.binomial(int(count),sampling_frac, size = int(length))
    bnd_choice = np.random.choice(bnd)
    index['count'] = bnd_choice

with open(fr'{PCR_DIR}/pcr_sampled.txt', "w") as f:
    f.write("sampled_count, sequence, length\n")
    for index in data:
        count = index['count']
        seq = index['seq']
        length = index['length']
        f.write(f'{count},{seq},{length}\n')

# PCR amplification
#num_cycles = int(input("Enter the desired PCR cycle to run: ")) # e.g 30
c1 = num_cycles // 3 # e.g 10
c2 = c1 + num_cycles // 3 # e.g 10 + 10 = 20
c3 = num_cycles # e.g 30

#max_yield = int(input("Enter the maximum pcr yield expected:")) # so new pool will be original molecules + pcr pool/yield

E0_i = 0.9  # initial efficiency
"""Calculating per oligo Efficiency in the sample based on [PRIMER BINDING, GC CONTENT] of the oligo"""
with open(fr'{PCR_DIR}/pcr_sampled.txt') as f:
    next(f)
    data = []
    for line in f:
        if "," in line.strip():
            count, seq, length = line.split(",")
            data.append({ 'count': int(count.strip()), 'seq': seq.strip(), 'length': length.strip(), 'eff': seq_check(E0_i,seq) })

    starting_molecules = sum(item['count'] for item in data)

    for c in range(0, c1): # on good conditions, exponential efficiency in c1

        for index_i in data:
            count_i = index_i['count']
            seq_i = index_i['seq']
            eff = index_i['eff']

            if not (0 <= eff <= 1):
                n_i = count_i
                delta = 0
            else:
                amp_i = amp_factor(eff)
                n_i = count_i * amp_i
                delta = n_i - count_i

            index_i['eff'] = eff # this could store -ve efficiency, but we are dealing with it
            index_i['count'] = n_i

    for c in range(c1, c2): # efficiency become somewhat linear--> slow linear decrease

        for index_i in data:
            count_i = index_i['count']
            seq_i = index_i['seq']
            eff = index_i['eff']
            eff = eff - (eff / 10)
            if not (0 <= eff <= 1):
                n_i = count_i
                delta = 0
            else:
                amp_i = amp_factor(eff)
                n_i = count_i * amp_i
                delta = n_i - count_i

            index_i['eff'] = eff
            index_i['count'] = n_i

    for c in range(c2, c3):  # efficiency plateaus

        for index_i in data:
            count_i = index_i['count']
            seq_i = index_i['seq']
            eff = index_i['eff']
            amp_i = 0.2 + eff  # very minimal amplification in plateau stage
            if amp_i < 1.0: # preventing counts decrements incase eff is less than 0.8
                n_i = count_i
                delta = 0
            else:
                n_i = count_i * amp_i
                delta = n_i - count_i

            index_i['eff'] = eff
            index_i['count'] = n_i


# Writing back the copy counts after PCR completion
with open(fr'{PCR_DIR}/pcr_complete.txt', "w") as f:
    f.write("sampled_count, sequence, length, efficiency_Remaining\n")
    for index in data:
        count = int(index['count'])
        seq = index['seq']
        length = index['length']
        eff = index['eff']
        f.write(f'{count},{seq},{length},{eff:.5f}\n')

#Mutating oligo's and making several variants per oligo
# around 5-10 variants randomly
# keeping the original oligo sequence but distributing the copy counts across variants

MUTATED_TEXT = []

with open(fr'{PCR_DIR}/pcr_complete.txt') as f:
    next(f)
    rows = [line.strip().split(",") for line in f if line.strip()]
    initial_copies, initial_lines, initial_length, eff = zip(*rows)

    seq_objs = [Error_simulation(seq, "pcr", attribute = mutation_attributes["1"],
                                     error_rate = err_rates["1"])
                    for seq in initial_lines
                    ]

CHANGED_TEXT = []
UN_CHANGED_TEXT = []
UN_CHANGED_TEXT_02 = []


VALVE_HIGH = 35
VALVE_MED = 30 # Knob for mutations
VALVE_LOW = 20
VALVE_NULL = 0

if args.c is not None:
    VALVE = args.c
    print(f"Using custom VALVE")
elif args.mut in ["0", "1", "2", "3"]:
    if args.mut == "0":
        VALVE = VALVE_NULL
    elif args.mut == "1":
        VALVE = VALVE_LOW
    elif args.mut == "2":
        VALVE = VALVE_MED
    else:
        VALVE = VALVE_HIGH
else:
    raise ValueError("Invalid mutation VALVE [try --help]")

print(f"Running with VALVE: {VALVE}")

if VALVE == 0:
    for sE in seq_objs:
        MUTATED_TEXT.append(sE.seq)

count = 0
while count < VALVE:
    MUTATED_TEXT.clear()
    for sE in seq_objs:
        #sE.reset_visited() # See important Notice for info
        sE.run_mutations()
        MUTATED_TEXT.append(sE.seq)
    count += 1

# giving each original oligo sequence 80% copy count, rest of the count will go to mutated variants and CHIMERAS
with open(fr'{PCR_DIR}/pcr_complete_2.txt', "w") as f:
    f.write("sampled_count, sequence, length\n")
    for copy_count, line, length in zip(initial_copies, initial_lines,initial_length):
        UN_CHANGED_TEXT.append((int(int(copy_count) * 0.80), line))
        f.write(f"{int(int(copy_count) * 0.80)},{line},{length}\n")

# giving each mutated sequence(sequences generated via Error_module.py) a 10% copy count of original sequence
with open(fr'{PCR_DIR}/pcr_CHANGED_POOL.txt', "w") as f:
    f.write("count, sequence, length\n")
    for copy_count, original_line, length in zip(initial_copies, initial_lines, initial_length):
        if original_line not in MUTATED_TEXT:
            mutated_line = original_line
            CHANGED_TEXT.append((int(int(copy_count) * 0.10),mutated_line)) # CORRECTIONNNNNN , handle deletion case
            f.write(f"{int(int(copy_count) * 0.10)},{line},{length}\n")
        else:
            UN_CHANGED_TEXT_02.append((int(int(copy_count) * 0.10),original_line)) # if a sequence is not mutated, adding the 10% count back to the original sequence


# Making different variants per oligo, distributing 10% copy_count
# NOTE: WE ARE DOING SUBSTITUTION ONLY as it is the predominant mutation error
allowed_bases = re.compile(r"[AGCT]", re.IGNORECASE)
LIST = []
for copy_count, line in zip(initial_copies, initial_lines):
    num_variant_line = random.randint(2,5)
    random_mut_num = random.randint(1,3) # <----- we can control the amount of SUBS we want
    copy_count_line = int(int(copy_count) * 0.10)
    #print(f"copy_count_line: {copy_count_line}")
    portions = np.random.multinomial(copy_count_line, [1 / num_variant_line] * num_variant_line) # distributes copy_count_line into x portions that sums upto the copy_count_line, GOOD!
    #print(f"Portions: {portions}")
    while num_variant_line > 0:
        line = line.upper()
        for i in range(0, random_mut_num + 1):
            random_base = random.choice(['A', 'G', 'C', 'T'])
            pos = random.randrange(len(line))
            line = line[:pos] + random_base + line[pos + 1:]
        num_variant_line -= 1
        LIST.append((portions[num_variant_line],line))


print(f"Length List = {len(LIST)+len(CHANGED_TEXT)+len(UN_CHANGED_TEXT)+len(UN_CHANGED_TEXT_02)}")
"""
LIST = [] # variants counts list with 10% of initial count
CHANGED_TEXT = [] # Mutated oligos with 10% of initial count
UN_CHANGED_TEXT = [] # original oligo with 80% of initial count
UN_CHANGED_TEXT_02 = [] # if not Mutated, 10% of initial count back to original
"""
with open(fr'{PCR_DIR}/pcr_pre_final.txt', "w") as f:
    f.write("count, sequence, length\n")
    for item in LIST:
        f.write(f"{item[0]},{item[1]},{len(item[1])}\n")
    for item in CHANGED_TEXT:
        f.write(f"{item[0]},{item[1]},{len(item[1])}\n")
    for item in UN_CHANGED_TEXT:
        f.write(f"{item[0]},{item[1]},{len(item[1])}\n")
    for item in UN_CHANGED_TEXT_02:
        f.write(f"{item[0]},{item[1]},{len(item[1])}\n")

# WORKING ON CHIMERAS NOW!!
# reducing total copy count a little (taking 5% from pool), and distributing amongst Chimeras
# DEFAULT: Chimeras are 5% of the total pcr reads in our simulator, Change the knob value below to increase/decrease their quantity
LIST_02 = []
LIST_03 = []
with open(fr'{PCR_DIR}/pcr_pre_final.txt') as f:
    next(f)
    rows = [line.strip().split(",") for line in f if line.strip()]
    copy_count, lines, _ = zip(*rows)
    num_lines = len(lines)

for copy, line in zip(copy_count, lines):
    copy = int(copy)
    if copy > 1000:
        LIST_02.append((copy,line))
    else:
        LIST_03.append((copy,line))

CHIMERAS_LIST = []
chimeras_variants = random.randint(10,30) # <--- Knob for number of Chimeras variants
chimeras_copy_count = 0
for i in range(len(LIST_02)):
    old_tuple = LIST_02[i]
    copy = int(old_tuple[0])
    new_copy_val = int(copy * 0.95)
    LIST_02[i] = (new_copy_val, old_tuple[1])
    chimeras_copy_count += ((5/100) * copy) #  <--- knob for Chimeras total counts (5% set here)

portions = np.random.multinomial(chimeras_copy_count, [1 / chimeras_variants] * chimeras_variants) # distributing chimeras_copy_count into x portions that sums upto chimeras_copy_count

while chimeras_variants > 0:
    seq_01, seq_02 = random.sample(initial_lines, 2)
    while len(seq_01) != len(seq_02):
        seq_01, seq_02 = random.sample(initial_lines, 2)
    pos = random.randrange(len(seq_01))
    new_chimeras = seq_01[:pos] + seq_02[pos:]
    chimeras_variants -= 1
    CHIMERAS_LIST.append((portions[chimeras_variants],new_chimeras))


with open(fr'{out_file}', "w") as f:
    f.write("count, sequence, length\n")
    for item in LIST_02:
        f.write(f"{item[0]},{item[1]},{len(item[1])}\n")
    for item in CHIMERAS_LIST:
        f.write(f"{item[0]},{item[1]},{len(item[1])}\n")
    for item in LIST_03:
        f.write(f"{item[0]},{item[1]},{len(item[1])}\n")

print(f"final file length: {len(LIST_02)} + {len(LIST_03)} + {len(CHIMERAS_LIST)}")

with open(fr'{out_file}') as f:
    next(f)
    rows = [line.strip().split(",") for line in f if line.strip()]
    copy_count, _, _ = zip(*rows)
    sum_copies_pcr_final = sum([int(x) for x in copy_count])

sum_initial_copies_pcr = sum([int(x) for x in initial_copies])

print(fr"starting molecules: {starting_molecules}")
print(f"initial_copies_pcr: {sum_initial_copies_pcr}")
print(f"sum_copies_pcr_final: {sum_copies_pcr_final}")
print(f"Diff:{sum_initial_copies_pcr - sum_copies_pcr_final} ")



if __name__ == "__main__":
    os.remove(fr'{PCR_DIR}/pcr_dropouts.txt')
    os.remove(fr'{PCR_DIR}/pcr_filtered.txt')
    os.remove(fr'{PCR_DIR}/pcr_sampled.txt')
    os.remove(fr'{PCR_DIR}/pcr_complete.txt')
    os.remove(fr'{PCR_DIR}/pcr_complete_2.txt')
    os.remove(fr'{PCR_DIR}/pcr_CHANGED_POOL.txt')
    os.remove(fr'{PCR_DIR}/pcr_pre_final.txt')

    print("Pcr.py run completed")
