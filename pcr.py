import os
import random
import re

import numpy as np
from fontTools.svgLib.path.parser import UPPERCASE

from GC_content import gc_error_probability
from Enzyme_Addition import primer_F, primer_R,orig_length,pF_length,pR_length
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


def amp_factor(eff_i):
    return 1 + eff_i

def seq_check(eff_i, sequence):
    num_primer_mismatch = 0
    num_primer_mismatch += sum(1 for a, b in zip(sequence[:pF_length], primer_F) if a != b)
    num_primer_mismatch += sum(1 for a, b in zip(sequence[-pR_length:], primer_R) if a != b)

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

    # dropping efficiency if mismatches at the primer sites
    if num_primer_mismatch < 1:
        eff_i = eff_i
    elif num_primer_mismatch < 2:
        eff_i *= 0.7
    elif num_primer_mismatch < 3:
        eff_i *= 0.4
    else:
        eff_i *= 0.0  # the oligo won't amplify if primer mutation > 3
    return eff_i



path = fr'{os.getcwd()}\dna-fountain'

# PCR pre-filtering --> removing non_specific amplicons as a cleanup for sequencing
with open(fr'{path}\storage_file_1.txt') as f:
    """
    #200823,AATGGTTTACCCATA
    #count = [200823], seq [AATGGTTTACCCATA]
    pairs = [line.strip().split(",") for line in f if line.strip()]
    count, seq = zip(*pairs)
    """
    data = []
    for line in f:
        if "," in line.strip():
            count, seq = line.split(",")
            data.append({'count': int(count.strip()), 'seq': seq.strip()})

dropouts = []
filtered_lines = []

for index in data:
    count = index['count']
    seq = index['seq']
    if len(seq) == orig_length:
        if seq[:pF_length] != primer_F or seq[-pR_length:] != primer_R:
            dropouts.append((count,seq))
        else:
            filtered_lines.append((count,seq))
    else:
        if abs(len(seq) - orig_length) >= 10:
            dropouts.append((count,seq))
        else:
            filtered_lines.append((count,seq))

with open(fr'{path}\pcr_dropouts.txt', "w") as f:
    f.write("count, sequence, length\n")
    for c, s in dropouts:
        l = len(s)
        f.write(f"{c},{s},{l}\n")

with open(fr'{path}\pcr_filtered.txt', "w") as f:
    f.write("count, sequence, length\n")
    for c, s in filtered_lines:
        l = len(s)
        f.write(f"{c},{s},{l}\n")


# PCR sampling
sampling_frac = float(input("Enter PCR sampling fraction (in %age): ")) / 100

with open(fr'{path}\pcr_filtered.txt') as f:
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

with open(fr'{path}\pcr_sampled.txt', "w") as f:
    f.write("sampled_count, sequence, length\n")
    for index in data:
        count = index['count']
        seq = index['seq']
        length = index['length']
        f.write(f'{count},{seq},{length}\n')

# PCR amplification
num_cycles = int(input("Enter the desired PCR cycle to run: ")) # e.g 30
c1 = num_cycles // 3 # e.g 10
c2 = c1 + num_cycles // 3 # e.g 10 + 10 = 20
c3 = num_cycles # e.g 30

max_yield = int(input("Enter the maximum pcr yield expected:")) # so new pool will be original molecules + pcr pool/yield

E0_i = 0.9  # initial efficiency
"""Calculating per oligo Efficiency in the sample based on [PRIMER BINDING, GC CONTENT] of the oligo"""
with open(fr'{path}\pcr_sampled.txt') as f:
    next(f)
    data = []
    for line in f:
        if "," in line.strip():
            count, seq, length = line.split(",")
            data.append({ 'count': int(count.strip()), 'seq': seq.strip(), 'length': length.strip(), 'eff': seq_check(E0_i,seq) })

    remaining_yield = max_yield

    for c in range(0, c1): # on good conditions, exponential efficiency in c1
        if remaining_yield <= 0:
            break

        new_copy_count = []
        requested_copies = []

        for index_i in data:
            #eff_01 = E0_i
            count_i = index_i['count']
            seq_i = index_i['seq']
            eff = index_i['eff']
            # Turning sequence check at every cycle off as it might cause unnecessary efficiency collapse, instead it is computed once at the start
            #eff = seq_check(eff, seq_i)
            if not (0 <= eff <= 1):
                n_i = count_i
                delta = 0
            else:
                amp_i = amp_factor(eff)
                n_i = count_i * amp_i
                delta = n_i - count_i
            requested_copies.append(delta)
            new_copy_count.append(n_i)
            index_i['eff'] = eff # this could store -ve efficiency, but we are dealing with it
        total_demand = sum(requested_copies)
        if total_demand == 0: # nothing to amplify in this cycle
            break
        if remaining_yield >= total_demand:
            for i, item in enumerate(data):
                item['count'] = int(new_copy_count[i])
            remaining_yield -= total_demand
        else: # if demand higher than what we have, we first check if yield remains then distribute proportionally, otherwise break from the cycle
            if remaining_yield:
                sum_amount_given = 0
                for i, item in enumerate(data):
                    fraction = requested_copies[i] / total_demand
                    amount_given = fraction * remaining_yield
                    sum_amount_given += amount_given
                    item['count'] += int(amount_given)
                remaining_yield -= sum_amount_given
            else:
                break

    for c in range(c1, c2): # efficiency become somewhat linear--> slow linear decrease
        if remaining_yield <= 0:
            break

        new_copy_count = []
        requested_copies = []

        for index_i in data:
            #eff_02 = eff_01
            count_i = index_i['count']
            seq_i = index_i['seq']
            eff = index_i['eff']
            eff = eff - (eff / 10)
            #eff = seq_check(eff, seq_i)
            if not (0 <= eff <= 1):
                n_i = count_i
                delta = 0
            else:
                amp_i = amp_factor(eff)
                n_i = count_i * amp_i
                delta = n_i - count_i
            requested_copies.append(delta)
            new_copy_count.append(n_i)
            index_i['eff'] = eff
        total_demand = sum(requested_copies)
        if total_demand == 0:
            break
        if remaining_yield >= total_demand:
            for i, item in enumerate(data):
                item['count'] = int(new_copy_count[i])
            remaining_yield -= total_demand
        else:
            if remaining_yield:
                sum_amount_given = 0
                for i, item in enumerate(data):
                    fraction = requested_copies[i] / total_demand
                    amount_given = fraction * remaining_yield
                    sum_amount_given += amount_given
                    item['count'] += int(amount_given)
                remaining_yield -= sum_amount_given
            else:
                break

    for c in range(c2, c3):  # efficiency plateaus
        if remaining_yield <= 0:
            break

        new_copy_count = []
        requested_copies = []

        for index_i in data:
            #eff_03 = 0
            count_i = index_i['count']
            seq_i = index_i['seq']
            eff = index_i['eff']
            #eff = seq_check(eff, seq_i)
            amp_i = 0.2 + eff  # very minimal amplification in plateau stage
            if amp_i < 1.0: # preventing counts decrements incase eff is less than 0.8
                n_i = count_i
                delta = 0
            else:
                n_i = count_i * amp_i
                delta = n_i - count_i

            requested_copies.append(delta)
            new_copy_count.append(n_i)
            index_i['eff'] = eff
        total_demand = sum(requested_copies)
        if total_demand == 0:
            break

        if remaining_yield >= total_demand:
            for i, item in enumerate(data):
                item['count'] = int(new_copy_count[i])
            remaining_yield -= total_demand
        else:
            if remaining_yield:
                sum_amount_given = 0
                for i, item in enumerate(data):
                    fraction = requested_copies[i] / total_demand
                    amount_given = int(fraction * remaining_yield)
                    sum_amount_given += amount_given
                    item['count'] += amount_given
                remaining_yield -= sum_amount_given
            else:
                break

# Writing back the copy counts after PCR completion
with open(fr'{path}\pcr_complete.txt', "w") as f:
    f.write("sampled_count, sequence, length, efficiency_Remaining\n")
    for index in data:
        count = index['count']
        seq = index['seq']
        length = index['length']
        eff = index['eff']
        f.write(f'{count},{seq},{length},{eff:.5f}\n')

#Mutating oligo's and making several variants per oligo
# around 5-10 variants randomly
# keeping the original oligo sequence but distributing the copy counts across variants

MUTATED_TEXT = []

with open(fr'{path}\pcr_complete.txt') as f:
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

count = 1
while count < 30:
    MUTATED_TEXT.clear()
    for sE in seq_objs:
        #sE.reset_visited() # See important Notice for info
        sE.run_mutations()
        MUTATED_TEXT.append(sE.seq)
    count += 1

# giving each original oligo sequence 80% copy count, rest of the count will go to mutated variants and CHIMERAS
with open(fr'{path}\pcr_complete_2.txt', "w") as f:
    f.write("sampled_count, sequence, length\n")
    for copy_count, line, length in zip(initial_copies, initial_lines,initial_length):
        UN_CHANGED_TEXT.append((int(int(copy_count) * 0.80), line))
        f.write(f"{int(int(copy_count) * 0.80)},{line},{length}\n")

# giving each mutated sequence(sequences generated via Error_module.py) a 10% copy count of original sequence
with open(fr'{path}\pcr_CHANGED_POOL.txt', "w") as f:
    f.write("count, sequence, length\n")
    for copy_count, line, length in zip(initial_copies, initial_lines, initial_length):
        if line not in MUTATED_TEXT:
            CHANGED_TEXT.append((int(int(copy_count) * 0.10),line))
            f.write(f"{int(int(copy_count) * 0.10)},{line},{length}\n")
        else:
            UN_CHANGED_TEXT_02.append((int(int(copy_count) * 0.10),line)) # if a sequence is not mutated, adding the 10% count back to the original sequence


# Making different variants per oligo, distributing 10% copy_count
# NOTE: WE ARE DOING SUBSTITUTION ONLY as it is the predominant mutation error
allowed_bases = re.compile(r"[AGCT]", re.IGNORECASE)
LIST = []
for copy_count, line in zip(initial_copies, initial_lines):
    num_variant_line = random.randint(2,5)
    random_mut_num = random.randint(1,3) # <----- we can control the amount of SUBS we want
    line = line.upper()
    copy_count_line = int(int(copy_count) * 0.10)
    #print(f"copy_count_line: {copy_count_line}")
    portions = np.random.multinomial(copy_count_line, [1 / num_variant_line] * num_variant_line) # distributes copy_count_line into x portions that sums upto the copy_count_line, GOOD!
    #print(f"Portions: {portions}")
    while num_variant_line > 0:
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
with open(fr'{path}\pcr_pre_final.txt', "w") as f:
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
with open(fr'{path}\pcr_pre_final.txt') as f:
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


with open(fr'{path}\pcr_final.txt',"w") as f:
    f.write("count, sequence, length\n")
    for item in LIST_02:
        f.write(f"{item[0]},{item[1]},{len(item[1])}\n")
    for item in CHIMERAS_LIST:
        f.write(f"{item[0]},{item[1]},{len(item[1])}\n")
    for item in LIST_03:
        f.write(f"{item[0]},{item[1]},{len(item[1])}\n")

print(f"final file length: {len(LIST_02)} + {len(LIST_03)} + {len(CHIMERAS_LIST)}")

with open(fr'{path}\pcr_final.txt') as f:
    next(f)
    rows = [line.strip().split(",") for line in f if line.strip()]
    copy_count, _, _ = zip(*rows)
    sum_copies_pcr_final = sum([int(x) for x in copy_count])

sum_initial_copies_pcr = sum([int(x) for x in initial_copies])
print(f"initial_copies_pcr: {sum_initial_copies_pcr}")
print(f"sum_copies_pcr_final: {sum_copies_pcr_final}")
print(f"Diff:{sum_initial_copies_pcr - sum_copies_pcr_final} ")


os.remove(fr'{path}\pcr_complete.txt')
os.remove(fr'{path}\pcr_complete_2.txt')
os.remove(fr'{path}\pcr_CHANGED_POOL.txt')
os.remove(fr'{path}\pcr_CHANGED_POOL_02.txt')
os.remove(fr'{path}\pcr_pre_final.txt')




if __name__ == "__main__":
    print("Done")