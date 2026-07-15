

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
        eff_i *= 0.60
    elif gc_error <= 0.60:
        eff_i *= 0.40
    else:
        eff_i *= 0.20
        #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
        """OPTIONAL: Local GC Content affect"""

    sim_pF = Levenshtein.ratio(primer_F, seq_pF)
    sim_pR = Levenshtein.ratio(primer_R, seq_pR)

    # dropping efficiency if mismatches at the primer sites
    if (sim_pF and sim_pR) < 0.95:
        eff_i *= 0.70
    elif (sim_pF and sim_pR) < 0.90:
        eff_i *= 0.60
    elif (sim_pF and sim_pR) < 0.80:
        eff_i *= 0.50
    else:
        eff_i = eff_i
    return eff_i


# PCR pre-filtering --> removing non_specific amplicons as a cleanup for sequencing
dtype = [
    ("pid", np.int64),
    ("count", np.int64),
    ("seq", object),
]

with open(in_file, "r") as f:
    next(f, None)

    np_data = np.array(
        [
            (pid, count, seq)
            for line in f
            if "," in line
            for parts in [line.strip().split(",", maxsplit=2)]
            for pid in [int(parts[0].strip())]
            for count in [int(parts[1].strip())]
            for seq in [parts[2].strip()]
            if count > 0
        ],
        dtype=dtype,
    )

    """
    #200823,AATGGTTTACCCATA
    #count = [200823], seq [AATGGTTTACCCATA]
    pairs = [line.strip().split(",") for line in f if line.strip()]
    count, seq = zip(*pairs)
    """


seqs = np_data["seq"]
num_records = len(np_data)

similarity_pF = np.fromiter(
    (
        Levenshtein.ratio(primer_F, seq[:pF_length])
        for seq in seqs
    ),
    dtype=np.float64,
    count=num_records,
)

similarity_pR = np.fromiter(
    (
        Levenshtein.ratio(primer_R, seq[-pR_length:])
        for seq in seqs
    ),
    dtype=np.float64,
    count=num_records,
)

length_difference = np.fromiter(
    (
        abs(len(seq) - orig_length)
        for seq in seqs
    ),
    dtype=np.int64,
    count=num_records,
)

# True means the row is a PCR dropout.
dropout_mask = (
    ((similarity_pF < 0.75) & (similarity_pR < 0.75))
    | (length_difference > 20)
)

dropouts = np_data[dropout_mask]
filtered_lines = np_data[~dropout_mask]






primer_dropout_mask = (
    (similarity_pF < 0.75)
    & (similarity_pR < 0.75)
)

length_dropout_mask = length_difference > 20

print("\n--- PCR PRE-FILTER DEBUG ---")
print(f"Input oligos:                 {len(np_data)}")
print(f"Primer dropouts only:         "
      f"{np.count_nonzero(primer_dropout_mask & ~length_dropout_mask)}")
print(f"Length dropouts only:         "
      f"{np.count_nonzero(length_dropout_mask & ~primer_dropout_mask)}")
print(f"Primer and length dropouts:   "
      f"{np.count_nonzero(primer_dropout_mask & length_dropout_mask)}")
print(f"Total dropouts:               {np.count_nonzero(dropout_mask)}")
print(f"Surviving filtered oligos:    {np.count_nonzero(~dropout_mask)}")
print("----------------------------\n")






with open(fr"{PCR_DIR}/pcr_dropouts.txt", "w") as f:
    f.write("parent_id,count,sequence\n")

    for row in dropouts:
        f.write(
            f"{row['pid']},{row['count']},{row['seq']}\n"
        )

if filtered_lines.size == 0:
    raise ValueError("PCR FAILED!!!")

with open(fr"{PCR_DIR}/pcr_filtered.txt", "w") as f:
    f.write("parent_id,count,sequence\n")

    for row in filtered_lines:
        f.write(
            f"{row['pid']},{row['count']},{row['seq']}\n"
        )


# PCR sampling
sampled_data = filtered_lines.copy()

sampled_data["count"] = np.random.binomial(
    n=sampled_data["count"],
    p=sampling_frac
)

with open(fr"{PCR_DIR}/pcr_sampled.txt", "w") as f:
    f.write("parent_id,count,sequence\n")

    for row in sampled_data:
        f.write(
            f"{row['pid']},{row['count']},{row['seq']}\n"
        )

# PCR amplification
#num_cycles = int(input("Enter the desired PCR cycle to run: ")) # e.g 30
c1 = num_cycles // 3 # e.g 10
c2 = c1 + num_cycles // 3 # e.g 10 + 10 = 20
c3 = num_cycles # e.g 30

E0_i = 0.9  # initial efficiency
"""Calculating per oligo Efficiency in the sample based on [PRIMER BINDING, GC CONTENT] of the oligo"""

#sampled_data = sampled_data[sampled_data["count"] > 0]
amp_dtype = [
    ("pid", np.int64),
    ("count", np.float64),
    ("seq", object),
    ("eff", np.float64),
]
data = np.empty(len(sampled_data), dtype=amp_dtype)
data["pid"] = sampled_data["pid"]
data["count"] = sampled_data["count"]
data["seq"] = sampled_data["seq"]

data["eff"] = np.fromiter(
    (
        seq_check(E0_i, seq)
        for seq in sampled_data["seq"]
    ),
    dtype=np.float64,
    count=len(sampled_data),
)

starting_molecules = sampled_data["count"].sum(dtype=np.int64)

for c in range(0, c1):  # on good conditions, exponential efficiency in c1

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

        index_i['eff'] = eff  # this could store -ve efficiency, but we are dealing with it
        index_i['count'] = n_i

for c in range(c1, c2):  # efficiency become somewhat linear--> slow linear decrease

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
        # amp_i = 0.1 + eff  # very minimal amplification in plateau stage
        amp_i = 0.4 + eff  # very minimal amplification in plateau stage
        if amp_i < 1.0:  # preventing counts decrements incase eff is less than 0.8
            n_i = count_i
            delta = 0
        else:
            n_i = count_i * amp_i
            delta = n_i - count_i

        index_i['eff'] = eff
        index_i['count'] = n_i

# Writing back the copy counts after PCR completion
with open(fr'{PCR_DIR}/pcr_complete.txt', "w") as f:
    f.write("parent_id, sampled_count, sequence,efficiency_Remaining\n")
    for index in data:
        pid = index['pid']
        count = int(index['count'])
        seq = index['seq']
        eff = index['eff']
        f.write(f'{pid},{count},{seq},{eff:.5f}\n')

#Mutating oligo's and making several variants per oligo
# around 5-10 variants randomly
# keeping the original oligo sequence but distributing the copy counts across variants

initial_pid = data["pid"].copy()
initial_copies = data["count"].astype(np.int64, copy=True)
initial_lines = data["seq"].copy()
initial_eff = data["eff"].copy()

seq_objs = np.fromiter(
    (
        Error_simulation(
            seq, "pcr", attribute=mutation_attributes["1"], error_rate=err_rates["1"],
        )
        for seq in initial_lines
    ),
    dtype=object,
    count=len(initial_lines),
)


# Initially contains the unchanged sequences.
# This also correctly handles VALVE == 0.
MUTATED_TEXT = initial_lines.copy()

CHANGED_TEXT = []
UN_CHANGED_TEXT = []
UN_CHANGED_TEXT_02 = []


VALVE_HIGH = 15
VALVE_MED = 10 # Knob for mutations
VALVE_LOW = 5
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

for _ in range(VALVE):
    for i, seq_obj in enumerate(seq_objs):
        seq_obj.run_mutations()
        MUTATED_TEXT[i] = seq_obj.seq


# giving each original oligo sequence 80% copy count, rest of the count will go to mutated variants and CHIMERAS
unchanged_counts = (initial_copies * 0.80).astype(np.int64)

unchanged_dtype = [
    ("pid", np.int64),
    ("count", np.int64),
    ("seq", object),
]
UN_CHANGED_TEXT = np.empty(
    len(initial_lines),
    dtype=unchanged_dtype,
)

UN_CHANGED_TEXT["pid"] = initial_pid
UN_CHANGED_TEXT["count"] = unchanged_counts
UN_CHANGED_TEXT["seq"] = initial_lines

with open(fr"{PCR_DIR}/pcr_complete_2.txt", "w") as f:
    f.write("parent_id,count,sequence\n")

    for row in UN_CHANGED_TEXT:
        f.write(
            f"{row['pid']}"
            f"{row['count']}"
            f"{row['seq']}"
            f"\n"
        )


# Giving 10% copy count to mutated sequences.
# If no mutation occurred, return that 10% to the original sequence.

ten_percent_counts = (initial_copies * 0.10).astype(np.int64)

changed_mask = MUTATED_TEXT != initial_lines
unchanged_mask = ~changed_mask

variant_dtype = [
    ("pid", np.int64),
    ("count", np.int64),
    ("seq", object),
]

# Sequences that actually changed
CHANGED_TEXT = np.empty(
    np.count_nonzero(changed_mask),
    dtype=variant_dtype,
)


CHANGED_TEXT["pid"] = initial_pid[changed_mask]
CHANGED_TEXT["count"] = ten_percent_counts[changed_mask]
CHANGED_TEXT["seq"] = MUTATED_TEXT[changed_mask]

# Sequences that did not change
UN_CHANGED_TEXT_02 = np.empty(
    np.count_nonzero(unchanged_mask),
    dtype=variant_dtype,
)


UN_CHANGED_TEXT_02["pid"] = initial_pid[unchanged_mask]
UN_CHANGED_TEXT_02["count"] = ten_percent_counts[unchanged_mask]
UN_CHANGED_TEXT_02["seq"] = initial_lines[unchanged_mask]

with open(fr"{PCR_DIR}/pcr_CHANGED_POOL.txt", "w") as f:
    f.write("parent_id,count,sequence\n")

    for row in CHANGED_TEXT:
        f.write(
            f"{row['pid']}"
            f"{row['count']}"
            f"{row['seq']}\n"
        )





# Making different variants per oligo, distributing 10% copy_count
# NOTE: WE ARE DOING SUBSTITUTION ONLY as it is the predominant mutation error

allowed_bases = re.compile(r"[AGCT]", re.IGNORECASE)

variant_records = []

for parent_id, copy_count, original_line in zip(
    initial_pid,
    initial_copies,
    initial_lines,
):
    num_variant_line = random.randint(2, 5)
    random_mut_num = random.randint(1, 3)  # <----- we can control the amount of SUBS we want

    copy_count_line = int(int(copy_count) * 0.10)

    portions = np.random.multinomial(
        copy_count_line,
        [1 / num_variant_line] * num_variant_line,
    ) # distributes copy_count_line into x portions that sums upto the copy_count_line, GOOD!

    line = original_line

    while num_variant_line > 0:
        line = line.upper()

        for _ in range(0, random_mut_num + 1):
            random_base = random.choice(["A", "G", "C", "T"])
            pos = random.randrange(len(line))

            line = (
                line[:pos]
                + random_base
                + line[pos + 1:]
            )

        num_variant_line -= 1

        variant_records.append(
            (
                int(parent_id),
                int(portions[num_variant_line]),
                line,
            )
        )

 # pid, count, seq
LIST = np.array(
    variant_records,
    dtype=variant_dtype,
)



print(f"Length List = {len(LIST)+len(CHANGED_TEXT)+len(UN_CHANGED_TEXT)+len(UN_CHANGED_TEXT_02)}")

"""
LIST = [] # variants counts list with 10% of initial count
CHANGED_TEXT = [] # Mutated oligos with 10% of initial count
UN_CHANGED_TEXT = [] # original oligo with 80% of initial count
UN_CHANGED_TEXT_02 = [] # if not Mutated, 10% of initial count back to original
"""
with open(fr"{PCR_DIR}/pcr_pre_final.txt", "w") as f:
    f.write("parent_id,count,sequence,length\n")

    for item in LIST:
        f.write(
            f"{item['pid']}"
            f"{item['count']}"
            f"{item['seq']}"
             f"\n"
        )

    for item in CHANGED_TEXT:
        f.write(
            f"{item['pid']}"
            f"{item['count']}"
            f"{item['seq']}"
             f"\n"
        )

    for item in UN_CHANGED_TEXT:
        f.write(
            f"{item['pid']}"
            f"{item['count']}"
            f"{item['seq']}"
             f"\n"
        )

    for item in UN_CHANGED_TEXT_02:
        f.write(
            f"{item['pid']}"
            f"{item['count']}"
            f"{item['seq']}"
            f"\n"
        )

# WORKING ON CHIMERAS NOW!!
# reducing total copy count a little (taking 5% from pool), and distributing amongst Chimeras
# DEFAULT: Chimeras are 5% of the total pcr reads in our simulator, Change the knob value below to increase/decrease their quantity

pre_final_data = np.concatenate(
    [
        LIST.astype(variant_dtype, copy=False),
        CHANGED_TEXT.astype(variant_dtype, copy=False),
        UN_CHANGED_TEXT.astype(variant_dtype, copy=False),
        UN_CHANGED_TEXT_02.astype(variant_dtype, copy=False),
    ]
)

high_copy_mask = pre_final_data["count"] > 1000

LIST_02 = pre_final_data[high_copy_mask].copy()
LIST_03 = pre_final_data[~high_copy_mask].copy()


CHIMERAS_LIST = []

chimeras_variants = random.randint(10, 30) # <--- Knob for number of Chimeras variants

original_high_counts = LIST_02["count"].copy()

LIST_02["count"] = np.fromiter(
    (
        int(int(copy) * 0.95)
        for copy in original_high_counts
    ),
    dtype=np.int64,
    count=len(original_high_counts),
)

chimeras_copy_count = int( sum(  (5 / 100) * int(copy)  for copy in original_high_counts  ) ) #  <--- knob for Chimeras total counts (5% set here)

portions = np.random.multinomial( chimeras_copy_count,  [1 / chimeras_variants] * chimeras_variants, ) # distributing chimeras_copy_count into x portions that sums upto chimeras_copy_count

chimera_sources = list(
    zip(
        initial_pid.tolist(),
        initial_lines.tolist(),
    )
)

chimera_records = []

while chimeras_variants > 0:
    (pid_01, seq_01), (pid_02, seq_02) = random.sample( chimera_sources, 2, )

    while len(seq_01) != len(seq_02):
        (pid_01, seq_01), (pid_02, seq_02) = random.sample( chimera_sources,2, )

    pos = random.randrange(len(seq_01))

    new_chimera = ( seq_01[:pos] + seq_02[pos:] )

    chimeras_variants -= 1

    chimera_records.append(
        (
            int(pid_01),  # PID of the first parent
            int(portions[chimeras_variants]),
            new_chimera,
        )
    )


chimera_dtype = [
    ("pid", np.int64),
    ("count", np.int64),
    ("seq", object),
]

CHIMERAS_LIST = np.array(
    chimera_records,
    dtype=chimera_dtype,
)














def safe_count_sum(array):
    """Use Python integers to avoid possible np.int64 overflow."""
    return sum(int(value) for value in array["count"])


# Total molecules available after PCR amplification
available_total = sum(int(value) for value in initial_copies)

# Pre-chimera allocation
original_sequence_total = safe_count_sum(UN_CHANGED_TEXT)

error_simulation_total = (
    safe_count_sum(CHANGED_TEXT)
    + safe_count_sum(UN_CHANGED_TEXT_02)
)

substitution_variant_total = safe_count_sum(LIST)

pre_final_total = (
    original_sequence_total
    + error_simulation_total
    + substitution_variant_total
)


# Chimera-stage totals
high_count_before_chimeras = sum(
    int(value) for value in original_high_counts
)

high_count_after_chimeras = safe_count_sum(LIST_02)
low_count_total = safe_count_sum(LIST_03)
chimera_total = safe_count_sum(CHIMERAS_LIST)

final_total = (
    high_count_after_chimeras
    + low_count_total
    + chimera_total
)


# Losses caused only by integer rounding
split_rounding_loss = available_total - pre_final_total

chimera_rounding_loss = (
    high_count_before_chimeras
    - high_count_after_chimeras
    - chimera_total
)

total_loss = available_total - final_total


print("\n--- MOLECULE COUNT AUDIT ---")
print(f"Available after PCR:          {available_total}")
print(f"Original-sequence bucket:     {original_sequence_total}")
print(f"Error-simulation bucket:      {error_simulation_total}")
print(f"Substitution-variant bucket:  {substitution_variant_total}")
print(f"Pre-chimera total:            {pre_final_total}")
print(f"Split rounding loss:          {split_rounding_loss}")
print()
print(f"High-count before chimeras:   {high_count_before_chimeras}")
print(f"High-count after chimeras:    {high_count_after_chimeras}")
print(f"Chimera total:                {chimera_total}")
print(f"Low-count total:              {low_count_total}")
print(f"Chimera rounding loss:        {chimera_rounding_loss}")
print()
print(f"Final total:                  {final_total}")
print(f"Total molecule loss:          {total_loss}")
print("----------------------------\n")


# Consistency checks
assert pre_final_total <= available_total, (
    "Pre-chimera distribution exceeded available molecules"
)

assert (
    high_count_after_chimeras + chimera_total
    <= high_count_before_chimeras
), "Chimera distribution exceeded the high-count pool"

assert final_total <= available_total, (
    "Final molecule count exceeded amplified molecules"
)

assert pre_final_total == (
    high_count_before_chimeras + low_count_total
), "High/low partition does not match the pre-final pool"













with open(out_file, "w") as f:
    f.write("parent_id,count,sequence\n")

    # High-copy sequences after the 5% reduction
    for item in LIST_02:
        f.write(
            f"{item['pid']},"
            f"{item['count']},"
            f"{item['seq']}"
            f"\n"
        )

    # Chimera sequences
    for item in CHIMERAS_LIST:
        f.write(
            f"{item['pid']},"
            f"{item['count']},"
            f"{item['seq']}"
            f"\n"
        )

    for item in LIST_03:
        f.write(
            f"{item['pid']},"
            f"{item['count']},"
            f"{item['seq']}"
            f"\n"
        )

#print(f"final file length: {len(LIST_02)} + {len(LIST_03)} + {len(CHIMERAS_LIST)}")

combined_array = np.concatenate((LIST_02, LIST_03, CHIMERAS_LIST))
print(f"final file length: {len(combined_array)}")

sum_copies_pcr_final = (
    LIST_02["count"].sum(dtype=np.int64)
    + CHIMERAS_LIST["count"].sum(dtype=np.int64)
    + LIST_03["count"].sum(dtype=np.int64)
)

sum_initial_copies_pcr = initial_copies.sum(dtype=np.int64)

print(f"starting molecules: {int(starting_molecules)}")
#print(f"initial_copies_pcr: {int(sum_initial_copies_pcr)}")
print(f"sum_copies_pcr_final: {int(sum_copies_pcr_final)}")
print(
    f"Diff: "
    f"{int(sum_initial_copies_pcr - sum_copies_pcr_final)}"
)


if __name__ == "__main__":
    #os.remove(fr'{PCR_DIR}/pcr_dropouts.txt')
    #os.remove(fr'{PCR_DIR}/pcr_filtered.txt')
    #os.remove(fr'{PCR_DIR}/pcr_sampled.txt')
    #os.remove(fr'{PCR_DIR}/pcr_complete.txt')
    #os.remove(fr'{PCR_DIR}/pcr_complete_2.txt')
    #os.remove(fr'{PCR_DIR}/pcr_CHANGED_POOL.txt')
    #os.remove(fr'{PCR_DIR}/pcr_pre_final.txt')

    print("Pcr.py run completed")
