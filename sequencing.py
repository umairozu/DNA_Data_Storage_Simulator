

# Error_rates taken from mesa
# https://github.com/umr-ds/mesa_dna_sim/blob/master/simulators/sequencing/sequencing_error.py

"""
Error rates origin

- 1 & 2 = Illumina     1 -> method == 'Single End', 2 -> method == 'Paired End'
- 3 & 4 = PacBio       3 -> method == 'CCS', 4 -> method == 'Subread'
- 5 & 6 = Nanopore     5 -> method == '1D' , 6 -> method == '2D'

mutation attributes origin
- 1 = Illumina
- 2 = PacBio
- 3 = Nanopore

"""
import argparse
import json
import os
import random

import numpy as np
from Error_module import Error_simulation
import Levenshtein
from Bio import Align
from Bio.Seq import Seq


err_rates = {"1": {"raw_rate": 0.0021, "substitution": 0.996, "deletion": 0.0024, "insertion": 0.0013},
             "2": {"raw_rate": 0.0032, "substitution": 0.997, "deletion": 0.0018, "insertion": 0.0011},
             #"1": {"raw_rate": 0.0010, "substitution": 0.995, "deletion": 0.003, "insertion": 0.003},# Illumia MiSeq SE, PE from Grass et al.
             "3": {"raw_rate": 0.02, "substitution": 0.75, "deletion": 0.20, "insertion": 0.05},
             "4": {"raw_rate": 0.14, "substitution": 0.37, "deletion": 0.21, "insertion": 0.42},
             "5": {"raw_rate": 0.2, "substitution": 0.48, "deletion": 0.37, "insertion": 0.15},
             "6": {"raw_rate": 0.13, "substitution": 0.41, "deletion": 0.36, "insertion": 0.23}}

mutation_attributes = {"1": {"deletion": {"position": {"random": 1},
                                          "pattern": {"G": 0.25, "C": 0.25, "A": 0.25, "T": 0.25}},
                             "insertion": {"position": {"random": 1},
                                           "pattern": {"G": 0.25, "C": 0.25, "A": 0.25, "T": 0.25}},
                             "substitution": {"pattern": {"A": {"G": 0.50, "T": 0.25, "C": 0.25},
                                                      "T": {"G": 0.50, "A": 0.25, "C": 0.25},
                                                      "C": {"G": 0.50, "A": 0.25, "T": 0.25},
                                                      "G": {"T": 0.50, "A": 0.25, "C": 0.25}}}},
                       "2": {"deletion": {"position": {"homopolymer": 0.85, "random": 0.15},
                                          "pattern": {"G": 0.35, "C": 0.35, "A": 0.15, "T": 0.15}},
                             "insertion": {"position": {"homopolymer": 0.85, "random": 0.15},
                                           "pattern": {"A": 0.35, "T": 0.35, "C": 0.15, "G": 0.15}},
                             "substitution": {"pattern": {"CG": {"CA": 0.5, "TG": 0.5}}}},
                       "3": {"deletion": {"position": {"homopolymer": 0.46, "random": 0.54},
                                          "pattern": {"G": 0.35, "C": 0.35, "A": 0.15, "T": 0.15}},
                             "insertion": {"position": {"homopolymer": 0.46, "random": 0.54},
                                           "pattern": {"A": 0.35, "T": 0.35, "C": 0.15, "G": 0.15}},
                             "substitution": {"pattern": {"TAG": "TGG", "TAC": "TGC"}}}}

BASE_DIR = fr'{os.getcwd()}'
os.makedirs(fr'{BASE_DIR}/sequencing', exist_ok=True)
PCR_DIR = fr'{os.getcwd()}/pcr'
SEQ_DIR = fr'{os.getcwd()}/sequencing'

metadata_path = fr'{BASE_DIR}/metadata.json'
with open(metadata_path, "r") as f:
    metadata = json.load(f)

primer_F = metadata["primer_F"]
primer_R = metadata["primer_R"]
orig_length = metadata["orig_length"]
pF_length = metadata["pF_length"]
pR_length = metadata["pR_length"]
payload_length = metadata["payload_length"]


def get_clean_segment(parent_aligned, child_aligned, start_pos, target_len):
    valid_chars = 0
    chars_consumed = 0

    for i in range(start_pos, len(parent_aligned)):
        if valid_chars == target_len:
            break

        chars_consumed += 1
        if parent_aligned[i] != '-':
            valid_chars += 1

    end_pos = start_pos + chars_consumed

    clean_segment = "".join(child_aligned[start_pos:end_pos]).replace('-', '')

    return clean_segment, end_pos


_p = argparse.ArgumentParser(description="sequencing.py run")
_p.add_argument("--type", type = lambda x: x if x == "1" else _p.error("type 1 for Illumina sequencing "),
                                        default = 1, help = "Sequencing technology: 1. Illumina")
_p.add_argument("--m", type = lambda x: x if x in ["1","2"] else _p.error("Illumina's sequencing methods \n1. Single-end, 2. Paired-end"),
                                        default = 2, help = "Illumina's sequencing methods \n1. Single-end, 2. Paired-end")
_p.add_argument("--s", type = lambda x: float(x) if 0.0 < float(x) <= 100.0 else _p.error("keep 0 < Sampling fraction <= 100"),
                                        default = "1.0" , help = "sampling fraction > 0.0")
_p.add_argument("--t", type = lambda x: int(x) if int(x) > 0 else _p.error("total reads needs to be > 0"),
                                            required = True, help = "total number of sequencing read required")
_p.add_argument("--rl", type = lambda x: int(x) if int(x) > 0 else _p.error("read length needs to be > 0"),
                                        required = True, help = "read length for sequencing")

_p.add_argument("--mut", default = "2" , help = "Mutation intensity (0-3)")
_p.add_argument("--c", type = lambda x: int(x) if  0 < int(x) <= 10 else _p.error("keep 0 < custom VALVE <= 10"),
                                        help = "Optional custom VALVE (0-10), it is basically a mutation knob" )
_p.add_argument("--in_file", required = True, help="Input file name (e.g  pcr_in)")
_p.add_argument("--order_file", required = True, help="original dna order file")


args = _p.parse_args()
mode = args.type
method = args.m
SAMPLING_FRAC = float(args.s) / 100.0
TARGET_READS = int(args.t)
READ_LENGTH = int(args.rl)
in_file = fr'{args.in_file}'
order_file = fr'{args.order_file}'

MUTATED_TEXT = []

# PCR input format:
# parent_id,count,sequence

input_dtype = [
    ("pid", np.int64),
    ("count", np.int64),
    ("seq", object),
]

input_records = []

with open(in_file) as f:
    next(f, None)

    for line_number, line in enumerate(f, start=2):
        line = line.strip()

        if not line or "," not in line:
            continue

        parts = line.split(",", maxsplit=2)

        if len(parts) != 3:
            raise ValueError(
                f"Invalid PCR input at line {line_number}: {line!r}"
            )

        pid_text, count_text, seq_text = parts

        input_records.append(
            (
                int(pid_text.strip()),
                int(count_text.strip()),
                seq_text.strip(),
            )
        )

data = np.array( input_records, dtype=input_dtype,)

POOL_COUNT = int( data["count"].sum(dtype=np.int64) )  # total molecule in input pcr file

# Collecting original oligo's given for order (input is .dna_order file from Enzyme_Addition.py class)

order_oligos = []

with open(order_file) as file:
    for line_number, line in enumerate(file, start=1):
        line = line.strip()

        if not line:
            continue

        parts = line.split(",", maxsplit=1)

        if len(parts) != 2:
            raise ValueError(
                f"Invalid order-file format at line "
                f"{line_number}: {line!r}"
            )

        pid_text, sequence_text = parts

        try:
            pid = int(pid_text.strip())
        except ValueError:
            raise ValueError(
                f"Invalid parent ID at order-file line "
                f"{line_number}: {pid_text!r}"
            ) from None

        sequence = sequence_text.strip()

        if not sequence:
            raise ValueError(
                f"Missing DNA sequence at order-file line "
                f"{line_number}"
            )

        # Enforce PID 0 = array position 0, PID 1 = position 1, etc.
        expected_pid = len(order_oligos)

        if pid != expected_pid:
            raise ValueError(
                f"Unexpected parent ID at order-file line "
                f"{line_number}: expected {expected_pid}, got {pid}"
            )

        order_oligos.append(sequence)


ORIGINAL_ORDER_OLIGOS = np.array(
    order_oligos,
    dtype=object,
)


invalid_pid_mask = (
    (data["pid"] < 0)
    | (data["pid"] >= len(ORIGINAL_ORDER_OLIGOS))
)

if np.any(invalid_pid_mask):
    invalid_pids = np.unique(
        data["pid"][invalid_pid_mask]
    )

    raise ValueError(
        "PCR output contains parent IDs that do not exist in "
        f"the order file. Order-file oligos: "
        f"{len(ORIGINAL_ORDER_OLIGOS)}, "
        f"invalid PID examples: "
        f"{invalid_pids[:10].tolist()}"
    )


parent_raw = ORIGINAL_ORDER_OLIGOS[data["pid"]]

similarities = np.fromiter(
    (
        Levenshtein.ratio(child_seq, parent_seq)
        for child_seq, parent_seq in zip(
        data["seq"],
        parent_raw,
    )
    ),
    dtype=np.float64,
    count=len(data),
)

parent_child_dtype = [
    ("count", np.int64),
    ("parent", object),
    ("child", object),
    ("similarity", np.float64),
]

PARENT_CHILD_INFO = np.empty(
    len(data),
    dtype=parent_child_dtype,
)
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# [(child's copy count, parent, child, parent-child similarity score)]

PARENT_CHILD_INFO["count"] = data["count"]

PARENT_CHILD_INFO["parent"] = np.fromiter(
    (parent.strip() for parent in parent_raw),
    dtype=object,
    count=len(parent_raw),
)

PARENT_CHILD_INFO["child"] = data["seq"]

PARENT_CHILD_INFO["similarity"] = similarities



#############################################


parent_pF_length = pF_length
parent_insert_length = payload_length
parent_pR_length = pR_length


aligner = Align.PairwiseAligner()
aligner.mode = "global"


child_info_dtype = [
    ("parent", object),
    ("count", np.int64),
    ("pF", object),
    ("insert", object),
    ("pR", object),
    ("similarity", np.float64),
]

CHILD_INFO = np.empty(
    len(PARENT_CHILD_INFO),
    dtype=child_info_dtype,
)

# Calculating child's pF, pR and insert given a parent --> for illumina SE and PE sequencing

for index, item in enumerate(PARENT_CHILD_INFO):
    child_copy = int(item["count"])
    parent_seq = item["parent"]
    child_seq = item["child"]
    parent_child_score = float(item["similarity"])

    if child_seq == parent_seq:
        pos1 = parent_pF_length
        pos2 = pos1 + parent_insert_length

        child_pF = child_seq[:pos1]
        child_insert = child_seq[pos1:pos2]
        child_pR = child_seq[pos2:]
    else:
        alignments = aligner.align(
            child_seq,
            parent_seq,
        )

        best_match = alignments[0]

        child = list(best_match[0])
        parent = list(best_match[1])


        child_pF, pos1 = get_clean_segment(
            parent,
            child,
            0,
            parent_pF_length,
        )

        child_insert, pos2 = get_clean_segment(
            parent,
            child,
            pos1,
            parent_insert_length,
        )

        child_pR = "".join(
            child[pos2:]
        ).replace("-", "")

    CHILD_INFO[index] = (
        parent_seq,
        child_copy,
        child_pF,
        child_insert,
        child_pR,
        parent_child_score,
    )

###############################################################


# Quality Check on child: 1) Alignment score > 0.80
# (2) parent-child length difference with +- 10 [optional] ---> consider lengths of mutated seqs, variants and chimeras,
# they might be filtered, so removing them not so good unless they are very short as compared to parent!!


alignment_THRESHOLD = 0.80
len_VALVE = 10

child_insert_lengths = np.fromiter(
    (
        len(child_insert)
        for child_insert in CHILD_INFO["insert"]
    ),
    dtype=np.int64,
    count=len(CHILD_INFO),
)

quality_mask = (
        (CHILD_INFO["similarity"] > alignment_THRESHOLD)
        & (
                np.abs(payload_length - child_insert_lengths)
                < len_VALVE
        )
)

LIST_01 = CHILD_INFO[quality_mask].copy()


###########################################


# Sequence input pool --> after sampling
list_02_dtype = [
    ("parent", object),
    ("count", np.int64),
    ("pF", object),
    ("insert", object),
    ("pR", object),
    ("length", np.int64),
]

LIST_02 = np.empty(
    len(LIST_01),
    dtype=list_02_dtype,
)

LIST_02["parent"] = LIST_01["parent"]
LIST_02["pF"] = LIST_01["pF"]
LIST_02["insert"] = LIST_01["insert"]
LIST_02["pR"] = LIST_01["pR"]


LIST_02["count"] = np.fromiter(
    (
        int(SAMPLING_FRAC * int(child_count))
        for child_count in LIST_01["count"]
    ),
    dtype=np.int64,
    count=len(LIST_01),
)


LIST_02["length"] = np.fromiter(
    (
        len(child_pF + child_insert + child_pR)
        for child_pF, child_insert, child_pR in zip(
        LIST_01["pF"],
        LIST_01["insert"],
        LIST_01["pR"],
    )
    ),
    dtype=np.int64,
    count=len(LIST_01),
)


##########################################################


# LIST_02 currently stores:
# (parent, sampled_count, child_pF, child_insert, child_pR, length)

# Use only positive sampled counts to construct multinomial weights.
positive_mask = LIST_02["count"] > 0
sequencing_pool = LIST_02[positive_mask]

weights = sequencing_pool["count"].astype(float)
weights /= weights.sum()

read_allocations = np.random.multinomial(
    TARGET_READS,
    weights,
)# allocating how many reads each template gets


# generating sequencing reads for the chosen method
READ_SE = []
READ_PE = []

for item, n_reads in zip(sequencing_pool, read_allocations):
    pF = item["pF"]
    insert = item["insert"]
    pR = item["pR"]

    insert_length = len(insert)
    n_reads = int(n_reads)

    if mode == "1":
        if method == "1":  # Illumina single-end
            if READ_LENGTH > insert_length:
                read1 = (insert + pR)[:READ_LENGTH]
            elif READ_LENGTH == insert_length:
                read1 = insert
            else:
                read1 = insert[:READ_LENGTH]

            # Equivalent to appending the same read n_reads times.
            READ_SE.extend([read1] * n_reads)

        else:  # Illumina paired-end
            if READ_LENGTH > insert_length:
                read1 = (insert + pR)[:READ_LENGTH]
                read2 = str(
                    Seq(
                        (pF + insert)[-READ_LENGTH:]
                    ).reverse_complement()
                )

            elif READ_LENGTH == insert_length:
                read1 = insert
                read2 = str(
                    Seq(insert).reverse_complement()
                )

            else:
                read1 = insert[:READ_LENGTH]
                read2 = str(
                    Seq(
                        insert[-READ_LENGTH:]
                    ).reverse_complement()
                )

            READ_SE.extend([read1] * n_reads)
            READ_PE.extend([read2] * n_reads)



#################################################

# writing and rereading Read_SE.txt / Read_PE.txt.

if method == "1":

    read_sequences = (
        read.strip()
        for read in READ_SE
    )

elif method == "2":

    read_sequences = (
        read.strip()
        for read_pair in zip(READ_SE, READ_PE)
        for read in read_pair
    )

else:
    raise ValueError("Unsupported method")


print(f"Length SE: {len(READ_SE)}")
print(f"Length PE: {len(READ_PE)}")


#################################################


# for fixed quality score
"""
def phred_char(q=30):
    return chr(q + 33)  # Q30 -> '?'

def qual_string(s, q=30):
    return phred_char(q) * len(s)
"""

# for start to end varying quality score, initially higher at q_start then goes down lower till q_end
"""def quality_score(read_length, q_start=35, q_end=28):
    if read_length == 1:
        return chr(q_start + 33)
    qs = []
    for i in range(read_length):
        q = round(q_start + (q_end - q_start) * i / (read_length - 1))
        qs.append(chr(q + 33))
    return "".join(qs)"""


def q_to_char(q):
    q = max(2, min(40, int(q)))
    return chr(q + 33)


def quality_score(read_length, read_type="R1"):
    r = random.random()
    if r < 0.85:
        base_q = random.randint(30, 38)  # 85% of the base reads
    elif r < 0.97:
        base_q = random.randint(20, 29)  # 12% of the base reads
    else:
        base_q = random.randint(8, 19)  # 3% of the base reads

    if read_type == "R2":
        base_q -= random.randint(0,
                                 3)  # R2 base reads quality are slightly less than R1, so decrease base_quality a little

    quals = []
    for i in range(read_length):
        q = base_q + random.randint(-2, 2)  # a small jitter [OPTIONAL]

        if i > 0.8 * read_length:  # slight tail decay in last 20% of read
            q -= random.randint(0, 4)

        quals.append(q_to_char(q))
    return "".join(quals)



# Making sequencing object + Applying sequencing errors
seq_objs = [
    Error_simulation(
        seq,
        "sequencing",
        attribute=mutation_attributes[mode],
        error_rate=err_rates[method],
    )
    for seq in read_sequences
]

#################################################

VALVE_MAP = {
    "0": 0,
    "1": 1,
    "2": 2,
    "3": 4,
}

if args.c is not None:
    VALVE = args.c
    print("Using custom VALVE")
else:
    try:
        VALVE = VALVE_MAP[args.mut]
    except KeyError:
        raise ValueError(
            "Invalid mutation VALVE [try --help]"
        ) from None

print(f"Running with VALVE: {VALVE}")


for _ in range(VALVE):
    for seq_obj in seq_objs:
        seq_obj.run_mutations()


MUTATED_TEXT = [
    seq_obj.seq
    for seq_obj in seq_objs
]

MUTATED_TEXT_FINAL = [
    "".join(read).replace(" ", "")
    for read in MUTATED_TEXT
]

if len(MUTATED_TEXT_FINAL) == 0:
    print(
        f"length MUTATED LIST FINAL: "
        f"{len(MUTATED_TEXT_FINAL)}"
    )
    raise ValueError("Sequencing Failure!!")

###################################################

if mode == "1":

    # Illumina single-end
    if method == "1":
        filepath = fr"{SEQ_DIR}/sequencing_R1.fastq"

        with open(filepath, "w") as f:
            for read_id, read in enumerate(
                MUTATED_TEXT_FINAL,
                start=1,
            ):
                qual1 = quality_score(
                    len(read),
                    "R1",
                )

                f.write(
                    f"@read{read_id}\n"
                    f"{read}\n"
                    f"+\n"
                    f"{qual1}\n"
                )

    # Illumina paired-end
    elif method == "2":
        r1_filepath = fr"{SEQ_DIR}/sequencing_R1.fastq"
        r2_filepath = fr"{SEQ_DIR}/sequencing_R2.fastq"

        with (
            open(r1_filepath, "w") as f1,
            open(r2_filepath, "w") as f2,
        ):
            read_id = 1

            for i in range(
                0,
                len(MUTATED_TEXT_FINAL),
                2,
            ):
                read1 = MUTATED_TEXT_FINAL[i]
                read2 = MUTATED_TEXT_FINAL[i + 1]

                qual1 = quality_score(
                    len(read1),
                    "R1",
                )

                qual2 = quality_score(
                    len(read2),
                    "R2",
                )

                f1.write(
                    f"@read{read_id}/1\n"
                    f"{read1}\n"
                    f"+\n"
                    f"{qual1}\n"
                )

                f2.write(
                    f"@read{read_id}/2\n"
                    f"{read2}\n"
                    f"+\n"
                    f"{qual2}\n"
                )

                read_id += 1

    else:
        raise ValueError("Unsupported method")
#################################################


if __name__ == "__main__":
    print("Sequencing.py run completed")
