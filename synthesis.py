import argparse
import json
import os
import random
import numpy as np
from Error_module import Error_simulation


err_rates = {
             "1": {"raw_rate": 0.0001, "substitution": 0.15, "deletion": 0.7, "insertion": 0.15}, #MutS
             "2": {"raw_rate": 0.000125, "substitution": 0.15, "deletion": 0.7, "insertion": 0.15}, #Consensus Shuffle
             "3": {"raw_rate": 0.00011, "substitution": 0.2, "deletion": 0.6, "insertion": 0.2}, #NGS-based error correction
             "4": {"raw_rate": 0.0017, "substitution": 0.2, "deletion": 0.6, "insertion": 0.2}, #High-temperature ligation/hybridization based error correction
             "5": {"raw_rate": 0.00125, "substitution": 0.2, "deletion": 0.6, "insertion": 0.2}, #ErrASE
             "6": {"raw_rate": 0.00033, "substitution": 0.2, "deletion": 0.6, "insertion": 0.2}, #Nuclease-based error correction
             "7": {"raw_rate": 0.000025, "substitution": 0.2, "deletion": 0.6, "insertion": 0.2}, #ErrASE
             "8": {"raw_rate": 0.0004, "substitution": 0.2, "deletion": 0.6, "insertion": 0.2}, #Oligo Hybridization based error correction
             "9": {"raw_rate": 0, "substitution": 0.3, "deletion": 0.3, "insertion": 0.3} #None
             }

mutation_attributes = {
                        #MutS
                        "1": {"deletion": {"position": {"homopolymer": 0.0, "random": 1}, "pattern": {"G": 0.2, "C": 0.2, "A": 0.4, "T": 0.2}},"insertion": {"position": {"homopolymer": 0, "random": 1},"pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},"substitution": {}},

                        #Consensus Shuffle
                        "2": {"deletion": {"position": {"homopolymer": 0.0, "random": 1}, "pattern": {"G": 0.2, "C": 0.2, "A": 0.4, "T": 0.2}},"insertion": {"position": {"homopolymer": 0, "random": 1},"pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},"substitution": {}},

                        #NGS-based error correction
                        "3": {"deletion": {"position": {"homopolymer": 0.0, "random": 1}, "pattern": {"G": 0.2, "C": 0.2, "A": 0.4, "T": 0.2}},"insertion": {"position": {"homopolymer": 0, "random": 1},"pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},"substitution": {}},

                        #High-temperature ligation/hybridization based error correction
                        "4": {"deletion": {"position": {"homopolymer": 0.0, "random": 1}, "pattern": {"G": 0.2, "C": 0.2, "A": 0.4, "T": 0.2}},"insertion": {"position": {"homopolymer": 0, "random": 1},"pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},"substitution": {"pattern": {"AAAA": {"ACGT":1.0}, "ACCC": {"ACGC":1.0}}}},

                        #ErrASE
                        "5": {"deletion": {"position": {"homopolymer": 0.0, "random": 1}, "pattern": {"G": 0.2, "C": 0.2, "A": 0.4, "T": 0.2}},"insertion": {"position": {"homopolymer": 0, "random": 1},"pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},"substitution": {}},

                        #Nuclease-based error correction
                        "6": {"deletion": {"position": {"homopolymer": 0.0, "random": 1}, "pattern": {"G": 0.2, "C": 0.2, "A": 0.4, "T": 0.2}},"insertion": {"position": {"homopolymer": 0, "random": 1},"pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},"substitution": {}},

                        #ErrASE
                        "7": {"deletion": {"pattern": {"A": 0.4, "C": 0.2, "G": 0.2, "T": 0.2}, "position": {"homopolymer": 0, "random": 1}}, "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}, "position": {"homopolymer": 0, "random": 1}}, "substitution": {"pattern": {}}},

                        #Oligo Hybridization based error correction
                        "8": {"deletion": {"pattern": {"A": 0.4, "C": 0.2, "G": 0.2, "T": 0.2}, "position": {"homopolymer": 0, "random": 1}}, "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}, "position": {"homopolymer": 0, "random": 1}}, "substitution": {"pattern": {}}},

                        #None
                        "9": {"deletion": {"position": {"homopolymer": 0.0, "random": 1}, "pattern": {"G": 0.2, "C": 0.2, "A": 0.4, "T": 0.2}},"insertion": {"position": {"homopolymer": 0, "random": 1},"pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},"substitution": {}}

                        }


BASE_DIR = fr'{os.getcwd()}'
os.makedirs(fr'{os.getcwd()}/synthesis', exist_ok=True)
SYNTHESIS_DIR = fr'{os.getcwd()}/synthesis'



metadata_path = fr'{BASE_DIR}/metadata.json'
with open(metadata_path, "r") as f:
    metadata = json.load(f)

orig_length = metadata["orig_length"]

def efficiency(L = orig_length, coup_eff = 0.999):
    return coup_eff ** (L - 1)


def copy_distribution(seqs, avg_copy_oligo=1e8, k=35):
    # k is the divergence parameter
    min_copy_count = 1  # if 0, then we are allowing dropouts
    #rng = np.random.default_rng(seed)
    p_mol = 6.022e11  # concentration

    p = k / (k + avg_copy_oligo)

    oligo_count = np.max(np.random.negative_binomial(k, p, size = len(seqs)))
    if oligo_count == 0:
        print(oligo_count)
    return oligo_count

def valve_range(value):
    v = int(value)
    if not (0 <= v <= 100):
        raise argparse.ArgumentTypeError("custom VALVE must be between 0 and 100.")
    return v


_p = argparse.ArgumentParser(description="synthesis.py run")
_p.add_argument("--mut", default = "2" , help = "Mutation intensity (0-3)")
_p.add_argument("--c", type = valve_range, help = "Optional custom VALVE (0-100), it is basically a mutation knob" )
_p.add_argument("--in_file", required=True, help="Input file name (e.g  test_file)")
_p.add_argument("--out_file", default="synthesis_output.txt",
                help="Output file name [default is 'synthesis_output'] ")

args = _p.parse_args()
in_file = fr'{args.in_file}'
out_file = fr'{args.out_file}'

VALVE_HIGH = 20
VALVE_MED = 10
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


with open(fr'{BASE_DIR}/{in_file}') as f:
    rows = [line.strip().split(",") for line in f if line.strip()]
    columns = list(zip(*rows))

    pid = list(columns[0])
    sequences = list(columns[1])

    np_pid_LINES = np.array(pid)
    np_LINES = np.array(sequences)

counts_list = [copy_distribution(line) for line in np_LINES]

np_COUNT = np.array(counts_list)

eff = efficiency()

np_COUNT2 = (np_COUNT * eff).astype(int)

with open(fr'{SYNTHESIS_DIR}/file1.txt', "w+") as f:
    for pid, count, line in zip(np_pid_LINES, np_COUNT2, np_LINES):
        f.write(f"{pid},{count},{line}\n")
    f.seek(0)
    rows_01 = [line.strip().split(",") for line in f if line.strip()]

seq_objs = [Error_simulation(seq, "synthesis", attribute=mutation_attributes["3"], error_rate=err_rates["3"])
                for seq in np_LINES
            ]

MUTATED_TEXT = []

count = 0
while count < VALVE:
    MUTATED_TEXT.clear()
    for sE in seq_objs:
        # sE.reset_visited() # See important Notice for info
        sE.run_mutations()
        MUTATED_TEXT.append(sE.seq)
    count += 1

rows_02 = []

if VALVE != "0":
    np_COUNT3 = (np_COUNT * (1 - eff)).astype(int)
    np_MUTATED = np.array(MUTATED_TEXT)

    with open(fr'{SYNTHESIS_DIR}/file2.txt', "w+") as f:
        for parent_id, mutated_count, original_seq, mutated_seq in zip(
                np_pid_LINES,
                np_COUNT3,
                np_LINES,
                np_MUTATED
        ):
            # Writing only sequences that are actually mutated
            if mutated_seq != original_seq:
                cleaned_mutated_seq = mutated_seq.split()[0]

                f.write(
                    f"{parent_id},{mutated_count},{cleaned_mutated_seq}\n"
                )

        f.seek(0)
        rows_02 = [
            line.strip().split(",", maxsplit=2)
            for line in f
            if line.strip()
        ]

combined_rows = rows_01 + rows_02
random.shuffle(combined_rows)

with open(fr'{SYNTHESIS_DIR}/{out_file}', "w") as f:
    f.write(f"parent_id, count, sequence" +"\n")
    for row in combined_rows:
        f.write(f"{row[0]},{row[1]},{row[2]}\n")

if __name__ == "__main__":

    #os.remove(fr'{SYNTHESIS_DIR}/file1.txt')
    #os.remove(fr'{SYNTHESIS_DIR}/file2.txt')

    """
    print(sum(np_COUNT2))
    print("   +   ")
    print(sum(np_COUNT3))
    print("------------")
    print(sum(np_COUNT2 + np_COUNT3))
    print(sum(np_COUNT))

    """
    print("Synthesis.py run completed")


