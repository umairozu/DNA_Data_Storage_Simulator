import argparse
import os
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

def copy_distribution(avg_copy_oligo=100e6, k=4):
    # k is the divergence parameter
    min_copy_count = 1  # if 0, then we are allowing dropouts
    #rng = np.random.default_rng(seed)
    p_mol = 6.022e11  # concentration

    p = k / (k + avg_copy_oligo)

    oligo_count = np.random.negative_binomial(k, p)
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
_p.add_argument("--out_file", default="synthesis_output",
                help="Output file name [default is 'synthesis_output'] ")

args = _p.parse_args()
in_file = fr'dna-fountain/{args.in_file}.tar.gz.dna_order.txt'
out_file = fr'{args.out_file}.txt'

VALVE_HIGH = 20
VALVE_MED = 15
VALVE_LOW = 10
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
    initial_lines = [line.strip() for line in f if line.strip()]
    seq_objs = [Error_simulation(seq, "synthesis", attribute=mutation_attributes["3"],
                                 error_rate=err_rates["3"])
                for seq in initial_lines
                ]

MUTATED_TEXT = []

if VALVE == 0:
    for sE in seq_objs:
        MUTATED_TEXT.append(sE.seq)

count = 1
while count < VALVE:
    MUTATED_TEXT.clear()
    for sE in seq_objs:
        # sE.reset_visited() # See important Notice for info
        sE.run_mutations()
        MUTATED_TEXT.append(sE.seq)
    count += 1

with open(fr'{SYNTHESIS_DIR}/synthesis_file_0.txt', "w") as f:
    f.write("\n".join(MUTATED_TEXT) + "\n")

NEW = []
with open(fr'{SYNTHESIS_DIR}/synthesis_file_0.txt') as f:
    for line in f:
        clean_lines = line.split()[0]
        NEW.append(clean_lines)

with open(fr'{SYNTHESIS_DIR}/synthesis_file_1.txt', "w") as f:
    f.write("\n".join(NEW) + "\n")

    # file read for copy:
with open(fr'{SYNTHESIS_DIR}/synthesis_file_1.txt') as f:
    lines = f.readlines()

copies = []
for line in lines:
    oligo_copies = copy_distribution()
    copies.append(oligo_copies)

with open(fr'{SYNTHESIS_DIR}/{out_file}', "w") as f:
    for copy, line in zip(copies, lines):
        f.write(f"{copy},{line}")

if __name__ == "__main__":

    os.remove(fr'{SYNTHESIS_DIR}/synthesis_file_0.txt')
    os.remove(fr'{SYNTHESIS_DIR}/synthesis_file_1.txt')

    print("Synthesis.py run completed")


