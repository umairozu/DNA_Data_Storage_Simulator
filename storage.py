import argparse
import math
import os
import re
import numpy as np
from Error_module import Error_simulation
from Arrhenius_decay import Arrhenius_decay

SEC_PER_WEEK = 7 * 24 * 3600
SEC_PER_YEAR = 365 * 24 * 3600

xlsx = "RawData.xlsx"
sim = Arrhenius_decay.from_xlsx(xlsx)

BASE_DIR = fr'{os.getcwd()}'
os.makedirs(fr'{BASE_DIR}/storage', exist_ok=True)
SYNTHESIS_DIR = fr'{os.getcwd()}/synthesis'
STORAGE_DIR = fr'{os.getcwd()}/storage'

def survival_distribution(copy_count, temp_C, encapsulated, week):
    remaining_frac = Arrhenius_decay.remaining_dna_frac(sim, temp_C, encapsulated, week)
    remaining_copies = copy_count * remaining_frac

    return remaining_copies

# Per base depurination calculation
def depurination(pH,temp_C, encapsulated):

    T = 273.15 + float(temp_C)
    if pH < 2.5:
        lgk = 14.6 - 0.707 * pH - (5.63e3 / T)
    else:
        lgk = 16.5 - 0.982 * pH - (5.85e3 / T)
    if encapsulated == "False":
        k = 10 ** lgk
    else:
        k = 9e-10 # still expecting some mutations (deletions) when protected

    return k

# ln N(t) = -kt + ln N0
# N(t) = e^(-kt + ln N0)
# N(t) = N0 * e^(-kt)

# Remaining intact fraction = N(t) / N0 = e^(-kt)
# Fraction lost to depurination = 1 - Remaining intact fraction
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# per base (A/G) mutation probability
# this will be the raw rate for mutation

def fraction_lost_depurination(pH, temp_C, encapsulated, week):
    t = SEC_PER_WEEK * week
    k = depurination(pH, temp_C,encapsulated)
    return 1 - math.exp(-k * t)


MUTATED_TEXT = []
MUTATED_COPY_COUNT = []

"""input_temp_C = float(input("Enter Temperature in Celsius: \n"))

input_pH = float(input("Enter pH value (specify good range for user): \n"))
while input_pH < 0.0:
    input_pH = float(input("Enter pH value (specify good range for user): \n"))

input_weeks = float(input("Enter number of weeks for storage: \n"))
while input_weeks < 0.0:
    input_weeks = float(input("Enter number of weeks for storage: \n"))

input_encapsulation = input("Enter if Encapsulated or not (True | False): \n").capitalize()
while input_encapsulation != "True" and input_encapsulation != "False":
    print("Invalid input! Please type True or False.")
    input_encapsulation = input("Enter if Encapsulated or not (True | False): \n").capitalize()"""


def valve_range(value):
    v = int(value)
    if not (0 <= v <= 14):
        raise argparse.ArgumentTypeError("pH must be between 0 and 14.")
    return v

_p = argparse.ArgumentParser(description="storage.py run")
_p.add_argument("--temp",required = True , help = "Temperature in Celsius")
_p.add_argument("--ph", default = "7", type = valve_range, help = "Enter pH value [Keep at 7 neutral for good results]" )
_p.add_argument("--week", type = lambda x: float(x) if float(x) > 0.0 else _p.error("Week must be > 0.0") ,
                                        required = True, help = "number of weeks for storage")
_p.add_argument("--encap", default = "0", type = lambda x: x if x in ["0","1"] else _p.error("True = 1, False = 0") ,
                                            help = "Specify if the DNA is zif protected, True = 1, False = 0 ")
_p.add_argument("--mut", default = "4" , help = "Mutation intensity (0-3)")
_p.add_argument("--c", type = lambda x: int(x) if  0 < int(x) <= 30 else _p.error("keep 0 < custom VALVE <= 30"),
                                        help = "Optional custom VALVE (0-30), it is basically a mutation knob" )
_p.add_argument("--in_file", required = True, help="Input file name (e.g  synthesis_in)")
_p.add_argument("--out_file", default = "storage_output.txt", help = "Output file name [default is 'storage_output'] ")

args = _p.parse_args()

input_temp_C = float(args.temp)
input_pH = float(args.ph)
input_weeks = float(args.week)

if args.encap == "0":
    input_encapsulation = "False"
else:
    input_encapsulation = "True"

in_file = fr'{args.in_file}'
out_file = fr'{STORAGE_DIR}/{args.out_file}'


raw_rate = fraction_lost_depurination(input_pH, input_temp_C, input_encapsulation, input_weeks)

print("\n--------\nSummary\n--------")
print(f"input Temp: {input_temp_C}")
print(f"input pH: {input_pH}")
print(f"input weeks: {input_weeks}")
print(f"Encapsulated: {input_encapsulation}")

print("k =", depurination(input_pH, input_temp_C, input_encapsulation), "raw_rate =", raw_rate)
print("\n")

err_rates = {
    "1": {"deletion": 0.33, "insertion": 0.33, "substitution": 0.34, "raw_rate": 0.000000021},  # D. melanogaster
    "2": {"deletion": 0.13, "insertion": 0.13, "substitution": 0.74, "raw_rate": 0.000000079},  # S cerevisiae
    "3": {"deletion": 0.08, "insertion": 0.08, "substitution": 0.84, "raw_rate": 0.000000317},  # E. coli
    "4": {"deletion": 0.06, "insertion": 0.06, "substitution": 0.88, "raw_rate": 0.000000000069},  # H sapiens
    "5": {"deletion": 0.025, "insertion": 0.025, "substitution": 0.95, "raw_rate": 0.0000000044},  # M musculus
    "6": {"deletion": 0.3333, "insertion": 0.3333, "substitution": 0.33340000000000003, "raw_rate": 0.0}, # Depurination at pH 7 and 293.15K
    "7": {"deletion": 1, "insertion": 0, "substitution": 0, "raw_rate": 0.0001231},  # Depurination at pH 8 and 253.15K
    "8": {"deletion": 1, "insertion": 0, "substitution": 0, "raw_rate": 1e-08}, # Erasure Channel with an error probability of 0.5 percent
    "9": {"deletion": 1, "insertion": 0, "substitution": 0, "raw_rate": 0.005},  # Depurination at pH 7 and 193.15K
    "10": {"deletion": 1, "insertion": 0, "substitution": 0, "raw_rate": 5.736e-15}, # White Gaussian Noise with an error probability of 0.5 percent
    "11": {"deletion": 0, "insertion": 0, "substitution": 1, "raw_rate": 0.005},  # Depurination at pH 8 and 193.15K
    "12": {"deletion": 1, "insertion": 0, "substitution": 0, "raw_rate": 5.98e-16},  # Depurination at pH 8 and 293.15K
    "13": {"deletion": 1, "insertion": 0, "substitution": 0, "raw_rate": 1.283e-05},  # Depurination at pH 7 and 253.15K
    "14": {"deletion": 1, "insertion": 0, "substitution": 0, "raw_rate": 9e-08},  # None
    "15": {"deletion": 1, "insertion": 0, "substitution": 0, "raw_rate": raw_rate}  # General
}

mutation_attributes = {
    # D. melanogaster
    "1": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                       "position": {"homopolymer": 0, "random": 1}},
          "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                        "position": {"homopolymer": 0, "random": 1}},
          "substitution": {"pattern": {}}},

    # S cerevisiae
    "2": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                       "position": {"homopolymer": 0, "random": 1}},
          "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                        "position": {"homopolymer": 0, "random": 1}},
          "substitution": {"pattern": {}}},

    # E. coli
    "3": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                       "position": {"homopolymer": 0, "random": 1}},
          "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                        "position": {"homopolymer": 0, "random": 1}},
          "substitution": {"pattern": {}}},

    # H sapiens
    "4": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                       "position": {"homopolymer": 0, "random": 1}},
          "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                        "position": {"homopolymer": 0, "random": 1}},
          "substitution": {"pattern": {}}},

    # M musculus
    "5": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                       "position": {"homopolymer": 0, "random": 1}},
          "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                        "position": {"homopolymer": 0, "random": 1}},
          "substitution": {"pattern": {}}},

    # Depurination at pH 7 and 293.15K
    "6": {"deletion": {"pattern": {"A": 0.5, "C": 0.0, "G": 0.5, "T": 0.0},
                       "position": {"homopolymer": 0.0, "random": 1.0}},
          "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                        "position": {"homopolymer": 0.5, "random": 0.5}},
          "substitution": {"pattern": {}}},

    # Depurination at pH 8 and 253.15K
    "7": {"deletion": {"pattern": {"A": 0.5, "C": 0.0, "G": 0.5, "T": 0.0},
                       "position": {"homopolymer": 0.0, "random": 1.0}},
          "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                        "position": {"homopolymer": 0.5, "random": 0.5}},
          "substitution": {"pattern": {}}},

    # Erasure Channel with an error probability of 0.5 percent
    "8": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                       "position": {"homopolymer": 0.0, "random": 1.0}},
          "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                        "position": {"homopolymer": 0.5, "random": 0.5}},
          "substitution": {"pattern": {}}},

    # Depurination at pH 7 and 193.15K
    "9": {"deletion": {"pattern": {"A": 0.5, "C": 0.0, "G": 0.5, "T": 0.0},
                       "position": {"homopolymer": 0.0, "random": 1.0}},
          "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                        "position": {"homopolymer": 0.5, "random": 0.5}},
          "substitution": {"pattern": {}}},

    # White Gaussian Noise with an error probability of 0.5 percent
    "10": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                        "position": {"homopolymer": 0.5, "random": 0.5}},
           "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                         "position": {"homopolymer": 0.5, "random": 0.5}},
           "substitution": {"pattern": {"A": {"T": 0.3333, "G": 0.33340000000000003, "C": 0.3333},
                                        "T": {"A": 0.3333, "G": 0.33340000000000003, "C": 0.3333},
                                        "C": {"G": 0.3333, "T": 0.33340000000000003, "A": 0.3333},
                                        "G": {"A": 0.3333, "T": 0.33340000000000003, "C": 0.3333}}}},

    # Depurination at pH 8 and 193.15K
    "11": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                        "position": {"homopolymer": 0.5, "random": 0.5}},
           "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                         "position": {"homopolymer": 0.5, "random": 0.5}},
           "substitution": {"pattern": {}}},

    # Depurination at pH 8 and 293.15K
    "12": {"deletion": {"pattern": {"A": 0.5, "C": 0.0, "G": 0.5, "T": 0.0},
                        "position": {"homopolymer": 0.0, "random": 1.0}},
           "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                         "position": {"homopolymer": 0.5, "random": 0.5}},
           "substitution": {"pattern": {}}},

    # Depurination at pH 7 and 253.15K
    "13": {"deletion": {"pattern": {"A": 0.5, "C": 0.0, "G": 0.5, "T": 0.0},
                        "position": {"homopolymer": 0.0, "random": 1.0}},
           "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                         "position": {"homopolymer": 0.5, "random": 0.5}},
           "substitution": {"pattern": {}}},

    # None
    "14": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                        "position": {"homopolymer": 0.5, "random": 0.5}},
           "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                         "position": {"homopolymer": 0.5, "random": 0.5}},
           "substitution": {"pattern": {}}},

    # Depurination General
    "15": {"deletion": {"pattern": {"A": 0.50, "C": 0.0, "G": 0.50, "T": 0.0},
                        "position": {"homopolymer": 0.0, "random": 1.0}},
           "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                         "position": {"homopolymer": 0.0, "random": 1.0}},
           "substitution": {"pattern": {}}}

}

with open(fr'{in_file}') as f:
    initial_lines = [line.split(",")[1].strip() for line in f if line.strip()]
    f.seek(0)
    initial_copies = [line.split(",")[0].strip() for line in f if line.strip()]
    seq_objs = [Error_simulation(seq, "storage", attribute=mutation_attributes["15"],
                                 error_rate=err_rates["15"])
                for seq in initial_lines
                ]
VALVE_HIGH = 6
VALVE_MED = 4
VALVE_LOW = 2
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
while count < VALVE:  # adjust the exit condition
    MUTATED_TEXT.clear()
    for sE in seq_objs:
        # sE.reset_visited() # See important Notice for info
        sE.run_mutations()
        MUTATED_TEXT.append(sE.seq)
    count += 1

for copy in initial_copies:
    MUTATED_COPY_COUNT.append(survival_distribution(int(copy), input_temp_C, input_encapsulation, input_weeks))

# Replacing copy count and mutated lines in file after degradation
with open(fr'{STORAGE_DIR}/storage_file_0.txt', "w") as f:
    for copy, line in zip(MUTATED_COPY_COUNT, MUTATED_TEXT):
        f.write(f"{copy:.0f},{line}\n")

FRAGMENT_COUNT = 0

# Removing spaces in a line and making it shorter instead of throwing the line away
with open(fr'{STORAGE_DIR}/storage_file_0.txt') as f:
    MUTATED_TEXT = [line.split(",")[1].replace(" ", "").strip() for line in f if line.strip()]

with open(out_file, "w") as f:
    f.write("count, sequence, length\n")
    for copy, line in zip(MUTATED_COPY_COUNT, MUTATED_TEXT):
        f.write(f"{copy:.0f},{line},{len(line)}\n")


if __name__ == "__main__":

    os.remove(fr'{STORAGE_DIR}/storage_file_0.txt')
    print("Storage.py run completed")