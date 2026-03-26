

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
import os
from Error_module import Error_simulation

BASE_DIR = fr'{os.getcwd()}\dna-fountain'

err_rates = {"1": {"raw_rate": 0.0021, "substitution": 0.81, "deletion": 0.0024, "insertion": 0.0013},
             "2": {"raw_rate": 0.0032, "substitution": 0.79, "deletion": 0.0018, "insertion": 0.0011},
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


MUTATED_TEXT = []

if __name__ == "__main__":

    mode = input("Enter Sequencing Mode:  --help [1. Illumina, 2. PacBio, 3. Nanopore ]\n")
    #assert mode in ["1","2","3"]
    if mode == "1":
        method = input("Enter method for illumina sequencing:  --help [1. Single-end, 2. Paired-end]\n")
        #assert method in ["1", "2"]
    elif mode == "2":
        method = input("Enter method for PacBio sequencing:  --help [3. CCS, 4. Subread]\n")
        #assert method in ["3", "4"]
    elif mode == "3":
        method = input("Enter method for Nanopore sequencing:  --help [5. 1D, 6. 2D]\n")
        #assert method in ["5", "6"]
    else:
        mode = "1"
        method = "1"
        print("Default sequencing method chosen [illumina single-end sequencing]")



    with open(fr'{BASE_DIR}\pcr_final.txt') as f:
        rows = [line.strip().split(",") for line in f if line.strip()]
        copy_count, lines, length = zip(*rows)
        seq_objs = [Error_simulation(seq, "sequencing", attribute = mutation_attributes[mode],
                                     error_rate = err_rates[method])
                    for seq in lines
                    ]

# sampling for sequencing
while True:
    user_input = input("Enter sequencing sampling fraction (0-100%): \n")

    try:
        value = float(user_input)
        if 0 <= value <= 100:
            sampling_frac = value / 100
            break
        else:
            print("Please enter a number between 0 and 100.")
    except ValueError:
        print("Invalid input. Please enter a numerical value (e.g., 5.5).")

# sequencing depth input
while True:
    user_input_input = input("Enter a number for sequencing depth (0-100%): \n")

    try:
        value = float(user_input)
        if 0 <= value <= 100:
            sequencing_depth = value / 100
            break
        else:
            print("Please enter a number between 0 and 100.")
    except ValueError:
        print("Invalid input. Please enter a numerical value (e.g., 5.5).")



with open(fr'{BASE_DIR}\pcr_final.txt') as f:
    next(f)
    data = []
    for line in f:
        if "," in line.strip():
            count, seq, length = line.split(",")
            data.append({ 'count': int(count.strip()), 'seq': seq.strip(), 'length': length.strip() })



    """    count = 1
    while count < 100:
        MUTATED_TEXT.clear()
        for sE in seq_objs:
            #sE.reset_visited() # See important Notice for info
            sE.run_mutations()
            MUTATED_TEXT.append(sE.seq)

        with open(fr'{os.getcwd()}\dna-fountain\turkish_anthem.tar.gz.dna_order.txt', "w") as f:
            f.write("\n".join(MUTATED_TEXT) + "\n")
        count += 1"""


    print("Sequencing.py run")