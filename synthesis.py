
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



MUTATED_TEXT = []

if __name__ == "__main__":

    with open(fr'{os.getcwd()}\dna-fountain\turkish_anthem.tar.gz.dna_order.txt') as f:
        initial_lines = [line.strip() for line in f if line.strip()]
        seq_objs = [Error_simulation(seq, "synthesis", attribute = mutation_attributes["3"],
                                     error_rate = err_rates["3"])
                    for seq in initial_lines
                    ]

    count = 1
    while count < 15:
        MUTATED_TEXT.clear()
        for sE in seq_objs:
            #sE.reset_visited() # See important Notice for info
            sE.run_mutations()
            MUTATED_TEXT.append(sE.seq)

        with open(fr'{os.getcwd()}\dna-fountain\synthesis_file_0.txt', "w") as f:
            f.write("\n".join(MUTATED_TEXT) + "\n")
        count += 1



    NEW = []
    with open(fr'{os.getcwd()}\dna-fountain\synthesis_file_0.txt') as f:
        for line in f:
            clean_lines = line.split()[0]
            NEW.append(clean_lines)

    with open(fr'{os.getcwd()}\dna-fountain\synthesis_file_1.txt', "w") as f:
        f.write("\n".join(NEW) + "\n")


    #file read for copy:
    with open(fr'{os.getcwd()}\dna-fountain\synthesis_file_1.txt') as f:
        lines = f.readlines()


    copies = []
    for line in lines:
        oligo_copies = copy_distribution()
        copies.append(oligo_copies)


    with open(fr'{os.getcwd()}\dna-fountain\synthesis_file_2.txt', "w") as f:
        for copy, line in zip(copies, lines):
            f.write(f"{copy},{line}")


