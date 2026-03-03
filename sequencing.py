

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

    with open(fr'{os.getcwd()}\dna-fountain\turkish_anthem.tar.gz.dna_order.txt') as f:
        initial_lines = [line.strip() for line in f if line.strip()]
        seq_objs = [Error_simulation(seq, "sequencing", attribute = mutation_attributes["1"],
                                     error_rate = err_rates["3"])
                    for seq in initial_lines
                    ]

    count = 1
    while count < 300:
        MUTATED_TEXT.clear()
        for sE in seq_objs:
            #sE.reset_visited() # See important Notice for info
            sE.run_mutations()
            MUTATED_TEXT.append(sE.seq)

        with open(fr'{os.getcwd()}\dna-fountain\turkish_anthem.tar.gz.dna_order.txt', "w") as f:
            f.write("\n".join(MUTATED_TEXT) + "\n")
        count += 1
