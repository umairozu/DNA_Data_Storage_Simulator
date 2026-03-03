
import os
from Error_module import Error_simulation


err_rates = {
             "1": {"deletion": 0.33, "insertion": 0.33, "substitution": 0.34, "raw_rate": 0.000000021},
             "2": {"deletion": 0.13, "insertion": 0.13, "substitution": 0.74, "raw_rate": 0.000000079},
             "3": {"deletion": 0.08, "insertion": 0.08, "substitution": 0.84, "raw_rate": 0.000000317},
             "4": {"deletion": 0.06, "insertion": 0.06, "substitution": 0.88, "raw_rate": 0.000000000069},
             "5": {"deletion": 0.025, "insertion": 0.025, "substitution": 0.95, "raw_rate": 0.0000000044},
             "6": {"deletion": 0.3333, "insertion": 0.3333, "substitution": 0.33340000000000003, "raw_rate": 0.0},
             "7": {"deletion": 1, "insertion": 0, "substitution": 0, "raw_rate": 0.0001231},
             "8": {"deletion": 1, "insertion": 0, "substitution": 0, "raw_rate": 1e-08},
             "9": {"deletion": 1, "insertion": 0, "substitution": 0, "raw_rate": 0.005},
             "10":{"deletion": 1, "insertion": 0, "substitution": 0, "raw_rate": 5.736e-15},
             "11":{"deletion": 0.0, "insertion": 0.0, "substitution": 1.0, "raw_rate": 0.005},
             "12":{"deletion": 1, "insertion": 0, "substitution": 0, "raw_rate": 5.98e-16},
             "13":{"deletion": 1, "insertion": 0, "substitution": 0, "raw_rate": 1.283e-05},
             "14":{"deletion": 1, "insertion": 0, "substitution": 0, "raw_rate": 9e-08}
            }

mutation_attributes = {
                        #D melanogaster
                        "1": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                           "position": {"homopolymer": 0, "random": 1}},
                              "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                            "position": {"homopolymer": 0, "random": 1}},
                              "substitution": {"pattern": {}}},

                        #S cerevisiae
                        "2": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                           "position": {"homopolymer": 0, "random": 1}},
                              "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                            "position": {"homopolymer": 0, "random": 1}},
                              "substitution": {"pattern": {}}},

                        #E coli
                        "3": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                           "position": {"homopolymer": 0, "random": 1}},
                              "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                            "position": {"homopolymer": 0, "random": 1}},
                              "substitution": {"pattern": {}}},

                        #H sapiens
                        "4": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                           "position": {"homopolymer": 0, "random": 1}},
                              "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                            "position": {"homopolymer": 0, "random": 1}},
                              "substitution": {"pattern": {}}},

                        #M musculus
                        "5": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                           "position": {"homopolymer": 0, "random": 1}},
                              "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                            "position": {"homopolymer": 0, "random": 1}},
                              "substitution": {"pattern": {}}},

                        #Depurination at pH 7 and 293.15K
                        "6": {"deletion": {"pattern": {"A": 0.5, "C": 0.0, "G": 0.5, "T": 0.0},
                                           "position": {"homopolymer": 0.0, "random": 1.0}},
                              "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                            "position": {"homopolymer": 0.5, "random": 0.5}},
                              "substitution": {"pattern": {}}},

                        #Depurination at pH 8 and 253.15K
                        "7": {"deletion": {"pattern": {"A": 0.5, "C": 0.0, "G": 0.5, "T": 0.0},
                                           "position": {"homopolymer": 0.0, "random": 1.0}},
                              "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                            "position": {"homopolymer": 0.5, "random": 0.5}},
                              "substitution": {"pattern": {}}},

                        #Erasure Channel with an error probability of 0.5 percent
                        "8": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                           "position": {"homopolymer": 0.0, "random": 1.0}},
                              "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                            "position": {"homopolymer": 0.5, "random": 0.5}},
                              "substitution": {"pattern": {}}},

                        #Depurination at pH 7 and 193.15K
                        "9":{"deletion": {"pattern": {"A": 0.5, "C": 0.0, "G": 0.5, "T": 0.0},
                                           "position": {"homopolymer": 0.0, "random": 1.0}},
                              "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                            "position": {"homopolymer": 0.5, "random": 0.5}},
                              "substitution": {"pattern": {}}},

                        #White Gaussian Noise with an error probability of 0.5 percent
                        "10":{"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                           "position": {"homopolymer": 0.5, "random": 0.5}},
                              "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                            "position": {"homopolymer": 0.5, "random": 0.5}},
                              "substitution": {"pattern": {"A": {"T": 0.3333, "G": 0.33340000000000003, "C": 0.3333},
                                                           "T": {"A": 0.3333, "G": 0.33340000000000003, "C": 0.3333},
                                                           "C": {"G": 0.3333, "T": 0.33340000000000003, "A": 0.3333},
                                                           "G": {"A": 0.3333, "T": 0.33340000000000003, "C": 0.3333}}}},

                        #Depurination at pH 8 and 193.15K
                        "11":{"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                           "position": {"homopolymer": 0.5, "random": 0.5}},
                              "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                            "position": {"homopolymer": 0.5, "random": 0.5}},
                              "substitution": {"pattern": {}}},

                        #Depurination at pH 8 and 293.15K
                        "12":{"deletion": {"pattern": {"A": 0.5, "C": 0.0, "G": 0.5, "T": 0.0},
                                           "position": {"homopolymer": 0.0, "random": 1.0}},
                              "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                            "position": {"homopolymer": 0.5, "random": 0.5}},
                              "substitution": {"pattern": {}}},

                        #Depurination at pH 7 and 253.15K
                        "13":{"deletion": {"pattern": {"A": 0.5, "C": 0.0, "G": 0.5, "T": 0.0},
                                           "position": {"homopolymer": 0.0, "random": 1.0}},
                              "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                            "position": {"homopolymer": 0.5, "random": 0.5}},
                              "substitution": {"pattern": {}}} ,

                        #None
                        "14": {"deletion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                           "position": {"homopolymer": 0.5, "random": 0.5}},
                              "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
                                            "position": {"homopolymer": 0.5, "random": 0.5}},
                              "substitution": {"pattern": {}}}
                        }


MUTATED_TEXT = []

if __name__ == "__main__":

    with open(fr'{os.getcwd()}\dna-fountain\turkish_anthem.tar.gz.dna_order.txt') as f:
        initial_lines = [line.strip() for line in f if line.strip()]
        seq_objs = [Error_simulation(seq, "storage", attribute = mutation_attributes["11"],
                                     error_rate = err_rates["11"])
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

