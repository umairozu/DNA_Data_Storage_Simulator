import os
from Error_module import Error_simulation

err_rates = {
             "1": {"raw_rate": 0.000043, "substitution": 0.99, "deletion": 0.01, "insertion": 0},
             "2": {"raw_rate": 0.0000024, "substitution": 1, "deletion": 0, "insertion": 0},
             "3": {"raw_rate": 0.0000028, "substitution": 1, "deletion": 0, "insertion": 0},
             "4": {"raw_rate": 0.0000026, "substitution": 0.84, "deletion": 0.08, "insertion": 0.08},
             "5": {"deletion": 0.3333, "substitution": 0.3333, "mismatch": 0.33340000000000003, "raw_rate": 0.0}
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


MUTATED_TEXT = []

if __name__ == "__main__":

    with open(fr'{os.getcwd()}\dna-fountain\turkish_anthem.tar.gz.dna_order.txt') as f:
        initial_lines = [line.strip() for line in f if line.strip()]
        seq_objs = [Error_simulation(seq, "pcr", attribute = mutation_attributes["1"],
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
