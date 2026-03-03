

err_rates = {
             "1": {"raw_rate": 0.0001, "mismatch": 0.15, "deletion": 0.7, "insertion": 0.15},
             "2": {"raw_rate": 0.000125, "mismatch": 0.15, "deletion": 0.7, "insertion": 0.15},
             "3": {"raw_rate": 0.00011, "mismatch": 0.2, "deletion": 0.6, "insertion": 0.2},
             "4": {"raw_rate": 0.0017, "mismatch": 0.2, "deletion": 0.6, "insertion": 0.2},
             "5": {"raw_rate": 0.00125, "mismatch": 0.2, "deletion": 0.6, "insertion": 0.2},
             "6": {"raw_rate": 0.00033, "mismatch": 0.2, "deletion": 0.6, "insertion": 0.2},
             "7": {"raw_rate": 0.000025, "mismatch": 0.2, "deletion": 0.6, "insertion": 0.2},
             "8": {"raw_rate": 0.0004, "mismatch": 0.2, "deletion": 0.6, "insertion": 0.2},
             "9": {"raw_rate": 0, "mismatch": 0.3, "deletion": 0.3, "insertion": 0.3}
             }

mutation_attributes = {
                        #MutS
                        "1": {"deletion": {"position": {"homopolymer": 0.0, "random": 1}, "pattern": {"G": 0.2, "C": 0.2, "A": 0.4, "T": 0.2}},"insertion": {"position": {"homopolymer": 0, "random": 1},"pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},"mismatch": {}},

                        #Consensus Shuffle
                        "2": {"deletion": {"position": {"homopolymer": 0.0, "random": 1}, "pattern": {"G": 0.2, "C": 0.2, "A": 0.4, "T": 0.2}},"insertion": {"position": {"homopolymer": 0, "random": 1},"pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},"mismatch": {}},

                        #NGS-based error correction
                        "3": {"deletion": {"position": {"homopolymer": 0.0, "random": 1}, "pattern": {"G": 0.2, "C": 0.2, "A": 0.4, "T": 0.2}},"insertion": {"position": {"homopolymer": 0, "random": 1},"pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},"mismatch": {}},

                        #High-temperature ligation/hybridization based error correction
                        "4": {"deletion": {"position": {"homopolymer": 0.0, "random": 1}, "pattern": {"G": 0.2, "C": 0.2, "A": 0.4, "T": 0.2}},"insertion": {"position": {"homopolymer": 0, "random": 1},"pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},"mismatch": {"pattern": {"AAAA": {"ACGT":1.0}, "ACCC": {"ACGC":1.0}}}},

                        #ErrASE
                        "5": {"deletion": {"position": {"homopolymer": 0.0, "random": 1}, "pattern": {"G": 0.2, "C": 0.2, "A": 0.4, "T": 0.2}},"insertion": {"position": {"homopolymer": 0, "random": 1},"pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},"mismatch": {}},

                        #Nuclease-based error correction
                        "6": {"deletion": {"position": {"homopolymer": 0.0, "random": 1}, "pattern": {"G": 0.2, "C": 0.2, "A": 0.4, "T": 0.2}},"insertion": {"position": {"homopolymer": 0, "random": 1},"pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},"mismatch": {}},

                        #ErrASE
                        "7": {"deletion": {"pattern": {"A": 0.4, "C": 0.2, "G": 0.2, "T": 0.2}, "position": {"homopolymer": 0, "random": 1}}, "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}, "position": {"homopolymer": 0, "random": 1}}, "mismatch": {"pattern": {}}},

                        #Oligo Hybridization based error correction
                        "8": {"deletion": {"pattern": {"A": 0.4, "C": 0.2, "G": 0.2, "T": 0.2}, "position": {"homopolymer": 0, "random": 1}}, "insertion": {"pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}, "position": {"homopolymer": 0, "random": 1}}, "mismatch": {"pattern": {}}},

                        #None
                        "9": {"deletion": {"position": {"homopolymer": 0.0, "random": 1}, "pattern": {"G": 0.2, "C": 0.2, "A": 0.4, "T": 0.2}},"insertion": {"position": {"homopolymer": 0, "random": 1},"pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}},"mismatch": {}}

                        }