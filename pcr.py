import os
import random
import numpy as np
from GC_content import gc_error_probability
from Enzyme_Addition import primer_F, primer_R,orig_length,pF_length,pR_length
from Error_module import Error_simulation


err_rates = {
             "1": {"raw_rate": 0.000043, "substitution": 0.99, "deletion": 0.01, "insertion": 0}, #Taq
             "2": {"raw_rate": 0.0000024, "substitution": 1, "deletion": 0, "insertion": 0}, #Pwo
             "3": {"raw_rate": 0.0000028, "substitution": 1, "deletion": 0, "insertion": 0}, #Pfu
             "4": {"raw_rate": 0.0000026, "substitution": 0.84, "deletion": 0.08, "insertion": 0.08}, #Phusion
             "5": {"deletion": 0.3333, "substitution": 0.33340000000000003, "insertion": 0.3333, "raw_rate": 0.0} #None
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


def amp_factor(eff_i):
    return 1 + eff_i

def seq_check(eff_i, sequence):
    num_primer_mismatch = 0
    num_primer_mismatch += sum(1 for a, b in zip(sequence[:pF_length], primer_F) if a != b)
    num_primer_mismatch += sum(1 for a, b in zip(sequence[-pR_length:], primer_R) if a != b)

    # dropping efficiency if GC deviates from 45% - 55% range
    gc_error = gc_error_probability(sequence)
    if 0.45 <= gc_error <= 0.55:
        pass
    else:
        eff_i *= 0.4
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
        """OPTIONAL: Local GC Content affect"""

    # dropping efficiency if mismatches at the primer sites
    if num_primer_mismatch < 1:
        eff_i = eff_i
    elif num_primer_mismatch < 2:
        eff_i *= 0.7
    elif num_primer_mismatch < 3:
        eff_i *= 0.4
    else:
        eff_i *= 0.0  # the oligo won't amplify if primer mutation > 3
    return eff_i



path = fr'{os.getcwd()}\dna-fountain'

# PCR pre-filtering --> removing non_specific amplicons as a cleanup for sequencing
with open(fr'{path}\storage_file_1.txt') as f:
    """
    #200823,AATGGTTTACCCATA
    #count = [200823], seq [AATGGTTTACCCATA]
    pairs = [line.strip().split(",") for line in f if line.strip()]
    count, seq = zip(*pairs)
    """
    data = []
    for line in f:
        if "," in line.strip():
            count, seq = line.split(",")
            data.append({'count': int(count.strip()), 'seq': seq.strip()})

dropouts = []
filtered_lines = []

for index in data:
    count = index['count']
    seq = index['seq']
    if len(seq) == orig_length:
        if seq[:pF_length] != primer_F or seq[-pR_length:] != primer_R:
            dropouts.append((count,seq))
        else:
            filtered_lines.append((count,seq))
    else:
        if abs(len(seq) - orig_length) >= 10:
            dropouts.append((count,seq))
        else:
            filtered_lines.append((count,seq))

with open(fr'{path}\pcr_dropouts.txt', "w") as f:
    f.write("count, sequence, length\n")
    for c, s in dropouts:
        l = len(s)
        f.write(f"{c},{s},{l}\n")

with open(fr'{path}\pcr_filtered.txt', "w") as f:
    f.write("count, sequence, length\n")
    for c, s in filtered_lines:
        l = len(s)
        f.write(f"{c},{s},{l}\n")


# PCR sampling
sampling_frac = float(input("Enter PCR sampling fraction (in %age): ")) / 100

with open(fr'{path}\pcr_filtered.txt') as f:
    next(f)
    data = []
    for line in f:
        if "," in line.strip():
            count, seq, length = line.split(",")
            data.append({ 'count': int(count.strip()), 'seq': seq.strip(), 'length': length.strip() })

for index in data:
    count = index['count']
    seq = index['seq']
    length = index['length']
    bnd = np.random.binomial(int(count),sampling_frac, size = int(length))
    bnd_choice = np.random.choice(bnd)
    index['count'] = bnd_choice

with open(fr'{path}\pcr_sampled.txt', "w") as f:
    f.write("sampled_count, sequence, length\n")
    for index in data:
        count = index['count']
        seq = index['seq']
        length = index['length']
        f.write(f'{count},{seq},{length}\n')

# PCR amplification
num_cycles = int(input("Enter the desired PCR cycle to run: ")) # e.g 30
c1 = num_cycles // 3 # e.g 10
c2 = c1 + num_cycles // 3 # e.g 10 + 10 = 20
c3 = num_cycles # e.g 30

max_yield = int(input("Enter the maximum pcr yield expected:")) # so new pool will be original molecules + pcr pool/yield

E0_i = 0.9  # initial efficiency
"""Calculating per oligo Efficiency in the sample based on [PRIMER BINDING, GC CONTENT] of the oligo"""
with open(fr'{path}\pcr_sampled.txt') as f:
    next(f)
    data = []
    for line in f:
        if "," in line.strip():
            count, seq, length = line.split(",")
            data.append({ 'count': int(count.strip()), 'seq': seq.strip(), 'length': length.strip(), 'eff': seq_check(E0_i,seq) })

    remaining_yield = max_yield

    for c in range(0, c1): # on good conditions, exponential efficiency in c1
        if remaining_yield <= 0:
            break

        new_copy_count = []
        requested_copies = []

        for index_i in data:
            #eff_01 = E0_i
            count_i = index_i['count']
            seq_i = index_i['seq']
            eff = index_i['eff']
            # Turning sequence check at every cycle off as it might cause unnecessary efficiency collapse, rather computing it once at the start
            #eff = seq_check(eff, seq_i)
            if not (0 <= eff <= 1):
                n_i = count_i
                delta = 0
            else:
                amp_i = amp_factor(eff)
                n_i = count_i * amp_i
                delta = n_i - count_i
            requested_copies.append(delta)
            new_copy_count.append(n_i)
            index_i['eff'] = eff # this could store -ve efficiency, but we are dealing with it
        total_demand = sum(requested_copies)
        if total_demand == 0: # nothing to amplify in this cycle
            break
        if remaining_yield >= total_demand:
            for i, item in enumerate(data):
                item['count'] = int(new_copy_count[i])
            remaining_yield -= total_demand
        else: # if demand higher than what we have, we first check if yield remains then distribute proportionally, otherwise break from the cycle
            if remaining_yield:
                sum_amount_given = 0
                for i, item in enumerate(data):
                    fraction = requested_copies[i] / total_demand
                    amount_given = fraction * remaining_yield
                    sum_amount_given += amount_given
                    item['count'] += int(amount_given)
                remaining_yield -= sum_amount_given
            else:
                break

    for c in range(c1, c2): # efficiency become somewhat linear--> slow linear decrease
        if remaining_yield <= 0:
            break

        new_copy_count = []
        requested_copies = []

        for index_i in data:
            #eff_02 = eff_01
            count_i = index_i['count']
            seq_i = index_i['seq']
            eff = index_i['eff']
            eff = eff - (eff / 10)
            #eff = seq_check(eff, seq_i)
            if not (0 <= eff <= 1):
                n_i = count_i
                delta = 0
            else:
                amp_i = amp_factor(eff)
                n_i = count_i * amp_i
                delta = n_i - count_i
            requested_copies.append(delta)
            new_copy_count.append(n_i)
            index_i['eff'] = eff
        total_demand = sum(requested_copies)
        if total_demand == 0:
            break
        if remaining_yield >= total_demand:
            for i, item in enumerate(data):
                item['count'] = int(new_copy_count[i])
            remaining_yield -= total_demand
        else:
            if remaining_yield:
                sum_amount_given = 0
                for i, item in enumerate(data):
                    fraction = requested_copies[i] / total_demand
                    amount_given = fraction * remaining_yield
                    sum_amount_given += amount_given
                    item['count'] += int(amount_given)
                remaining_yield -= sum_amount_given
            else:
                break

    for c in range(c2, c3):  # efficiency plateaus
        if remaining_yield <= 0:
            break

        new_copy_count = []
        requested_copies = []

        for index_i in data:
            #eff_03 = 0
            count_i = index_i['count']
            seq_i = index_i['seq']
            eff = index_i['eff']
            #eff = seq_check(eff, seq_i)
            amp_i = 0.2 + eff  # very minimal amplification in plateau stage
            if amp_i < 1.0: # preventing counts decrements incase eff is less than 0.8
                n_i = count_i
                delta = 0
            else:
                n_i = count_i * amp_i
                delta = n_i - count_i

            requested_copies.append(delta)
            new_copy_count.append(n_i)
            index_i['eff'] = eff
        total_demand = sum(requested_copies)
        if total_demand == 0:
            break

        if remaining_yield >= total_demand:
            for i, item in enumerate(data):
                item['count'] = int(new_copy_count[i])
            remaining_yield -= total_demand
        else:
            if remaining_yield:
                sum_amount_given = 0
                for i, item in enumerate(data):
                    fraction = requested_copies[i] / total_demand
                    amount_given = int(fraction * remaining_yield)
                    sum_amount_given += amount_given
                    item['count'] += amount_given
                remaining_yield -= sum_amount_given
            else:
                break

"""
    with open(path) as f:
        initial_lines = [line.strip() for line in f if line.strip()]
        seq_objs = [Error_simulation(seq, "pcr", attribute = mutation_attributes["1"],
                                     error_rate = err_rates["1"])
                    for seq in initial_lines
                    ]

    count = 1
    while count < 30:
        MUTATED_TEXT.clear()
        for sE in seq_objs:
            #sE.reset_visited() # See important Notice for info
            sE.run_mutations()
            MUTATED_TEXT.append(sE.seq)

        with open(fr'{os.getcwd()}\dna-fountain\storage_file_1.txt', "w") as f:
            f.write("\n".join(MUTATED_TEXT) + "\n")
        count += 1
"""

MUTATED_TEXT = []

if __name__ == "__main__":
    print("HI")