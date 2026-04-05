

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

import numpy as np

from Enzyme_Addition import pF_length, pR_length, payload_length
from Error_module import Error_simulation
import Levenshtein
from Bio import Align
from Bio.Seq import Seq


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

BASE_DIR = fr'{os.getcwd()}\dna-fountain'
os.makedirs(fr'{BASE_DIR}\sequencing', exist_ok=True)
PCR_DIR = fr'{os.getcwd()}\dna-fountain\pcr'
SEQ_DIR = fr'{os.getcwd()}\dna-fountain\sequencing'



def get_clean_segment(parent_aligned, child_aligned, start_pos, target_len):
    valid_chars = 0
    chars_consumed = 0

    for i in range(start_pos, len(parent_aligned)):
        if valid_chars == target_len:
            break

        chars_consumed += 1
        if parent_aligned[i] != '-':
            valid_chars += 1

    end_pos = start_pos + chars_consumed

    clean_segment = "".join(child_aligned[start_pos:end_pos]).replace('-', '')

    return clean_segment, end_pos


MUTATED_TEXT = []

if __name__ == "__main__":

    mode = input("Enter Sequencing Mode:  --help [1. Illumina, 2. PacBio, 3. Nanopore ]\n")
    assert mode in ["1","2","3"]
    if mode == "1":
        method = input("Enter method for illumina sequencing:  --help [1. Single-end, 2. Paired-end]\n")
        assert method in ["1", "2"]
    elif mode == "2":
        method = input("Enter method for PacBio sequencing:  --help [3. CCS, 4. Subread]\n")
        assert method in ["3", "4"]
    elif mode == "3":
        method = input("Enter method for Nanopore sequencing:  --help [5. 1D, 6. 2D]\n")
        assert method in ["5", "6"]
    else:
        mode = "1"
        method = "1"
        print("Default sequencing method chosen [illumina single-end sequencing]")

    # sampling for sequencing
    while True:
        user_input = input("Enter sequencing sampling fraction (0-100%): \n")

        try:
            value = float(user_input)
            if 0 <= value <= 100:
                SAMPLING_FRAC = value / 100
                break
            else:
                print("Please enter a number between 0 and 100.")
        except ValueError:
            print("Invalid input. Please enter a numerical value (e.g., 5.5).")

    # target_reads & READ_LENGTH input
    while True:
        user_input02 = input("Enter a number for total reads required (e.g 100K, 10M etc): \n")
        user_input03 = input("Enter a number for Read Length (e.g 120, 150 etc): \n")
        try:
            value02 = int(user_input02)
            value03 = int(user_input03)
            if value02 and value03:
                TARGET_READS = value02
                READ_LENGTH = value03
                break
            else:
                print("Please enter a number greater than 0.")
        except ValueError:
            print("Invalid input. Please enter a numerical value (e.g., 5.5).")


    POOL_COUNT = 0 # total molecule in input pcr file
    with open(fr'{PCR_DIR}\pcr_final.txt') as f:
        next(f)
        data = []
        for line in f:
            if "," in line.strip():
                count, seq, length = line.split(",")
                POOL_COUNT += int(count)
                data.append({ 'count': int(count.strip()), 'seq': seq.strip(), 'length': length.strip() })


    #Collecting original oligo's given for order
    ORIGINAL_ORDER_OLIGOS = []
    with open(fr'{BASE_DIR}\original_order_file.txt') as file:
        for line in file:
            ORIGINAL_ORDER_OLIGOS.append(line)

    #Collecting best matching Parent oligo given a sequence
    PARENT_CHILD_INFO = []
    for item in data:
        count = item['count']
        seq = item['seq']
        length = item['length']

        best_match = min(ORIGINAL_ORDER_OLIGOS, key = lambda x: Levenshtein.distance(seq,x))
        similarity = Levenshtein.ratio(seq, best_match)

        PARENT_CHILD_INFO.append((count, best_match.strip(), seq.strip(), similarity))
        #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # [(child's copy count, parent, child, parent-child similarity score)]

        #print(f"seq 0: {seq} \nseq 1: {best_match}\nSimilarity: {similarity * 100:.2f}%")

    with open(fr'{SEQ_DIR}\PARENT_CHILD_INFO.txt',"w") as f:
        f.write("child copy count, parent, child\n")
        for item in PARENT_CHILD_INFO:
            f.write(f"{item}\n")

    parent_pF_length = pF_length
    parent_insert_length = payload_length
    parent_pR_length = pR_length

    #Calculating child's pF, pR and insert given a parent --> for illumina SE and PE sequencing
    CHILD_INFO = []
    for item in PARENT_CHILD_INFO:

        child_copy = item[0]
        parent_seq = item[1]
        child_seq = item[2]
        parent_child_score = item[3]

        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        alignments = aligner.align(child_seq, parent_seq)

        best_match = alignments[0]

        alignment_lines = str(best_match).split("\n")
        parent = list(alignment_lines[2])
        child = list(alignment_lines[0])

        child_pF, pos1 = get_clean_segment(parent, child, 0, parent_pF_length)
        #print(f"Child_pF: {child_pF}")

        child_insert, pos2 = get_clean_segment(parent, child, pos1, parent_insert_length)
        #print(f"Child_insert: {child_insert}")

        child_pR = "".join(child[pos2:]).replace('-', '')
        #print(f"Child_pR: {child_pR}")

        CHILD_INFO.append((parent_seq,child_copy,child_pF,child_insert,child_pR, parent_child_score))
        #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # [(parent, child's copy count, child's pF, child's insert, child's pR, parent_child_score)]

        #print("-----......-----")


    with open(fr'{SEQ_DIR}\CHILD_INFO.txt',"w") as f:
        f.write("parent, child copy count, child pF, child insert, child pR, score\n")
        for item in CHILD_INFO:
            f.write(f"{item}\n")

    #Quality Check on child: 1) Alignment score > 0.80
    # (2) parent-child length difference with +- 10 [optional] ---> consider lengths of mutated seqs, variants and chimeras,
    # they might be filtered, so removing them not so good unless they are very short as compared to parent!!
    LIST_01 = []
    alignment_THRESHOLD = 0.80
    len_VALVE = 10
    with open(fr'{SEQ_DIR}\CHILD_INFO.txt') as f:
        next(f)
        for item in CHILD_INFO:
            child_seq = item[3] # child's insert
            score = float(item[5])
            #if (score > 0.80) and (abs(payload_length - len(child_seq)) < 5):(MUCH STRICTER)
            if (score > alignment_THRESHOLD) and (abs(payload_length - len(child_seq)) < len_VALVE): # parent-child alignment is okay!!
                LIST_01.append(item)
            else:
                pass

    with open(fr'{SEQ_DIR}\pre_sequencing_filter.txt',"w") as f:
        f.write("parent, child copy count, child pF, child insert, child pR, score\n")
        for item in LIST_01:
            f.write(f"{item}\n")

    #print(len(data))
    #print(len(PARENT_CHILD_INFO))
    #print(len(new_LIST))


    #Sequence input pool --> after sampling
    LIST_02 = []
    with open(fr'{SEQ_DIR}\sequencing_sampled_pool.txt', "w") as f:
        f.write("parent, child copy count, child pF, child insert, child pR, score, length\n")
        for item in LIST_01:
            parent = item[0]
            child_count = item[1]
            child_pF = item[2]
            child_insert = item[3]
            child_pR = item[4]
            score = item[5]

            seq = child_pF + child_insert + child_pR
            length = len(seq)
            count = int(SAMPLING_FRAC * child_count)

            LIST_02.append((parent,count,child_pF,child_insert,child_pR,length))
            f.write(f"{parent},{count},{child_pF},{child_insert},{child_pR},{length}\n")


    # LIST_02 currently stores:
    # (parent, sampled_count, child_pF, child_insert, child_pR, length)
    templates = []
    weights = []

    for item in LIST_02:
        sampled_count = int(item[1])
        if sampled_count <= 0:
            continue

        pF = item[2]
        insert = item[3]
        pR = item[4]

        templates.append((pF, insert, pR))
        weights.append(sampled_count)

    weights = np.array(weights, dtype=float)
    weights = weights / weights.sum()

    # allocating how many reads each template gets
    read_allocations = np.random.multinomial(TARGET_READS, weights)


    # generating sequencing reads for the chosen method
    READ_SE = []
    READ_PE = []
    with open(fr'{SEQ_DIR}\sequencing_sampled_pool.txt') as f:
        next(f)
        for item, n_reads in zip(LIST_02, read_allocations):

            pF = item[2]
            insert = item[3]
            pR = item[4]
            sequence = pF + insert + pR

            insert_length = len(insert)
            for _ in range(n_reads):
                if mode == '1':
                    if method == '1': # illumina's SE sequencing
                        if READ_LENGTH > insert_length:
                            read1 = (insert + pR)[:READ_LENGTH]
                        elif READ_LENGTH == insert_length:
                            read1 = insert
                        else:
                            read1 = insert[:READ_LENGTH]
                        READ_SE.append(read1)
                    else:           # illumina's PE sequencing
                        if READ_LENGTH > insert_length:
                            read1 = (insert + pR)[:READ_LENGTH]
                            read2 = str(Seq((pF + insert)[-READ_LENGTH:]).reverse_complement())
                        elif READ_LENGTH == insert_length:
                            read1 = insert
                            read2 = str(Seq(insert).reverse_complement())
                        else:
                            read1 = insert[:READ_LENGTH]
                            read2 = str(Seq(insert[-READ_LENGTH:]).reverse_complement())
                        READ_SE.append(read1)
                        READ_PE.append(read2)

    # Writing reads back to their respective fastq files:
    chosen_file = ''
    if method == '1': # For READ_SE:
        with open(fr'{SEQ_DIR}\Read_SE.txt', "w") as file:
            for read_id, read in enumerate(READ_SE):
                file.write(f"@read{read_id}\n{read}\n+\n!!!!!!\n")
        chosen_file = fr'{SEQ_DIR}\Read_SE.txt'
    elif method == '2': # For READ_PE:
        with open(fr'{SEQ_DIR}\Read_PE.txt', "w") as file:
            read_id = 1
            for r1, r2 in zip(READ_SE, READ_PE):
                file.write(f"@read{read_id}/1\n{r1}\n+\n!!!!!!\n@read{read_id}/2\n{r2}\n+\n!!!!!!\n")
                read_id += 1
        chosen_file = fr'{SEQ_DIR}\Read_PE.txt'
    else:
        raise ValueError("Unsupported method")

    print(f"Length SE: {len(READ_SE)}")
    print(f"Length PE: {len(READ_PE)}")


    # for fixed quality score
    """
    def phred_char(q=30):
        return chr(q + 33)  # Q30 -> '?'

    def qual_string(s, q=30):
        return phred_char(q) * len(s)
    """

    def quality_profile(read_length, q_start=35, q_end=28):
        if read_length == 1:
            return chr(q_start + 33)
        qs = []
        for i in range(read_length):
            q = round(q_start + (q_end - q_start) * i / (read_length - 1))
            qs.append(chr(q + 33))
        return "".join(qs)

    # Making sequencing object + Applying sequencing errors
    with open(chosen_file) as f:
        lines = [line.strip() for i, line in enumerate(f) if i % 4 == 1]
        seq_objs = [Error_simulation(seq, "sequencing", attribute = mutation_attributes[mode],
                                     error_rate = err_rates[method])
                    for seq in lines
                    ]

    VALVE = 10
    filename = "sequencing_SE.fastq" if method == '1' else "sequencing_PE.fastq"
    filepath = fr'{SEQ_DIR}\{filename}'
    if mode == '1':
        with open(filepath, "w") as f:
            for _ in range(1, VALVE):
                MUTATED_TEXT = []
                for sE in seq_objs:
                    sE.run_mutations()
                    MUTATED_TEXT.append(sE.seq)

            MUTATED_TEXT_FINAL = []
            for read in MUTATED_TEXT:
                clean_read = "".join(read).replace(" ","")
                MUTATED_TEXT_FINAL.append(clean_read)

            # Writing fast.q file with mutated reads
            if method == '1':
                for read_id, read in enumerate(MUTATED_TEXT_FINAL, start = 1):
                    qual1 = quality_profile(len(read), 35, 28)
                    f.write(f"@read{read_id}\n{read}\n+\n{qual1}\n")
            else:
                read_id = 1
                for i in range(0, len(MUTATED_TEXT_FINAL), 2):

                    read1 = MUTATED_TEXT_FINAL[i]
                    read2 = MUTATED_TEXT_FINAL[i + 1]

                    qual1 = quality_profile(len(read), 35, 28)
                    qual2 = quality_profile(len(read2), 33, 24)

                    f.write(f"@read{read_id}/1\n{read1}\n+\n{qual1}\n@read{read_id}/2\n{read2}\n+\n{qual2}\n")
                    read_id += 1

        # Writing Additional R1.fastq and R2.fastq separately as well
        if method == '2':
            with open(fr"{SEQ_DIR}\\sequencing_R1.fastq", "w") as f1, open(fr"{SEQ_DIR}\\sequencing_R2.fastq", "w") as f2:
                read_id = 1
                for i in range(0, len(MUTATED_TEXT_FINAL), 2):

                    read1 = MUTATED_TEXT_FINAL[i]
                    read2 = MUTATED_TEXT_FINAL[i + 1]

                    qual1 = quality_profile(len(read), 35, 28)
                    qual2 = quality_profile(len(read2), 33, 24)

                    f1.write(f"@read{read_id}/1\n{read1}\n+\n{qual1}\n")
                    f2.write(f"@read{read_id}/2\n{read2}\n+\n{qual2}\n")
                    read_id += 1

    print("Sequencing.py run")