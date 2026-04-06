import os


# SOURCE MESA Simulator
# (https://github.com/umr-ds/mesa_dna_sim/blob/master/simulators/error_sources/undesired_subsequences.py)

# More undesired motifs could be added later if needed,
# "De Bruijn Trim rotation graph encoding for reliable DNA storage" mentions some undesired motif's in section 3.3
undesired_sequences = [
    {
        "sequence": "TATAAA",
        "error_prob": "100.0",
        "description": "Eukaryotic promotor recognition motif https://doi.org/10.1016/0022-2836(90)90223-9"
    },
    {
        "sequence": "TTGACA",
        "error_prob": "100.0",
        "description": "Prokaryotic promoter recognition motif https://doi.org/10.1016/j.jmb.2011.01.018"
    },
    {
        "sequence": "TGTATAATG",
        "error_prob": "100.0",
        "description": "Prokaryotic promoter recognition motif https://doi.org/10.1016/j.jmb.2011.01.018"
    },
    {
        "sequence": "GCCACCATGG",
        "error_prob": "100.0",
        "description": "Eukaryotic ribosomal binding site https://doi.org/10.1016/0022-2836(87)90418-9"
    },
    {
        "sequence": "ACCACCATGG",
        "error_prob": "100.0",
        "description": "Eukaryotic ribosomal binding site https://doi.org/10.1016/0022-2836(87)90418-9"
    },
    {
        "sequence": "AATAAA",
        "error_prob": "100.0",
        "description": "Eukaryotic polyadenylation signal https://doi.org/10.1016/0092-8674(87)90292-3"
    },
    {
        "sequence": "TTGTGTGTTG",
        "error_prob": "100.0",
        "description": "Eukaryotic polyadenylation signal https://doi.org/10.1016/0092-8674(87)90292-3"
    },
    {
        "sequence": "ATAACTTCGTATAGCATACATTATACGAAGTTAT",
        "error_prob": "100.0",
        "description": "loxP https://doi.org/10.1016/j.jbiotec.2016.06.033"
    },
    {
        "sequence": "ATAACTTCGTATAGCATACATTATACGAACGGTA",
        "error_prob": "100.0",
        "description": "loxR https://doi.org/10.1016/j.jbiotec.2016.06.033"
    },
    {
        "sequence": "TACCGTTCGTATAGCATACATTATACGAAGTTAT",
        "error_prob": "100.0",
        "description": "loxL https://doi.org/10.1016/j.jbiotec.2016.06.033"
    },
    {
        "sequence": "TACCGTTCGTATAGCATACATTATACGAACGGTA",
        "error_prob": "100.0",
        "description": "loxLR https://doi.org/10.1016/j.jbiotec.2016.06.033"
    },
    {
        "sequence": "TACCGTTCGTATATGGTATTATATACGAAGTTAT",
        "error_prob": "100.0",
        "description": "lox1R https://doi.org/10.1016/j.jbiotec.2016.06.033"
    },
    {
        "sequence": "TACCGTTCGTATATTCTATCTTATACGAAGTTAT",
        "error_prob": "100.0",
        "description": "lox2R https://doi.org/10.1016/j.jbiotec.2016.06.033"
    },
    {
        "sequence": "TACCGTTCGTATAGGATACTTTATACGAAGTTAT",
        "error_prob": "100.0",
        "description": "lox3R https://doi.org/10.1016/j.jbiotec.2016.06.033"
    },
    {
        "sequence": "TACCGTTCGTATATACTATACTATACGAAGTTAT",
        "error_prob": "100.0",
        "description": "lox4R https://doi.org/10.1016/j.jbiotec.2016.06.033"
    },
    {
        "sequence": "TACCGTTCGTATACTATAGCCTATACGAAGTTAT",
        "error_prob": "100.0",
        "description": "lox5R https://doi.org/10.1016/j.jbiotec.2016.06.033"
    },
    {
        "sequence": "ATAACTTCGTATATGGTATTATATACGAACGGTA",
        "error_prob": "100.0",
        "description": "Lox1L https://doi.org/10.1016/j.jbiotec.2016.06.033"
    },
    {
        "sequence": "ATAACTTCGTATAGTATACCTTATACGAAGTTAT",
        "error_prob": "100.0",
        "description": "loxN https://doi.org/10.1038/nature06293"
    },
    {
        "sequence": "ATAACTTCGTATAGTATACATTATACGAAGTTAT",
        "error_prob": "100.0",
        "description": "loxP 511 https://doi.org/10.1093/nar/14.5.2287"
    },
    {
        "sequence": "ATAACTTCGTATAGTACACATTATACGAAGTTAT",
        "error_prob": "100.0",
        "description": "lox 5171 https://doi.org/10.1007/978-981-10-3874-733"
    },
    {
        "sequence": "GCATACAT",
        "error_prob": "100.0",
        "description": "Lox site spacer loxP WT https://doi.org/10.1007/978-981-10-3874-733"
    },
    {
        "sequence": "TGGTATTA",
        "error_prob": "100.0",
        "description": "Lox site spacer lox1 https://doi.org/10.1186/1471-2164-7-73"
    },
    {
        "sequence": "TTCTATCT",
        "error_prob": "100.0",
        "description": "Lox site spacer lox2 https://doi.org/10.1186/1471-2164-7-73"
    },
    {
        "sequence": "GGATACTT",
        "error_prob": "100.0",
        "description": "Lox site spacer lox3 https://doi.org/10.1186/1471-2164-7-73"
    },
    {
        "sequence": "TACTATAC",
        "error_prob": "100.0",
        "description": "Lox site spacer lox4 https://doi.org/10.1186/1471-2164-7-73"
    },
    {
        "sequence": "CTATAGCC",
        "error_prob": "100.0",
        "description": "Lox site spacer lox5 https://doi.org/10.1186/1471-2164-7-73"
    },
    {
        "sequence": "AGGTATGC",
        "error_prob": "100.0",
        "description": "Lox site spacer lox6 https://doi.org/10.1007/978-981-10-3874-733"
    },
    {
        "sequence": "TTGTATGG",
        "error_prob": "100.0",
        "description": "Lox site spacer lox7 https://doi.org/10.1186/1471-2164-7-73"
    },
    {
        "sequence": "GGATAGTA",
        "error_prob": "100.0",
        "description": "Lox site spacer lox8 https://doi.org/10.1186/1471-2164-7-73"
    },
    {
        "sequence": "GTGTATTT",
        "error_prob": "100.0",
        "description": "Lox site spacer lox9 https://doi.org/10.1186/1471-2164-7-73"
    },
    {
        "sequence": "GGTTACGG",
        "error_prob": "100.0",
        "description": "Lox site spacer lox10 https://doi.org/10.1186/1471-2164-7-73"
    },
    {
        "sequence": "TTTTAGGT",
        "error_prob": "100.0",
        "description": "Lox site spacer lox11 https://doi.org/10.1186/1471-2164-7-73"
    },
    {
        "sequence": "GTATACCT",
        "error_prob": "100.0",
        "description": "Lox site spacer loxN https://doi.org/10.1038/nature06293"
    },
    {
        "sequence": "GTACACAT",
        "error_prob": "100.0",
        "description": "Lox site spacer loxP 5171 https://doi.org/10.1007/978-981-10-3874-733"
    },
    {
        "sequence": "GAAGAC",
        "error_prob": "100.0",
        "description": "BbsI"
    },
    {
        "sequence": "GGTCTC",
        "error_prob": "100.0",
        "description": "BsaI"
    },
    {
        "sequence": "CGTCTC",
        "error_prob": "100.0",
        "description": "BsmBI"
    },
    {
        "sequence": "GCTCTTC",
        "error_prob": "100.0",
        "description": "BspQI"
    },
    {
        "sequence": "GCGATG",
        "error_prob": "100.0",
        "description": "BtgZI"
    },
    {
        "sequence": "CGTCTC",
        "error_prob": "100.0",
        "description": "Esp3I"
    },
    {
        "sequence": "GCTCTTC",
        "error_prob": "100.0",
        "description": "SapI"
    },
    {
        "sequence": "CTCGTAGACTGCGTACCA",
        "error_prob": "100.0",
        "description": "Adapter F https://doi.org/10.1186/1746-4811-8-32"
    },
    {
        "sequence": "GACGATGAGTCCTGAGTA",
        "error_prob": "100.0",
        "description": "Adapter R https://doi.org/10.1186/1746-4811-8-32"
    },
    {
        "sequence": "GGTTCCACGTAAGCTTCC",
        "error_prob": "100.0",
        "description": "H1 (HindIII) https://doi.org/10.1016/j.jbiotec.2003.08.005"
    },
    {
        "sequence": "GCGATTACCCTGTACACC",
        "error_prob": "100.0",
        "description": "B4 (BsrGI) https://doi.org/10.1016/j.jbiotec.2003.08.005"
    },
    {
        "sequence": "GCCAGTACATCAATTGCC",
        "error_prob": "100.0",
        "description": "M3 (MfeI) https://doi.org/10.1016/j.jbiotec.2003.08.005"
    },
    {
        "sequence": "AAATAT",
        "error_prob": "100.0",
        "description": "reversed - Eukaryotic promotor recognition motif https://doi.org/10.1016/0022-2836(90)90223-9"
    },
    {
        "sequence": "ACAGTT",
        "error_prob": "100.0",
        "description": "reversed - Prokaryotic promoter recognition motif https://doi.org/10.1016/j.jmb.2011.01.018"
    },
    {
        "sequence": "GTAATATGT",
        "error_prob": "100.0",
        "description": "reversed - Prokaryotic promoter recognition motif https://doi.org/10.1016/j.jmb.2011.01.018"
    },
    {
        "sequence": "GGTACCACCG",
        "error_prob": "100.0",
        "description": "reversed - Eukaryotic ribosomal binding site https://doi.org/10.1016/0022-2836(87)90418-9"
    },
    {
        "sequence": "GGTACCACCA",
        "error_prob": "100.0",
        "description": "reversed - Eukaryotic ribosomal binding site https://doi.org/10.1016/0022-2836(87)90418-9"
    },
    {
        "sequence": "AAATAA",
        "error_prob": "100.0",
        "description": "reversed - Eukaryotic polyadenylation signal https://doi.org/10.1016/0092-8674(87)90292-3"
    },
    {
        "sequence": "GTTGTGTGTT",
        "error_prob": "100.0",
        "description": "reversed -Eukaryotic polyadenylation signal https://doi.org/10.1016/0092-8674(87)90292-3"
    },
    {
        "sequence": "TATTGAAGCATATTACATACGATATGCTTCAATA",
        "error_prob": "100.0",
        "description": "reversed - loxP https://doi.org/10.1016/j.jbiotec.2016.06.033"
    },
    {
        "sequence": "ATGGCAAGCATATTACATACGATATGCTTCAATA",
        "error_prob": "100.0",
        "description": "reversed - loxR https://doi.org/10.1016/j.jbiotec.2016.06.033"
    },
    {
        "sequence": "TATTGAAGCATATTACATACGATATGCTTGCCAT",
        "error_prob": "100.0",
        "description": "reversed - loxL https://doi.org/10.1016/j.jbiotec.2016.06.033"
    },
    {
        "sequence": "ATGGCAAGCATATTACATACGATATGCTTGCCAT",
        "error_prob": "100.0",
        "description": "reversed - loxLR https://doi.org/10.1016/j.jbiotec.2016.06.033"
    },
    {
        "sequence": "TATTGAAGCATATATTATGGTATATGCTTGCCAT",
        "error_prob": "100.0",
        "description": "reversed - lox1R https://doi.org/10.1016/j.jbiotec.2016.06.033"
    }
]

# REMOVING those lines from the file that may have any of the above undesired sequences
"""IDEA*
Come back after decoding module is started, check instead of removing the entire line, what if we just
remove that particular sequence from the line and see if the line/Oligo is valid or not??

- I guess Oligo will become invalid ---> if undesired sequence is removed from the center or the line (in between Primers),
    then the strand will break and PCR won't work
- otherwise if a seq is removed from amongst the adapters or the seed then the identification of chunks and Oligo itself is not possible

    Think around this concept!!
"""


input_file = "Turkish_anthem.tar.gz.dna_order.txt"
path = fr'{os.getcwd()}\dna-fountain\{input_file}'
with open(path, "r") as file:
    lines = file.readlines()
    sequences = [item["sequence"] for item in undesired_sequences]
    file.seek(0)

    filtered_lines = []
    undesired_seq_count = 0
    for line in lines:
        found = False
        for seq in sequences:
            if seq in line:
                found = True
                undesired_seq_count += 1

        if not found:
            filtered_lines.append(line)

    # alternate way of writing
    #filtered_lines = [line for line in lines if not any(seq in line for seq in sequences)]

    print(f"undesired sequences count: {undesired_seq_count}")

    with open(path,"w+") as f:
        f.writelines(filtered_lines)


if __name__ == "__main__":
    print(len(undesired_sequences))


