
import os


# --------------------
# constants / paths
# --------------------
BASE_DIR = fr"{os.getcwd()}\dna-fountain"
INPUT_PATH = fr"{BASE_DIR}\turkish_anthem.tar.gz.dna"
ORDER_PATH = fr"{BASE_DIR}\turkish_anthem.tar.gz.dna_order.txt"
FASTA_PATH = fr"{BASE_DIR}\final_file.FASTA"

primer_F = "TGGCTCATTT"
primer_R = "ATAAATGACC"


# --------------------
# read helpers
# --------------------
def read_payload_sequences(path):
    lines = []
    with open(path) as file:
        for line in file:
            if line and line[0] != ">":
                lines.append(line.strip())
    return lines


def read_identifiers(path):
    identifiers = []
    with open(path) as file:
        for line in file:
            if line and line[0] == ">":
                identifiers.append(line)
    return identifiers


# --------------------
# write helpers
# --------------------
def write_dna_order_file(path, oligos):
    with open(path, "w") as file:
        for line in oligos:
            file.write(line + "\n")


# final FASTA file creation
# writing final Oligo sequence [Primer_F,RS,Payload,seed,cutting_site,binding_site] back with their
# identifiers as FASTA file for later use (Debugging/ decode etc.)
# --------------------
# build helpers
# --------------------
def build_oligos(payloads):
    return [primer_F + seq + primer_R for seq in payloads]

def write_final_fasta(path, identifiers, oligos):
    with open(path, "w") as file:
        for ident, oligo in zip(identifiers, oligos):
            file.write(ident)
            file.write(oligo + "\n")


# --------------------
# metadata for imports
# --------------------
_payloads = read_payload_sequences(INPUT_PATH)
_oligos = build_oligos(_payloads)

orig_length = len(_oligos[0]) if _oligos else 0
pF_length = len(primer_F)
pR_length = len(primer_R)


# --------------------
# runner
# --------------------
def prepare_enzyme_files():
    identifiers = read_identifiers(INPUT_PATH)
    write_dna_order_file(ORDER_PATH, _oligos)
    write_final_fasta(FASTA_PATH, identifiers, _oligos)


if __name__ == "__main__":
    prepare_enzyme_files()
    print("Enzyme_Addition.py run")

    # Users own primer's Input taken and added to the file
    """
    pattern = re.compile(r'^[ACGTacgt]{20}$')

    primer_F = input("Enter Forward primer: \n").upper()
    while not pattern.fullmatch(primer_F):
        primer_F = input("Enter Forward primer: \n").upper()


    primer_R = input("Enter Reverse primer: \n").upper()
    while not pattern.fullmatch(primer_R):
        primer_R = input("Enter Reverse primer: \n").upper()


    print("Valid Inputs --> Proceed forward")


    with open(new_path) as file:
        lines_2 = file.readlines()

    with open(new_path,"w+") as file:
        for line in lines_2:
            line = primer_F + line.strip() + primer_R
            file.writelines(line + "\n")
    """