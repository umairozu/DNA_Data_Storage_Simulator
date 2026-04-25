import argparse
import json
import os
import re

from GC_content import global_gc_content

# --------------------
# constants / paths
# --------------------
BASE_DIR = fr"{os.getcwd()}"
INPUT_PATH = fr"{BASE_DIR}/dna-fountain/test_file.tar.gz.dna"
ORDER_PATH = fr"{BASE_DIR}/dna-fountain/test_file.tar.gz.dna_order.txt"
FASTA_PATH = fr"{BASE_DIR}/dna-fountain/final_file.FASTA"
ORIGINAL_SEQS_PATH = fr"{BASE_DIR}/ORIGINAL_SEQS.txt"

# Users own primer's Input taken and added to the file
_p = argparse.ArgumentParser(description="Enzyme_Addition.py run")

_p.add_argument("--pf", default="ATAAATGACCTGCCGTGCAA", help="Forward Primer")
_p.add_argument("--pr", default="ACCGATTGTGAAATGAGCCA", help="Reverse Primer")
_p.add_argument("--gc_min", type=float, default=0.45, help="Min gc content")
_p.add_argument("--gc_max", type=float, default=0.55, help="Max gc content")

args = _p.parse_args()

def validate_primer(primer_seq, min_gc, max_gc, label):

    if not re.fullmatch(r'^[ACGTacgt]+$', primer_seq):
        _p.error(f"{label} contains invalid characters!")

    gc_content = global_gc_content(primer_seq) / 100
    print(f"{label} GC content: {gc_content}")

    if not (min_gc <= gc_content <= max_gc):
        _p.error(f"{label} ({primer_seq}) fails GC bounds ({min_gc}-{max_gc})")

    return primer_seq.upper()

primer_F = validate_primer(args.pf, args.gc_min, args.gc_max, "Forward Primer")
primer_R = validate_primer(args.pr, args.gc_min, args.gc_max, "Reverse Primer")


print("Valid Inputs --> You may proceed to synthesis.py")

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
payload_length = orig_length - (pF_length + pR_length)


# --------------------
# runner
# --------------------
def save_metadata():
    metadata = {
        "primer_F": primer_F,
        "primer_R": primer_R,
        "orig_length": orig_length,
        "pF_length": pF_length,
        "pR_length": pR_length,
        "payload_length": payload_length
    }
    with open(fr"{BASE_DIR}/metadata.json", "w") as f:
        json.dump(metadata, f, indent=5)

def prepare_enzyme_files():
    identifiers = read_identifiers(INPUT_PATH)
    write_dna_order_file(ORDER_PATH, _oligos)
    write_final_fasta(FASTA_PATH, identifiers, _oligos)
    save_metadata()

def fasta_to_txt(fasta_path, output_path):
    sequences = []
    with open(fasta_path) as file:
        for line in file:
            if '>' not in line:
                sequences.append(line)

    with open(output_path, 'w') as file:
        for item in sequences:
            file.write(fr'{item}')


if __name__ == "__main__":
    
    prepare_enzyme_files()
    fasta_to_txt(INPUT_PATH,ORIGINAL_SEQS_PATH)
    print("Enzyme_Addition.py run")