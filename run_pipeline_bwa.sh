#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")"

mkdir -p cutadapt_files bwa_files

echo "Starting Bioinformatics Pipeline (Cutadapt -> BWA -> DNA Fountain)..."

deactivate 2>/dev/null || true
conda deactivate 2>/dev/null || true

echo "Activating conda environment: dna_tools..."
source ~/miniforge3/etc/profile.d/conda.sh
conda activate dna_tools

# Pick whichever aligner is installed.
if command -v bwa-mem2 >/dev/null 2>&1; then
    BWA_BIN="bwa-mem2"
elif command -v bwa >/dev/null 2>&1; then
    BWA_BIN="bwa"
else
    echo "Error: neither bwa-mem2 nor bwa is available in PATH"
    exit 1
fi

echo "Using aligner: ${BWA_BIN}"

echo "Step 1: Trimming adapters with cutadapt..."

if [ -s sequencing/sequencing_R2.fastq ]; then
    MODE="PE"
    echo "Detected paired-end input (R1 + R2)..."
    cutadapt \
      -a ACCGATTGTGAAATGAGCCAX \
      -A TTGCACGGCAGGTCATTTATX \
      -e 0.1 -m 21 \
      -o cutadapt_files/trimmed_R1.fastq \
      -p cutadapt_files/trimmed_R2.fastq \
      sequencing/sequencing_R1.fastq sequencing/sequencing_R2.fastq
else
    MODE="SE"
    echo "Detected single-end input (R1 only)..."
    cutadapt \
      -a ACCGATTGTGAAATGAGCCAX \
      -e 0.1 -m 21 \
      -o cutadapt_files/trimmed_R1.fastq \
      sequencing/sequencing_R1.fastq
fi

# Find the designed oligo library.
if [ -f dna-fountain/test_file.tar.gz.dna_order.txt ]; then
    ORDER_FILE="dna-fountain/test_file.tar.gz.dna_order.txt"
elif [ -f test_files/test_file.tar.gz.dna_order.txt ]; then
    ORDER_FILE="test_files/test_file.tar.gz.dna_order.txt"
else
    echo "Error: test_file.tar.gz.dna_order.txt  not found in ./dna-fountain or project root"
    exit 1
fi

echo "Step 2: Building BWA reference from ${ORDER_FILE}..."
awk '{print ">oligo_" NR "\n" $0}' "$ORDER_FILE" > bwa_files/oligo_ref.fa
"${BWA_BIN}" index bwa_files/oligo_ref.fa

echo "Step 3: Aligning trimmed reads with ${BWA_BIN}..."

if [ "$MODE" = "PE" ]; then
    echo "Running paired-end alignment..."
    if [ "$BWA_BIN" = "bwa-mem2" ]; then
        bwa-mem2 mem -t 4 bwa_files/oligo_ref.fa \
            cutadapt_files/trimmed_R1.fastq cutadapt_files/trimmed_R2.fastq \
            > bwa_files/aln.sam
    else
        bwa mem -t 4 bwa_files/oligo_ref.fa \
            cutadapt_files/trimmed_R1.fastq cutadapt_files/trimmed_R2.fastq \
            > bwa_files/aln.sam
    fi
else
    echo "Running single-end alignment..."
    if [ "$BWA_BIN" = "bwa-mem2" ]; then
        bwa-mem2 mem -t 4 bwa_files/oligo_ref.fa \
            cutadapt_files/trimmed_R1.fastq \
            > bwa_files/aln.sam
    else
        bwa mem -t 4 bwa_files/oligo_ref.fa \
            cutadapt_files/trimmed_R1.fastq \
            > bwa_files/aln.sam
    fi
fi

echo "Step 4: Converting BWA alignments into decoder_input.seq..."
python - <<'PY'
from collections import Counter

PRIMER_F_LEN = 20
PRIMER_R_LEN = 20

refs = {}
name = None
with open('bwa_files/oligo_ref.fa') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            name = line[1:]
        else:
            refs[name] = line

counts = Counter()
with open('bwa_files/aln.sam') as f:
    for line in f:
        if line.startswith('@'):
            continue
        fields = line.rstrip('\n').split('\t')
        flag = int(fields[1])
        rname = fields[2]
        if rname == '*':
            continue
        if flag & 4:      # unmapped
            continue
        if flag & 256:    # secondary alignment
            continue
        if flag & 2048:   # supplementary alignment
            continue
        counts[rname] += 1

with open('bwa_files/decoder_input.seq', 'w') as out:
    kept = 0
    for rname, _ in counts.most_common():
        seq = refs[rname]
        payload = seq[PRIMER_F_LEN:len(seq)-PRIMER_R_LEN]
        if payload and 'N' not in payload:
            out.write(payload + '\n')
            kept += 1

print(f'Wrote {kept} sequences to bwa_files/decoder_input.seq')
PY

echo "Step 5: Switching to DNA Fountain environment..."
conda deactivate
source venv_sim/bin/activate

echo "Step 6: Moving to dna-fountain directory..."
cd dna-fountain/

echo "Step 7: Calculating decoder segment count (-n) from test_file.tar.gz..."
N_VALUE=$(python - <<'PY'
import os, math
print(math.ceil(os.path.getsize('../test_files/test_file.tar.gz') / 16))
PY
)
echo "Calculated -n value: ${N_VALUE}"

echo "Step 8: Running DNA Fountain Decoder..."
python decode.py -f ../bwa_files/decoder_input.seq -n "${N_VALUE}" -d 4 --size 16 -m 3 --gc 0.05 --rs 5 --delta 0.05 --c_dist 0.1 --out ../test_files/decoded_test_file.tar.gz

echo "--------------------------------------"
echo "PIPELINE COMPLETE: Results are in dna-fountain/decoded_test_file.tar.gz"
