#!/bin/bash

set -euo pipefail

mkdir -p cutadapt_files pear_files

echo "Starting Bioinformatics Pipeline..."

deactivate 2>/dev/null || true
conda deactivate 2>/dev/null || true

echo "Activating conda environment: dna_tools..."
source ~/miniforge3/etc/profile.d/conda.sh
conda activate dna_tools

echo "Step 1: Trimming adapters with cutadapt..."
cutadapt -a ACCGATTGTGAAATGAGCCAX -A TTGCACGGCAGGTCATTTATX -e 0.1 -m 21 -o cutadapt_files/trimmed_R1.fastq -p cutadapt_files/trimmed_R2.fastq sequencing/sequencing_R1.fastq sequencing/sequencing_R2.fastq

echo "Step 2: Merging paired-end reads with PEAR..."
pear -f cutadapt_files/trimmed_R1.fastq -r cutadapt_files/trimmed_R2.fastq -o pear_files/merged

echo "Step 3: Moving to pear_files directory..."
cd pear_files

echo "Step 4: Filtering for 100bp sequences..."
awk '(NR%4==2 && length($0)==100)' merged.assembled.fastq > merged.full.seq

echo "Step 5: Counting and sorting unique sequences..."
sort merged.full.seq | uniq -c | sed 's/^ *//' | sort -rn -k1,1 > merged.counts

echo "Step 6: Cleaning up sequences for decoder..."
awk '{print $2}' merged.counts | grep -v 'N' > decoder_input.seq

echo "Step 7: Switching to DNA Fountain environment..."
cd ..
conda deactivate

cd "$(dirname "$0")"

source venv_sim/bin/activate

echo "Step 8: Moving to dna-fountain directory..."
cd dna-fountain/

echo "Step 9: Calculating decoder segment count (-n) from test_file.tar.gz..."
N_VALUE=$(python - <<'PY'
import os, math
print(math.ceil(os.path.getsize("../test_files/test_file.tar.gz") / 16))
PY
)
echo "Calculated -n value: ${N_VALUE}"

echo "Step 10: Running DNA Fountain Decoder..."
python decode.py -f ../pear_files/decoder_input.seq -n "${N_VALUE}" -d 4 --size 16 -m 3 --gc 0.05 --rs 5 --delta 0.05 --out ../test_files/decoded_test_file.tar.gz
