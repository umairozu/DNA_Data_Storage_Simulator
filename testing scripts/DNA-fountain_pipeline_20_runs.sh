#!/bin/bash
set -euo pipefail


REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$REPO_DIR")"
cd "$ROOT_DIR"

PY="$ROOT_DIR/venv_sim/bin/python"

INPUT="$REPO_DIR/input_files/dna-fountain-input-files.tar.gz"
ENCODED="$REPO_DIR/encoded_files/dna-fountain-input-files.tar.gz.dna"
ORDER="${ENCODED}.dna_order.txt"

CUT="$REPO_DIR/cutadapt_files"
PEAR="$REPO_DIR/pear_files"
REC="$REPO_DIR/recovered_files"
LOG="$REPO_DIR/run_logs"
CSV="$REPO_DIR/run_results_20.csv"

mkdir -p "$REPO_DIR/encoded_files" "$CUT" "$PEAR" "$REC" "$LOG" \
         "$ROOT_DIR/synthesis" "$ROOT_DIR/pcr" "$ROOT_DIR/sequencing"

# ---------------- DNA Fountain parameters ----------------
CHUNK_SIZE=32
MAX_HOMOPOLYMER=3
GC=0.05
RS=2
DELTA=0.001
C_DIST=0.025
HEADER_SIZE=4
STOP_OLIGOS=72000

# 4 nt per encoded byte.
DNA_LEN=$((4 * (CHUNK_SIZE + HEADER_SIZE + RS)))

# ---------------- Simulator parameters ----------------
SYN_MUT=1

PCR_SAMPLE=3.33
PCR_CYCLES=10
PCR_POLYMERASE=1
PCR_MUT=2
PCR_CHIMERA_RATE=5
PCR_VARIANT_CAP=3

SEQ_TYPE=1
SEQ_MODE=2
SEQ_SAMPLE=0.5
SEQ_READS=750000
SEQ_READ_LEN=150
SEQ_MUT=1

RUNS=20
BASE_SEED=42


echo "=== DNA Fountain encode ==="

"$PY" "$REPO_DIR/encode.py" \
    --file_in "$INPUT" \
    --size "$CHUNK_SIZE" \
    -m "$MAX_HOMOPOLYMER" \
    --gc "$GC" \
    --rs "$RS" \
    --delta "$DELTA" \
    --c_dist "$C_DIST" \
    --stop "$STOP_OLIGOS" \
    --out "$ENCODED"

ENCODED_COUNT=$(awk 'NF && $0 !~ /^>/ {count++} END {print count+0}' "$ENCODED")
ORIGINAL_SIZE=$(stat -c '%s' "$INPUT")
CHUNK_NUM=$(( (ORIGINAL_SIZE + CHUNK_SIZE - 1) / CHUNK_SIZE ))

echo "Encoded oligos : $ENCODED_COUNT"
echo "Chunk count -n : $CHUNK_NUM"
echo "DNA length     : $DNA_LEN nt"


echo "=== Primer addition ==="

ENCODED_REL="$(realpath --relative-to="$ROOT_DIR" "$ENCODED")"

"$PY" "$ROOT_DIR/Enzyme_Addition.py" \
    --pf GTTCAGAGTTCTACAGTCCGACGATC \
    --pr TGGAATTCTCGGGTGCCAAGG \
    --gc_min 0.45 \
    --gc_max 0.60 \
    --in_file "$ENCODED_REL" \
    --out_file "$ORDER"

[[ -s "$ORDER" ]] || { echo "ERROR: order file not created"; exit 1; }

ORDER_REL="$(realpath --relative-to="$ROOT_DIR" "$ORDER")"

# Cutadapt + PEAR environment.
source ~/miniforge3/etc/profile.d/conda.sh
conda deactivate 2>/dev/null || true
conda activate dna_tools

printf 'iteration,seed,encoded_oligos,total_read_pairs,full_length_reads,decoder_sequences,pcpm,dropout_pct\n' > "$CSV"


for ((i=1; i<=RUNS; i++)); do

    seed=$((BASE_SEED + i - 1))
    id=$(printf "run_%02d" "$i")

    echo
    echo "============================================================"
    echo "$id / $RUNS   seed=$seed"
    echo "============================================================"

    SYN_NAME="fountain_${id}_synthesis_OUT.txt"
    PCR_NAME="fountain_${id}_pcr_OUT.txt"

    SYN_OUT="$ROOT_DIR/synthesis/$SYN_NAME"
    PCR_OUT="$ROOT_DIR/pcr/$PCR_NAME"

    R1="$ROOT_DIR/sequencing/fountain_${id}_R1.fastq"
    R2="$ROOT_DIR/sequencing/fountain_${id}_R2.fastq"

    CDIR="$CUT/$id"
    PDIR="$PEAR/$id"
    RECOVERED="$REC/${id}_decoded.tar.gz"

    mkdir -p "$CDIR" "$PDIR"

    # ---------------- Synthesis ----------------
    echo "[1/7] Synthesis"

    SYN_CMD=(
        "$PY" "$ROOT_DIR/synthesis.py"
        --mut "$SYN_MUT"
        --in_file "$ORDER_REL"
        --out_file "$SYN_NAME"
    )
    "${SYN_CMD[@]}"

    if [[ ! -s "$SYN_OUT" ]]; then
        echo "ERROR: synthesis.py did not create:"
        echo "  $SYN_OUT"
        exit 1
    fi

    # ---------------- PCR ----------------
    echo "[2/7] PCR"

    "$PY" "$ROOT_DIR/pcr.py" \
        --s "$PCR_SAMPLE" \
        --n "$PCR_CYCLES" \
        --polymerase "$PCR_POLYMERASE" \
        --mut "$PCR_MUT" \
        --chimera_rate "$PCR_CHIMERA_RATE" \
        --variant_cap "$PCR_VARIANT_CAP" \
        --seed "$seed" \
        --in_file "$SYN_OUT" \
        --out_file "$PCR_NAME"

    # ---------------- Sequencing ----------------
    echo "[3/7] Sequencing"

    rm -f "$ROOT_DIR/sequencing/sequencing_R1.fastq" \
          "$ROOT_DIR/sequencing/sequencing_R2.fastq"

    SEQ_CMD=(
        "$PY" "$ROOT_DIR/sequencing.py"
        --type "$SEQ_TYPE"
        --m "$SEQ_MODE"
        --s "$SEQ_SAMPLE"
        --t "$SEQ_READS"
        --rl "$SEQ_READ_LEN"
        --mut "$SEQ_MUT"
        --in_file "$PCR_OUT"
        --order_file "$ORDER"
    )
    "${SEQ_CMD[@]}"

    if [[ ! -s "$ROOT_DIR/sequencing/sequencing_R1.fastq" || \
          ! -s "$ROOT_DIR/sequencing/sequencing_R2.fastq" ]]; then
        echo "ERROR: sequencing.py did not create both paired FASTQ files."
        exit 1
    fi

    mv "$ROOT_DIR/sequencing/sequencing_R1.fastq" "$R1"
    mv "$ROOT_DIR/sequencing/sequencing_R2.fastq" "$R2"

    # ---------------- Cutadapt ----------------
    echo "[4/7] Cutadapt"

    # Matches the primers added above.
    cutadapt \
        -a TGGAATTCTCGGGTGCCAAGGX \
        -A GATCGTCGGACTGTAGAACTCTGAACX \
        -e 0.1 -m 21 \
        -o "$CDIR/trimmed_R1.fastq" \
        -p "$CDIR/trimmed_R2.fastq" \
        "$R1" "$R2"

    # ---------------- PEAR ----------------
    echo "[5/7] PEAR"

    pear \
        -f "$CDIR/trimmed_R1.fastq" \
        -r "$CDIR/trimmed_R2.fastq" \
        -o "$PDIR/merged"

    # ---------------- Decoder input ----------------
    echo "[6/7] Decoder input + PCPM/dropout"

    FULL="$PDIR/merged.full.seq"
    COUNTS="$PDIR/merged.counts"
    DECODER="$PDIR/decoder_input.seq"

    # size=32, header=4, RS=2 -> 152 nt.
    awk -v len="$DNA_LEN" '
        NR%4==2 && length($0)==len && $0 ~ /^[ACGT]+$/ {print $0}
    ' "$PDIR/merged.assembled.fastq" > "$FULL"

    sort "$FULL" | uniq -c | sed 's/^ *//' | sort -rn -k1,1 > "$COUNTS"
    awk '{print $2}' "$COUNTS" | grep -v 'N' > "$DECODER"

    TOTAL_READ_PAIRS=$(awk 'END {print int(NR/4)}' "$R1")
    FULL_LENGTH_READS=$(wc -l < "$FULL")
    DECODER_COUNT=$(wc -l < "$DECODER")

    read -r PCPM DROPOUT <<< "$(
        "$PY" - "$ENCODED" "$FULL" "$TOTAL_READ_PAIRS" <<'PY'
import sys
from collections import Counter

encoded_file, observed_file, total_reads = sys.argv[1], sys.argv[2], int(sys.argv[3])

with open(encoded_file) as f:
    encoded = [x.strip() for x in f if x.strip() and not x.startswith(">")]

with open(observed_file) as f:
    observed = Counter(x.strip() for x in f if x.strip())

perfect_counts = [observed[x] for x in encoded]
perfect_total = sum(perfect_counts)
observed_oligos = sum(x > 0 for x in perfect_counts)

pcpm = (perfect_total * 1_000_000.0) / (total_reads * len(encoded))
dropout = 100.0 * (len(encoded) - observed_oligos) / len(encoded)

print(f"{pcpm:.2f} {dropout:.2f}")
PY
    )"

    echo "Full-length reads : $FULL_LENGTH_READS"
    echo "Decoder sequences : $DECODER_COUNT"
    echo "PCPM              : $PCPM"
    echo "Dropout           : $DROPOUT%"

    # ---------------- Decode ----------------
    echo "[7/7] Decode"

    rm -f "$RECOVERED"

    if (( DECODER_COUNT > 0 )); then
        "$PY" "$REPO_DIR/decode.py" \
            -f "$DECODER" \
            -n "$CHUNK_NUM" \
            -d "$HEADER_SIZE" \
            --size "$CHUNK_SIZE" \
            -m "$MAX_HOMOPOLYMER" \
            --gc "$GC" \
            --rs "$RS" \
            --delta "$DELTA" \
            --c_dist "$C_DIST" \
            --out "$RECOVERED"
    fi


    printf '%d,%d,%d,%d,%d,%d,%d,%s,%s,%s\n' \
        "$i" "$seed" "$ENCODED_COUNT" "$SEQ_READS" \
        "$TOTAL_READ_PAIRS" "$FULL_LENGTH_READS" "$DECODER_COUNT" \
        "$PCPM" "$DROPOUT" >> "$CSV"

done

echo
echo "============================================================"
echo "DONE"
echo "Results: $CSV"
echo "============================================================"
column -s, -t "$CSV" 2>/dev/null || cat "$CSV"

