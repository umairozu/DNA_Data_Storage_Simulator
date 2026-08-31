#!/bin/bash
set -euo pipefail

# ============================================================
# DNA Fountain robustness benchmark
#   6 channel conditions x 3 replicates = 18 runs
#
# Conditions:
#   Perfect          : synthesis=0, storage=0, PCR=0, sequencing=0
#   Synthesis only   : synthesis=1, storage=0, PCR=0, sequencing=0
#   Storage only     : synthesis=0, storage=2, PCR=0, sequencing=0
#   PCR only         : synthesis=0, storage=0, PCR=2, sequencing=0
#   Sequencing only  : synthesis=0, storage=0, PCR=0, sequencing=1
#   Combined         : synthesis=1, storage=2, PCR=2, sequencing=1
#
# DNA Fountain encoding and primer addition are run once because
# they are unchanged across channel conditions and replicates.
# ============================================================

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$REPO_DIR")"
cd "$ROOT_DIR"

SCHEME="DNA Fountain"
SIM_PYTHON="$ROOT_DIR/venv_sim/bin/python"

INPUT_FILE="$REPO_DIR/input_files/Mona_Lisa.jpg"
ENCODED_FILE="$REPO_DIR/encoded_files/dna_fountain_encoded.dna"
ORDER_FILE="${ENCODED_FILE}.dna_order.txt"

CUTADAPT_DIR="$REPO_DIR/cutadapt_files"
PEAR_DIR="$REPO_DIR/pear_files"
RECOVERED_DIR="$REPO_DIR/recovered_files"
LOG_DIR="$REPO_DIR/run_logs"

RESULTS_TSV="$REPO_DIR/run_results.tsv"
RESULTS_CSV="$REPO_DIR/run_results.csv"

mkdir -p \
    "$REPO_DIR/encoded_files" \
    "$RECOVERED_DIR" \
    "$CUTADAPT_DIR" \
    "$PEAR_DIR" \
    "$LOG_DIR" \
    "$ROOT_DIR/synthesis" \
    "$ROOT_DIR/storage" \
    "$ROOT_DIR/pcr" \
    "$ROOT_DIR/sequencing"

# ------------------------------------------------------------
# DNA Fountain codec parameters.
# ------------------------------------------------------------
CHUNK_SIZE=23
MAX_HOMOPOLYMER=3
GC=0.05
RS=2
DELTA=0.001
C_DIST=0.025
ALPHA=0.25
HEADER_SIZE=4
MAX_HAMMING=0

# ------------------------------------------------------------
# Fixed simulator/channel settings.
# ------------------------------------------------------------
TEMPS=(10)
PH=7
RH=50
WEEKS_LIST=(0)
ENCAP=1

STORAGE_DAMAGE_MODEL="depurination"
STORAGE_DAMAGE_FATE="variants"
STORAGE_ENCAP_DAMAGE_FACTOR="bulk_ratio"
STORAGE_VARIANT_CAP=8
STORAGE_SURVIVAL_MODE="binomial"

PCR_SAMPLE=1
PCR_CYCLES=30
PCR_POLYMERASE=1
PCR_CHIMERA_RATE=5
PCR_VARIANT_CAP=3

SEQ_TYPE=1
SEQ_MODE=2
SEQ_SAMPLE=1
SEQ_COVERAGES=(1 2 5 10)
SEQ_READ_LEN=1000

NUM_REPLICATES=3
BASE_SEED=42


CONDITIONS=(
    #"Perfect|perfect|0|0|0|0"
    #"Synthesis only|synthesis_only|1|0|0|0"
    #"Storage only|storage_only|0|2|0|0"
    #"PCR only|pcr_only|0|0|2|0"
    #"Sequencing only|sequencing_only|0|0|0|1"
    #"Combined|combined|1|2|2|1"
    "Sequencing coverage|sequencing_coverage|0|0|0|0"
)

# ------------------------------------------------------------
# Basic checks
# ------------------------------------------------------------
for required in \
    "$SIM_PYTHON" \
    "$INPUT_FILE" \
    "$REPO_DIR/encode.py" \
    "$REPO_DIR/decode.py" \
    "$ROOT_DIR/Enzyme_Addition.py" \
    "$ROOT_DIR/synthesis.py" \
    "$ROOT_DIR/storage.py" \
    "$ROOT_DIR/pcr.py" \
    "$ROOT_DIR/sequencing.py" \
    "$ROOT_DIR/RawData.xlsx"
do
    if [[ ! -e "$required" ]]; then
        echo "ERROR: Required path not found: $required" >&2
        exit 1
    fi
done

if (( NUM_REPLICATES < 1 )); then
    echo "ERROR: NUM_REPLICATES must be at least 1." >&2
    exit 1
fi

if (( ${#SEQ_COVERAGES[@]} < 1 )); then
    echo "ERROR: SEQ_COVERAGES must contain at least one value." >&2
    exit 1
fi

for cov in "${SEQ_COVERAGES[@]}"; do
    if (( cov < 1 )); then
        echo "ERROR: All SEQ_COVERAGES values must be at least 1." >&2
        exit 1
    fi
done


echo "============================================================"
echo "DNA Fountain sequencing-coverage benchmark"
echo "Conditions      : ${#CONDITIONS[@]}"
echo "Coverage levels : ${#SEQ_COVERAGES[@]} (${SEQ_COVERAGES[*]}x)"
echo "Replicates      : $NUM_REPLICATES"
echo "Total runs      : $((${#CONDITIONS[@]} * ${#SEQ_COVERAGES[@]} * NUM_REPLICATES))"
echo "Coverage sweep  : ${SEQ_COVERAGES[*]} read pairs / encoded oligo"
echo "Seeds      : $BASE_SEED .. $((BASE_SEED + NUM_REPLICATES - 1))"
echo "============================================================"


echo
echo "Step 1: DNA Fountain encoding..."

"$SIM_PYTHON" "$REPO_DIR/encode.py" \
    --file_in "$INPUT_FILE" \
    --size "$CHUNK_SIZE" \
    -m "$MAX_HOMOPOLYMER" \
    --gc "$GC" \
    --rs "$RS" \
    --delta "$DELTA" \
    --c_dist "$C_DIST" \
    --alpha "$ALPHA" \
    --out "$ENCODED_FILE"

ENCODED_COUNT=$(awk 'NF && $0 !~ /^>/ {count++} END {print count+0}' "$ENCODED_FILE")
if (( ENCODED_COUNT == 0 )); then
    echo "ERROR: DNA Fountain encoder produced no oligos." >&2
    exit 1
fi

ORIGINAL_SIZE=$(stat -c '%s' "$INPUT_FILE")
CHUNK_NUM=$(( (ORIGINAL_SIZE + CHUNK_SIZE - 1) / CHUNK_SIZE ))
DNA_LEN=$((4 * (CHUNK_SIZE + HEADER_SIZE + RS)))

echo "Encoded oligos      : $ENCODED_COUNT"
echo "Input chunks        : $CHUNK_NUM"
echo "Expected oligo length: $DNA_LEN nt"
echo "Coverage sweep      : ${SEQ_COVERAGES[*]}x"


echo
echo "Step 2: Adding primers..."

ENCODED_REL="$(realpath --relative-to="$ROOT_DIR" "$ENCODED_FILE")"
ORDER_REL="$(realpath --relative-to="$ROOT_DIR" "$ORDER_FILE")"

"$SIM_PYTHON" "$ROOT_DIR/Enzyme_Addition.py" \
    --pf GTTCAGAGTTCTACAGTCCGACGATC \
    --pr TGGAATTCTCGGGTGCCAAGG \
    --gc_min 0.45 \
    --gc_max 0.60 \
    --in_file "$ENCODED_REL" \
    --out_file "$ORDER_FILE"

# Activate the environment containing Cutadapt and PEAR once.
source ~/miniforge3/etc/profile.d/conda.sh
conda deactivate 2>/dev/null || true
conda activate dna_tools

# ============================================================
# Results header
# ============================================================
printf 'scheme\tcondition\treplicate\tseed\tsynthesis_mut\tstorage_mut\tpcr_mut\tsequencing_mut\ttemp_c\tph\trh_pct\tweeks\tencap\tbulk_remaining_pct\tdamaged_molecule_pct\tpcr_cycles\tseq_coverage\trequested_read_pairs\ttotal_read_pairs\tencoded_oligos\tdecode_success\tpcpm\tdropout_pct\tdecoder_input_count\tdecode_status\n' \
> "$RESULTS_TSV"

run_one() {
    local condition_name="$1"
    local condition_safe="$2"
    local replicate="$3"
    local syn_mut="$4"
    local storage_mut="$5"
    local pcr_mut="$6"
    local seq_mut="$7"
    local temp="$8"
    local weeks="$9"
    local seq_coverage="${10}"

    local seq_reads=$((ENCODED_COUNT * seq_coverage))
    local seed=$((BASE_SEED + replicate - 1))
    local temp_safe="${temp//./p}"
    local weeks_safe="${weeks//./p}"
    local run_id="dna_fountain_${condition_safe}_T${temp_safe}_W${weeks_safe}_COV${seq_coverage}_rep${replicate}"
    local run_log="$LOG_DIR/${run_id}.log"
    local storage_stage_log="$LOG_DIR/${run_id}_storage_stage.log"

    local synth_basename="${run_id}_synthesis_OUT.txt"
    local storage_basename="${run_id}_storage_OUT.txt"
    local pcr_basename="${run_id}_pcr_OUT.txt"

    local synth_out="$ROOT_DIR/synthesis/$synth_basename"
    local storage_out="$ROOT_DIR/storage/$storage_basename"
    local pcr_out="$ROOT_DIR/pcr/$pcr_basename"

    local seq_r1_default="$ROOT_DIR/sequencing/sequencing_R1.fastq"
    local seq_r2_default="$ROOT_DIR/sequencing/sequencing_R2.fastq"
    local seq_r1="$ROOT_DIR/sequencing/${run_id}_R1.fastq"
    local seq_r2="$ROOT_DIR/sequencing/${run_id}_R2.fastq"
    local seq_marker="$ROOT_DIR/sequencing/.${run_id}_before_sequencing"

    local trimmed_r1="$CUTADAPT_DIR/${run_id}_trimmed_R1.fastq"
    local trimmed_r2="$CUTADAPT_DIR/${run_id}_trimmed_R2.fastq"
    local pear_prefix="$PEAR_DIR/${run_id}_merged"
    local full_seq="$PEAR_DIR/${run_id}_merged.full.seq"
    local counts_file="$PEAR_DIR/${run_id}_merged.counts"
    local decoder_input="$PEAR_DIR/${run_id}_decoder_input.seq"
    local recovered_file="$RECOVERED_DIR/${run_id}_recovered.jpg"

    {
        echo
        echo "============================================================"
        echo "RUN: $run_id"
        echo "Condition: $condition_name"
        echo "Replicate: $replicate/$NUM_REPLICATES"
        echo "Seed: $seed"
        echo "Storage: ${temp} C, ${weeks} week(s), RH=${RH}%, encap=${ENCAP}"
        echo "Sequencing coverage: ${seq_coverage}x (${seq_reads} read pairs)"
        echo "mut = synthesis:$syn_mut storage:$storage_mut PCR:$pcr_mut sequencing:$seq_mut"
        echo "============================================================"

        # ---------------- Synthesis ----------------
        echo "[1/8] Synthesis"
        synth_cmd=(
            "$SIM_PYTHON" "$ROOT_DIR/synthesis.py"
            --mut "$syn_mut"
            --in_file "$ORDER_REL"
            --out_file "$synth_basename"
        )
 	"${synth_cmd[@]}"

        # ---------------- Storage ----------------
        echo "[2/8] Storage"
        "$SIM_PYTHON" "$ROOT_DIR/storage.py" \
            --temp "$temp" \
            --ph "$PH" \
            --rh "$RH" \
            --week "$weeks" \
            --encap "$ENCAP" \
            --mut "$storage_mut" \
            --damage_model "$STORAGE_DAMAGE_MODEL" \
            --damage_fate "$STORAGE_DAMAGE_FATE" \
            --encap_damage_factor "$STORAGE_ENCAP_DAMAGE_FACTOR" \
            --variant_cap "$STORAGE_VARIANT_CAP" \
            --survival_mode "$STORAGE_SURVIVAL_MODE" \
            --seed "$seed" \
            --arrhenius_xlsx "$ROOT_DIR/RawData.xlsx" \
            --in_file "$synth_out" \
            --out_file "$storage_basename" \
            | tee "$storage_stage_log"

        local bulk_remaining_pct damaged_molecule_pct

        bulk_remaining_pct="$(
            awk -F: '/Bulk remaining fraction:/ {
                value=$2
                gsub(/[%[:space:]]/, "", value)
                print value
                exit
            }' "$storage_stage_log"
        )"

        damaged_molecule_pct="$(
            awk -F: '/Damaged fraction of survivors:/ {
                value=$2
                gsub(/[%[:space:]]/, "", value)
                print value
                exit
            }' "$storage_stage_log"
        )"

        if [[ -z "$bulk_remaining_pct" ]]; then
            echo "ERROR: Could not parse Bulk remaining fraction from storage.py output." >&2
            exit 1
        fi

        if [[ -z "$damaged_molecule_pct" ]]; then
            damaged_molecule_pct="0.00000000"
        fi

        echo "Bulk remaining    : ${bulk_remaining_pct}%"
        echo "Damaged molecules : ${damaged_molecule_pct}%"

        # ---------------- PCR ----------------
        echo "[3/8] PCR"
        "$SIM_PYTHON" "$ROOT_DIR/pcr.py" \
            --s "$PCR_SAMPLE" \
            --n "$PCR_CYCLES" \
            --polymerase "$PCR_POLYMERASE" \
            --mut "$pcr_mut" \
            --chimera_rate "$PCR_CHIMERA_RATE" \
            --variant_cap "$PCR_VARIANT_CAP" \
            --seed "$seed" \
            --in_file "$storage_out" \
            --out_file "$pcr_basename"

        # ---------------- Sequencing ----------------
        echo "[4/8] Sequencing"

        rm -f "$seq_r1_default" "$seq_r2_default"
        rm -f "$ROOT_DIR/sequencing/sequencing_R1.fq" \
              "$ROOT_DIR/sequencing/sequencing_R2.fq"

        rm -f "$seq_marker"
        touch "$seq_marker"

        seq_cmd=(
            "$SIM_PYTHON" "$ROOT_DIR/sequencing.py"
            --type "$SEQ_TYPE"
            --m "$SEQ_MODE"
            --s "$SEQ_SAMPLE"
            --t "$seq_reads"
            --rl "$SEQ_READ_LEN"
            --mut "$seq_mut"
            --in_file "$pcr_out"
            --order_file "$ORDER_FILE"
        )
	"${seq_cmd[@]}"


        local produced_r1=""
        local produced_r2=""

        if [[ -s "$seq_r1_default" ]]; then
            produced_r1="$seq_r1_default"
        fi
        if [[ -s "$seq_r2_default" ]]; then
            produced_r2="$seq_r2_default"
        fi

        if [[ -z "$produced_r1" ]]; then
            produced_r1="$(
                find "$ROOT_DIR" -maxdepth 3 -type f -newer "$seq_marker" \
                    \( -iname '*R1*.fastq' -o -iname '*R1*.fq' \) \
                    -printf '%T@\t%p\n' 2>/dev/null \
                | sort -nr \
                | head -n 1 \
                | cut -f2-
            )"
        fi

        if [[ -z "$produced_r2" ]]; then
            produced_r2="$(
                find "$ROOT_DIR" -maxdepth 3 -type f -newer "$seq_marker" \
                    \( -iname '*R2*.fastq' -o -iname '*R2*.fq' \) \
                    -printf '%T@\t%p\n' 2>/dev/null \
                | sort -nr \
                | head -n 1 \
                | cut -f2-
            )"
        fi

        if [[ -z "$produced_r1" || -z "$produced_r2" || \
              ! -s "$produced_r1" || ! -s "$produced_r2" ]]; then
            echo "ERROR: sequencing.py completed, but both paired FASTQ files could not be located." >&2
            echo "Files created/updated by sequencing.py:" >&2
            find "$ROOT_DIR" -maxdepth 3 -type f -newer "$seq_marker" \
                -printf '  %p\n' 2>/dev/null | sort >&2 || true
            rm -f "$seq_marker"
            exit 1
        fi

        echo "Sequencing R1 source: $produced_r1"
        echo "Sequencing R2 source: $produced_r2"

        mv -f "$produced_r1" "$seq_r1"
        mv -f "$produced_r2" "$seq_r2"
        rm -f "$seq_marker"

        if [[ ! -s "$seq_r1" || ! -s "$seq_r2" ]]; then
            echo "ERROR: Failed to preserve paired FASTQ files for $run_id" >&2
            exit 1
        fi

        # ---------------- Cutadapt ----------------
        echo "[5/8] Cutadapt"
        cutadapt \
            -a TGGAATTCTCGGGTGCCAAGGX \
            -A GATCGTCGGACTGTAGAACTCTGAACX \
            -e 0.1 \
            -m 21 \
            -o "$trimmed_r1" \
            -p "$trimmed_r2" \
            "$seq_r1" \
            "$seq_r2"

        # ---------------- PEAR ----------------
        echo "[6/8] PEAR"
        rm -f "${pear_prefix}.assembled.fastq" \
              "${pear_prefix}.discarded.fastq" \
              "${pear_prefix}.unassembled.forward.fastq" \
              "${pear_prefix}.unassembled.reverse.fastq"

        pear \
            -f "$trimmed_r1" \
            -r "$trimmed_r2" \
            -o "$pear_prefix"

        # ---------------- Decoder input + metrics ----------------
        echo "[7/8] Decoder input + metrics"

        #   1. retain valid, full-length assembled sequences;
        #   2. count identical sequences;
        #   3. pass one copy of each unique sequence to decode.py.
        awk -v len="$DNA_LEN" '
            NR%4==2 &&
            length($0)==len &&
            $0 ~ /^[ACGT]+$/ {
                print $0
            }
        ' "${pear_prefix}.assembled.fastq" > "$full_seq"

        local full_length_reads
        full_length_reads=$(wc -l < "$full_seq")

        sort "$full_seq" \
            | uniq -c \
            | sed 's/^ *//' \
            | sort -rn -k1,1 \
            > "$counts_file"

        awk '$2 !~ /N/ {print $2}' "$counts_file" > "$decoder_input"

        local decoder_count
        local total_read_pairs
        decoder_count=$(wc -l < "$decoder_input")
        total_read_pairs=$(awk 'END {print int(NR/4)}' "$seq_r1")

        local pcpm dropout_pct
        read -r pcpm dropout_pct <<< "$(
            "$SIM_PYTHON" - "$ENCODED_FILE" "$full_seq" "$total_read_pairs" <<'PY'
import sys
from collections import Counter

encoded_file, observed_file, total_reads = sys.argv[1], sys.argv[2], int(sys.argv[3])

with open(encoded_file) as f:
    encoded = [line.strip() for line in f if line.strip() and not line.startswith(">")]

with open(observed_file) as f:
    observed = Counter(line.strip() for line in f if line.strip())

if not encoded or total_reads <= 0:
    print("0.000000 100.000000")
    raise SystemExit

perfect_counts = [observed[oligo] for oligo in encoded]
perfect_total = sum(perfect_counts)
observed_oligos = sum(count > 0 for count in perfect_counts)

# perfect-call abundance normalized by total read-pair depth and
# by the number of encoded oligos.
pcpm = (perfect_total * 1_000_000.0) / (total_reads * len(encoded))
dropout_pct = 100.0 * (len(encoded) - observed_oligos) / len(encoded)

print(f"{pcpm:.6f} {dropout_pct:.6f}")
PY
        )"

        echo "Full-length reads   : $full_length_reads"
        echo "Decoder input count : $decoder_count"
        echo "PCPM                : $pcpm"
        echo "Dropout              : $dropout_pct%"

        # ---------------- Decode ----------------
        echo "[8/8] DNA Fountain decode"
        local decode_status="DECODE_FAILED"
        local decode_success="FAILURE"
        rm -f "$recovered_file"

        if (( decoder_count > 0 )); then
            if "$SIM_PYTHON" "$REPO_DIR/decode.py" \
                -f "$decoder_input" \
                --header_size "$HEADER_SIZE" \
                --rs "$RS" \
                --delta "$DELTA" \
                --c_dist "$C_DIST" \
                -n "$CHUNK_NUM" \
                -m "$MAX_HOMOPOLYMER" \
                --gc "$GC" \
                --max_hamming "$MAX_HAMMING" \
                --size "$CHUNK_SIZE" \
                --out "$recovered_file"
            then
                if [[ -f "$recovered_file" ]]; then
                    truncate -s "$ORIGINAL_SIZE" "$recovered_file"

                    if cmp -s "$INPUT_FILE" "$recovered_file"; then
                        decode_status="SUCCESS"
                        decode_success="SUCCESS"
                    else
                        decode_status="DIFFERENT"
                    fi
                fi
            fi
        else
            decode_status="NO_DECODER_INPUT"
        fi

        echo "Decode success: $decode_success"
        echo "Decode status : $decode_status"

        printf '%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%s\n' \
            "$SCHEME" \
            "$condition_name" \
            "$replicate" \
            "$seed" \
            "$syn_mut" \
            "$storage_mut" \
            "$pcr_mut" \
            "$seq_mut" \
            "$temp" \
            "$PH" \
            "$RH" \
            "$weeks" \
            "$ENCAP" \
            "$bulk_remaining_pct" \
            "$damaged_molecule_pct" \
            "$PCR_CYCLES" \
            "$seq_coverage" \
            "$seq_reads" \
            "$total_read_pairs" \
            "$ENCODED_COUNT" \
            "$decode_success" \
            "$pcpm" \
            "$dropout_pct" \
            "$decoder_count" \
            "$decode_status" \
            >> "$RESULTS_TSV"

    } 2>&1 | tee "$run_log"
}

# ============================================================
# Run all conditions and replicates sequentially
# ============================================================
for temp in "${TEMPS[@]}"; do
    for weeks in "${WEEKS_LIST[@]}"; do
        for condition in "${CONDITIONS[@]}"; do
            IFS='|' read -r condition_name condition_safe syn_mut storage_mut pcr_mut seq_mut <<< "$condition"

            for seq_coverage in "${SEQ_COVERAGES[@]}"; do
                for ((replicate=1; replicate<=NUM_REPLICATES; replicate++)); do
                    run_one \
                        "$condition_name" \
                        "$condition_safe" \
                        "$replicate" \
                        "$syn_mut" \
                        "$storage_mut" \
                        "$pcr_mut" \
                        "$seq_mut" \
                        "$temp" \
                        "$weeks" \
                        "$seq_coverage"
                done
            done
        done
    done
done

# TSV is the primary results file; also create a CSV copy.
"$SIM_PYTHON" - "$RESULTS_TSV" "$RESULTS_CSV" <<'PY'
import csv
import sys

src, dst = sys.argv[1], sys.argv[2]
with open(src, newline="") as fin, open(dst, "w", newline="") as fout:
    reader = csv.reader(fin, delimiter="\t")
    writer = csv.writer(fout)
    writer.writerows(reader)
PY

echo
echo "============================================================"
echo "BENCHMARK COMPLETE"
echo "TSV results: $RESULTS_TSV"
echo "CSV results: $RESULTS_CSV"
echo "Run logs   : $LOG_DIR"
echo "============================================================"

