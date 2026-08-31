#!/bin/bash
set -euo pipefail

# ============================================================
# Goldman test script
# ============================================================

DNA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$DNA_DIR")"
cd "$ROOT_DIR"

SCHEME="Goldman"
SIM_PYTHON="$ROOT_DIR/venv_sim/bin/python"

# Original Goldman codec integration.
GOLDMAN_REPO="${GOLDMAN_REPO:-$HOME/Desktop/CS 401_402/DNA_Data_Storage_Simulator/DNA}"
GOLDMAN_PYTHON="${GOLDMAN_PYTHON:-$HOME/miniforge3/envs/goldman_py2/bin/python2.7}"
GOLDMAN_DNA_SCRIPT="$GOLDMAN_REPO/dna/dna.py"

goldman_encode() {
    local input="$1"
    local output="$2"

    input="$(realpath "$input")"

    local output_dir
    output_dir="$(dirname "$output")"
    mkdir -p "$output_dir"
    output="$(realpath "$output_dir")/$(basename "$output")"

    echo "============================================"
    echo "Goldman encoding"
    echo "============================================"
    echo "Input : $input"
    echo "Output: $output"
    echo

    "$GOLDMAN_PYTHON" "$GOLDMAN_DNA_SCRIPT" -s "$input"

    local zip_file="${input}.splitted.zip"

    if [[ ! -f "$zip_file" ]]; then
        echo "ERROR: Goldman encoder did not create:" >&2
        echo "  $zip_file" >&2
        return 1
    fi

    python3 - "$zip_file" "$output" <<'PY'
import re
import sys
import zipfile

zip_path = sys.argv[1]
output_path = sys.argv[2]

def natural_key(name):
    parts = re.split(r"(\d+)", name)
    return [int(p) if p.isdigit() else p for p in parts]

with zipfile.ZipFile(zip_path, "r") as zf:
    names = sorted(zf.namelist(), key=natural_key)
    sequences = []

    for name in names:
        sequence = zf.read(name).decode("ascii").strip().upper()
        if sequence:
            sequences.append(sequence)

with open(output_path, "w") as out:
    for sequence in sequences:
        out.write(sequence + "\n")

print("Number of Goldman oligos:", len(sequences))

if sequences:
    lengths = [len(x) for x in sequences]
    print("Minimum oligo length:", min(lengths))
    print("Maximum oligo length:", max(lengths))
PY

    echo
    echo "Goldman encoding finished."
    echo "DNA oligos:"
    echo "  $output"
}

goldman_decode() {
    local input="$1"
    local output="$2"

    input="$(realpath "$input")"

    local output_dir
    output_dir="$(dirname "$output")"
    mkdir -p "$output_dir"
    output="$(realpath "$output_dir")/$(basename "$output")"

    echo "============================================"
    echo "Goldman decoding"
    echo "============================================"
    echo "Input : $input"
    echo "Output: $output"
    echo

    local tmpdir
    tmpdir="$(mktemp -d)"

    local zip_file="$tmpdir/goldman.splitted.zip"
    local decoded_file="$tmpdir/goldman.decoded"

    if ! python3 - "$input" "$zip_file" <<'PY'
import sys
import zipfile

input_path = sys.argv[1]
zip_path = sys.argv[2]

sequences = []

with open(input_path, "r") as f:
    for line in f:
        sequence = line.strip().upper()

        if not sequence:
            continue

        invalid = set(sequence) - set("ACGT")

        if invalid:
            raise ValueError(
                "Invalid DNA symbols detected: %s"
                % sorted(invalid)
            )

        sequences.append(sequence)

width = max(1, len(str(len(sequences))))

with zipfile.ZipFile(
        zip_path,
        "w",
        compression=zipfile.ZIP_DEFLATED
    ) as zf:
    for i, sequence in enumerate(sequences):
        filename = "fragment.%0*d" % (width, i)
        zf.writestr(filename, sequence)

print("Sequences supplied to Goldman decoder:", len(sequences))
PY
    then
        rm -rf "$tmpdir"
        return 1
    fi

    if ! "$GOLDMAN_PYTHON" "$GOLDMAN_DNA_SCRIPT" -j "$zip_file"; then
        rm -rf "$tmpdir"
        return 1
    fi

    if [[ ! -f "$decoded_file" ]]; then
        echo "ERROR: Goldman decoder failed." >&2
        rm -rf "$tmpdir"
        return 1
    fi

    cp "$decoded_file" "$output"
    rm -rf "$tmpdir"

    echo
    echo "Goldman decoding finished."
    echo "Recovered file:"
    echo "  $output"
}

INPUT_FILE="$DNA_DIR/input_files/Mona_Lisa.jpg"
ENCODED_FILE="$DNA_DIR/encoded_files/goldman_encoded.dna"
ORDER_FILE="${ENCODED_FILE}.dna_order.txt"

CUTADAPT_DIR="$DNA_DIR/cutadapt_files"
PEAR_DIR="$DNA_DIR/pear_files"
RECOVERED_DIR="$DNA_DIR/recovered_files"
LOG_DIR="$DNA_DIR/run_logs"

RESULTS_TSV="$DNA_DIR/run_results.tsv"
RESULTS_CSV="$DNA_DIR/run_results.csv"

mkdir -p \
    "$DNA_DIR/encoded_files" \
    "$RECOVERED_DIR" \
    "$CUTADAPT_DIR" \
    "$PEAR_DIR" \
    "$LOG_DIR" \
    "$ROOT_DIR/synthesis" \
    "$ROOT_DIR/storage" \
    "$ROOT_DIR/pcr" \
    "$ROOT_DIR/sequencing"

# ------------------------------------------------------------
# Fixed simulator/channel settings.
# Only the four MUT selectors change between benchmark conditions.
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
SEQ_READ_LEN=150

# Goldman reconstruction: top abundant valid variants used per embedded index.
GOLDMAN_CONSENSUS_TOP_K=16

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
    "$GOLDMAN_PYTHON" \
    "$GOLDMAN_DNA_SCRIPT" \
    "$INPUT_FILE" \
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
    echo "ERROR: SEQ_COVERAGES must contain at least one coverage value." >&2
    exit 1
fi

for cov in "${SEQ_COVERAGES[@]}"; do
    if (( cov < 1 )); then
        echo "ERROR: All SEQ_COVERAGES values must be at least 1." >&2
        exit 1
    fi
done

echo "============================================================"
echo "Goldman sequencing-coverage benchmark"
echo "Conditions      : ${#CONDITIONS[@]}"
echo "Coverage levels : ${#SEQ_COVERAGES[@]} (${SEQ_COVERAGES[*]}x)"
echo "Replicates      : $NUM_REPLICATES"
echo "Total runs      : $((${#CONDITIONS[@]} * ${#SEQ_COVERAGES[@]} * NUM_REPLICATES))"
echo "Coverage sweep  : ${SEQ_COVERAGES[*]} read pairs / encoded oligo"
echo "Seeds      : $BASE_SEED .. $((BASE_SEED + NUM_REPLICATES - 1))"
echo "============================================================"


echo
echo "Step 1: Goldman encoding..."
goldman_encode "$INPUT_FILE" "$ENCODED_FILE"

ENCODED_COUNT=$(awk 'NF {count++} END {print count+0}' "$ENCODED_FILE")
if (( ENCODED_COUNT == 0 )); then
    echo "ERROR: Goldman encoder produced no oligos." >&2
    exit 1
fi

ENCODED_LENGTHS=$(
    awk 'NF {print length($0)}' "$ENCODED_FILE" \
        | sort -nu \
        | paste -sd, -
)

echo "Encoded oligos   : $ENCODED_COUNT"
echo "Payload length(s): $ENCODED_LENGTHS nt"
echo "Coverage sweep   : ${SEQ_COVERAGES[*]}x"


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
    local run_id="goldman_${condition_safe}_T${temp_safe}_W${weeks_safe}_COV${seq_coverage}_rep${replicate}"
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
    local recon_summary="$PEAR_DIR/${run_id}_goldman_reconstruction_summary.tsv"
    local index_support="$PEAR_DIR/${run_id}_goldman_index_support.tsv"
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
                v=$2; gsub(/[%[:space:]]/, "", v); print v; exit
            }' "$storage_stage_log"
        )"

        damaged_molecule_pct="$(
            awk -F: '/Damaged fraction of survivors:/ {
                v=$2; gsub(/[%[:space:]]/, "", v); print v; exit
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

        # ---------------- Goldman reconstruction + metrics ----------------
        echo "[7/8] Goldman reconstruction + metrics"

        "$SIM_PYTHON" - \
            "${pear_prefix}.assembled.fastq" \
            "$full_seq" \
            "$counts_file" \
            "$decoder_input" \
            "$recon_summary" \
            "$index_support" \
            "$ENCODED_COUNT" \
            "$GOLDMAN_CONSENSUS_TOP_K" <<'PY'
import heapq
import sys
from collections import Counter, defaultdict

(
    fastq_file,
    full_out,
    counts_out,
    decoder_out,
    summary_out,
    support_out,
    expected_count_s,
    top_k_s,
) = sys.argv[1:9]

EXPECTED_COUNT = int(expected_count_s)
TOP_K = int(top_k_s)
EXPECTED_LEN = 117
EXPECTED_ID = "00"
BASES = "ACGT"
BASE_TO_INT = {b: i for i, b in enumerate(BASES)}

DNA_INV = {
    "A": {"T": "0", "G": "1", "C": "2"},
    "C": {"A": "0", "T": "1", "G": "2"},
    "G": {"C": "0", "A": "1", "T": "2"},
    "T": {"G": "0", "C": "1", "A": "2"},
}

DNA_TABLE = {
    "A": ["C", "G", "T"],
    "C": ["G", "T", "A"],
    "G": ["T", "A", "C"],
    "T": ["A", "C", "G"],
}

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def base3_12(value):
    if value < 0:
        raise ValueError("negative Goldman index")
    if value == 0:
        s = "0"
    else:
        digits = []
        n = value
        while n:
            digits.append(str(n % 3))
            n //= 3
        s = "".join(reversed(digits))
    if len(s) > 12:
        raise ValueError("Goldman index exceeds 12 trits: %d" % value)
    return s.zfill(12)


def checksum(ID, i3):
    return str(
        (
            int(ID[0])
            + int(i3[0])
            + int(i3[2])
            + int(i3[4])
            + int(i3[6])
            + int(i3[8])
            + int(i3[10])
        )
        % 3
    )


def normalize_orientation(seq):
    candidates = []
    for oriented in (seq, revcomp(seq)):
        if oriented[0] in "AT" and oriented[-1] in "CG":
            candidates.append(oriented)
    if len(candidates) != 1:
        return None
    return candidates[0]


def decode_index(seq):
    if len(seq) != EXPECTED_LEN or set(seq) - set(BASES):
        return None, "length_or_symbol"

    oriented = normalize_orientation(seq)
    if oriented is None:
        return None, "outer_orientation"

    if any(oriented[j] == oriented[j - 1] for j in range(1, EXPECTED_LEN)):
        return None, "homopolymer"

    core = oriented[1:116]
    payload = core[:-15]
    ix = core[-15:]

    try:
        trits = DNA_INV[ix[0]][payload[-1]]
        for j in range(1, 15):
            trits += DNA_INV[ix[j]][ix[j - 1]]
    except KeyError:
        return None, "illegal_transition"

    ID = trits[:2]
    i3 = trits[2:14]
    P = trits[14]

    if ID != EXPECTED_ID:
        return None, "wrong_file_id"
    if P != checksum(ID, i3):
        return None, "checksum"

    index = int(i3, 3)
    if not (0 <= index < EXPECTED_COUNT):
        return None, "index_out_of_range"

    return (index, payload, oriented), None


def constrained_consensus(weighted_payloads):
    if not weighted_payloads:
        return None

    L = len(weighted_payloads[0][0])
    votes = [[0, 0, 0, 0] for _ in range(L)]

    for payload, weight in weighted_payloads:
        if len(payload) != L:
            continue
        for pos, base in enumerate(payload):
            votes[pos][BASE_TO_INT[base]] += weight

    neg_inf = -10**30
    prev_score = [votes[0][b] for b in range(4)]
    back = bytearray(L * 4)

    for pos in range(1, L):
        cur_score = [neg_inf] * 4
        for b in range(4):
            best_prev = None
            best_score = neg_inf
            for pb in range(4):
                if pb == b:
                    continue
                score = prev_score[pb]
                if score > best_score or (score == best_score and (best_prev is None or pb < best_prev)):
                    best_score = score
                    best_prev = pb
            cur_score[b] = best_score + votes[pos][b]
            back[pos * 4 + b] = best_prev
        prev_score = cur_score

    last = max(range(4), key=lambda b: (prev_score[b], -b))
    path = [last]
    for pos in range(L - 1, 0, -1):
        last = back[pos * 4 + last]
        path.append(last)
    path.reverse()
    return "".join(BASES[b] for b in path)


def canonical_fragment(stored_payload, index):
    i3 = base3_12(index)
    IX = EXPECTED_ID + i3 + checksum(EXPECTED_ID, i3)

    ix_chars = [DNA_TABLE[stored_payload[-1]][int(IX[0])]]
    for trit in IX[1:]:
        ix_chars.append(DNA_TABLE[ix_chars[-1]][int(trit)])
    ix = "".join(ix_chars)

    first = stored_payload[0]
    if first == "A":
        prefix = "T"
    elif first == "T":
        prefix = "A"
    else:
        prefix = "A"

    if ix[-1] == "C":
        suffix = "G"
    elif ix[-1] == "G":
        suffix = "C"
    else:
        suffix = "G"

    seq = prefix + stored_payload + ix + suffix
    if len(seq) != EXPECTED_LEN:
        raise AssertionError("canonical Goldman fragment is not 117 nt")
    return seq


# A. Retain full-length merged reads and count exact variants.
counts = Counter()
full_length_reads = 0

with open(fastq_file) as f, open(full_out, "w") as full_handle:
    for line_no, line in enumerate(f, start=1):
        if line_no % 4 != 2:
            continue
        seq = line.strip().upper()
        if len(seq) != EXPECTED_LEN or set(seq) - set(BASES):
            continue
        full_handle.write(seq + "\n")
        counts[seq] += 1
        full_length_reads += 1

with open(counts_out, "w") as out:
    for seq, count in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
        out.write("%d %s\n" % (count, seq))

# B. Validate embedded Goldman indices and retain top variants per index.
heaps = defaultdict(list)
index_support = [0] * EXPECTED_COUNT
index_unique = [0] * EXPECTED_COUNT
rejection = Counter()
valid_reads = 0
valid_unique = 0

for seq, count in counts.items():
    parsed, reason = decode_index(seq)
    if parsed is None:
        rejection[reason] += count
        continue

    index, payload, oriented = parsed
    valid_reads += count
    valid_unique += 1
    index_support[index] += count
    index_unique[index] += 1

    item = (count, oriented, payload)
    heap = heaps[index]
    if len(heap) < TOP_K:
        heapq.heappush(heap, item)
    elif item > heap[0]:
        heapq.heapreplace(heap, item)

# C. Weighted within-index consensus.
local_s5 = {}
for index, heap in heaps.items():
    weighted = [(payload, count) for count, _seq, payload in heap]
    payload = constrained_consensus(weighted)
    if payload is None:
        continue
    local_s5[index] = revcomp(payload) if index % 2 == 1 else payload

missing_before_overlap = EXPECTED_COUNT - len(local_s5)

# D. Global consensus using the inherent 75-nt Goldman overlaps.
global_len = 25 * EXPECTED_COUNT + 75
votes = bytearray(global_len * 4)
coverage = bytearray(global_len)

for index, payload in local_s5.items():
    start = 25 * index
    for offset, base in enumerate(payload):
        pos = start + offset
        votes[pos * 4 + BASE_TO_INT[base]] += 1
        coverage[pos] += 1

uncovered_positions = sum(1 for c in coverage if c == 0)
complete = uncovered_positions == 0 and EXPECTED_COUNT > 0

global_s5 = None
if complete:
    L = global_len
    neg_inf = -10**30
    prev_score = [votes[b] for b in range(4)]
    back = bytearray(L * 4)

    for pos in range(1, L):
        cur_score = [neg_inf] * 4
        base_offset = pos * 4
        for b in range(4):
            best_prev = None
            best_score = neg_inf
            for pb in range(4):
                if pb == b:
                    continue
                score = prev_score[pb]
                if score > best_score or (score == best_score and (best_prev is None or pb < best_prev)):
                    best_score = score
                    best_prev = pb
            cur_score[b] = best_score + votes[base_offset + b]
            back[base_offset + b] = best_prev
        prev_score = cur_score

    last = max(range(4), key=lambda b: (prev_score[b], -b))
    path = bytearray(L)
    path[L - 1] = last
    for pos in range(L - 1, 0, -1):
        last = back[pos * 4 + last]
        path[pos - 1] = last
    global_s5 = "".join(BASES[b] for b in path)

# E. Regenerate one canonical Goldman fragment for every expected index.
reconstructed_count = 0
with open(decoder_out, "w") as out:
    if global_s5 is not None:
        for index in range(EXPECTED_COUNT):
            start = 25 * index
            s5_payload = global_s5[start:start + 100]
            stored_payload = revcomp(s5_payload) if index % 2 == 1 else s5_payload
            fragment = canonical_fragment(stored_payload, index)

            parsed, reason = decode_index(fragment)
            if parsed is None or parsed[0] != index:
                raise AssertionError(
                    "self-validation failed for reconstructed Goldman index %d (%s)"
                    % (index, reason)
                )

            out.write(fragment + "\n")
            reconstructed_count += 1

with open(support_out, "w") as out:
    out.write("index\tvalid_read_support\tvalid_unique_variants\tlocal_consensus\n")
    for index in range(EXPECTED_COUNT):
        out.write(
            "%d\t%d\t%d\t%d\n"
            % (
                index,
                index_support[index],
                index_unique[index],
                1 if index in local_s5 else 0,
            )
        )

summary = {
    "expected_fragments": EXPECTED_COUNT,
    "full_length_reads": full_length_reads,
    "unique_full_length_sequences": len(counts),
    "index_valid_reads": valid_reads,
    "index_valid_unique_sequences": valid_unique,
    "indices_with_local_consensus": len(local_s5),
    "missing_indices_before_overlap": missing_before_overlap,
    "global_consensus_length": global_len,
    "uncovered_global_positions": uncovered_positions,
    "reconstructed_fragments": reconstructed_count,
    "status": "COMPLETE" if reconstructed_count == EXPECTED_COUNT else "INCOMPLETE",
}

for reason in sorted(rejection):
    summary["rejected_%s_reads" % reason] = rejection[reason]

with open(summary_out, "w") as out:
    out.write("metric\tvalue\n")
    for key, value in summary.items():
        out.write("%s\t%s\n" % (key, value))

print("--- GOLDMAN RECONSTRUCTION AUDIT ---")
for key, value in summary.items():
    print("%-36s %s" % (key + ":", value))
print("------------------------------------")
PY

        local full_length_reads
        local decoder_count
        local total_read_pairs
        full_length_reads=$(wc -l < "$full_seq")
        decoder_count=$(wc -l < "$decoder_input")
        total_read_pairs=$(awk 'END {print int(NR/4)}' "$seq_r1")

        local pcpm dropout_pct
        read -r pcpm dropout_pct <<< "$(
            "$SIM_PYTHON" - "$ENCODED_FILE" "$full_seq" "$total_read_pairs" <<'PY'
import sys
from collections import Counter

encoded_file, observed_file, total_reads = sys.argv[1], sys.argv[2], int(sys.argv[3])

with open(encoded_file) as f:
    encoded = [line.strip().upper() for line in f if line.strip()]

with open(observed_file) as f:
    observed = Counter(line.strip().upper() for line in f if line.strip())

if not encoded or total_reads <= 0:
    print("0.000000 100.000000")
    raise SystemExit

perfect_counts = [observed[oligo] for oligo in encoded]
perfect_total = sum(perfect_counts)
observed_oligos = sum(count > 0 for count in perfect_counts)

pcpm = (perfect_total * 1_000_000.0) / (total_reads * len(encoded))
dropout_pct = 100.0 * (len(encoded) - observed_oligos) / len(encoded)

print(f"{pcpm:.2f} {dropout_pct:.2f}")
PY
        )"

        echo "Full-length reads   : $full_length_reads"
        echo "Decoder input count : $decoder_count"
        echo "PCPM                : $pcpm"
        echo "Dropout             : $dropout_pct%"

        # ---------------- Decode ----------------
        echo "[8/8] Goldman decode"
        local decode_status="DECODE_FAILED"
        local decode_success="FAILURE"
        rm -f "$recovered_file"

        # invoke Goldman decode only when reconstruction produced exactly one fragment per encoded index.
        if (( decoder_count == ENCODED_COUNT )); then
            if goldman_decode "$decoder_input" "$recovered_file"; then
                if [[ -f "$recovered_file" ]] && cmp -s "$INPUT_FILE" "$recovered_file"; then
                    decode_status="SUCCESS"
                    decode_success="SUCCESS"
                elif [[ -f "$recovered_file" ]]; then
                    decode_status="DIFFERENT"
                fi
            fi
        else
            decode_status="INCOMPLETE_RECONSTRUCTION"
            echo "Goldman decoder skipped: reconstructed $decoder_count / $ENCODED_COUNT fragments." >&2
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
# Run all sequencing coverage levels x replicates sequentially
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


