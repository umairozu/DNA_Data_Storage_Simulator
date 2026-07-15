#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")"

mkdir -p cutadapt_files bwa_files

echo "Starting Bioinformatics Pipeline (Cutadapt -> BWA -> covered-window consensus -> DNA Fountain)..."

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

echo "Step 1: Trimming adapters/primers with cutadapt..."

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
    echo "Error: test_file.tar.gz.dna_order.txt not found in ./dna-fountain or ./test_files"
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

echo "Step 4: Building read-derived covered-window consensus from BWA alignments..."
python - <<'PY'
import os
import re
from collections import Counter, defaultdict

# Simulator-friendly defaults:
# - Do not trust/use FASTQ base quality; each aligned base gets one vote.
# - Do not reject reads by MAPQ unless user explicitly sets CONS_MIN_MAPQ.
# - Do not trim fixed primer lengths after consensus. Cutadapt already handled that.
MIN_MAPQ = int(os.environ.get("CONS_MIN_MAPQ", "0"))
MIN_DEPTH = int(os.environ.get("CONS_MIN_DEPTH", "3"))
MIN_FRAC = float(os.environ.get("CONS_MIN_FRAC", "0.60"))

# 0 = auto. If you know the decoder oligo length, set for stricter output, e.g.:
#   CONS_EXPECT_LEN=100 ./run_pipeline_bwa_consensus_autocrop.sh
EXPECT_LEN = int(os.environ.get("CONS_EXPECT_LEN", "0"))

# Safety: no fixed trimming by default. This script writes the consensus of the
# covered read window, because reads have already been cutadapt-trimmed before BWA.
TRIM_LEFT = int(os.environ.get("CONS_TRIM_LEFT", "0"))
TRIM_RIGHT = int(os.environ.get("CONS_TRIM_RIGHT", "0"))
if TRIM_LEFT or TRIM_RIGHT:
    print(f"WARNING: fixed trimming requested: CONS_TRIM_LEFT={TRIM_LEFT}, CONS_TRIM_RIGHT={TRIM_RIGHT}")

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
BASES = "ACGT"
SKIP_FLAGS = 0x4 | 0x100 | 0x800  # unmapped, secondary, supplementary


def read_fasta(path):
    refs = {}
    name = None
    parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    refs[name] = "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name is not None:
            refs[name] = "".join(parts).upper()
    return refs


def oligo_sort_key(name):
    m = re.search(r"(\d+)$", name)
    return int(m.group(1)) if m else name


def best_consensus_window(cons, depths, min_depth, expect_len=0):
    """Return (start, end) for the decoder region.

    We do NOT output the whole reference oligo when reads only cover the cutadapt-trimmed
    payload region. We output the largest confident covered interval, or the best fixed
    length interval if CONS_EXPECT_LEN is set.
    """
    good = [(b != "N" and depths[i] >= min_depth) for i, b in enumerate(cons)]
    n = len(cons)

    if expect_len and expect_len > 0:
        if expect_len > n:
            return None
        best = None
        best_score = None
        for s in range(0, n - expect_len + 1):
            e = s + expect_len
            # Prefer windows with zero bad positions, then highest total depth.
            bad = sum(1 for i in range(s, e) if not good[i])
            total_depth = sum(depths[s:e])
            score = (bad, -total_depth)
            if best_score is None or score < best_score:
                best_score = score
                best = (s, e)
        if best is None:
            return None
        s, e = best
        if any(not good[i] for i in range(s, e)):
            return None
        return best

    # No expected length: choose the longest contiguous high-confidence block.
    best_s = best_e = None
    i = 0
    while i < n:
        while i < n and not good[i]:
            i += 1
        s = i
        while i < n and good[i]:
            i += 1
        e = i
        if e > s:
            if best_s is None or (e - s, sum(depths[s:e])) > (best_e - best_s, sum(depths[best_s:best_e])):
                best_s, best_e = s, e
    if best_s is None:
        return None
    return best_s, best_e


refs = read_fasta("bwa_files/oligo_ref.fa")
if not refs:
    raise SystemExit("No reference oligos found in bwa_files/oligo_ref.fa")

votes = {name: [Counter() for _ in seq] for name, seq in refs.items()}
depths = {name: [0 for _ in seq] for name, seq in refs.items()}
reads_per_oligo = Counter()
aligned_bases_per_oligo = Counter()
skipped = Counter()
used_records = 0

with open("bwa_files/aln.sam") as f:
    for line in f:
        if not line or line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            skipped["malformed"] += 1
            continue

        try:
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3]) - 1
            mapq = int(fields[4])
            cigar = fields[5]
            seq = fields[9].upper()
        except Exception:
            skipped["malformed"] += 1
            continue

        if flag & SKIP_FLAGS:
            skipped["unmapped_secondary_or_supplementary"] += 1
            continue
        if rname == "*" or rname not in refs:
            skipped["unknown_reference"] += 1
            continue
        if mapq < MIN_MAPQ:
            skipped["low_mapq"] += 1
            continue
        if cigar == "*":
            skipped["no_cigar"] += 1
            continue

        rpos = pos
        qpos = 0
        ref_len = len(refs[rname])
        added = 0

        for n_s, op in CIGAR_RE.findall(cigar):
            n = int(n_s)
            if op in "M=X":
                for i in range(n):
                    rp = rpos + i
                    qp = qpos + i
                    if 0 <= rp < ref_len and qp < len(seq):
                        b = seq[qp]
                        if b in BASES:
                            votes[rname][rp][b] += 1
                            depths[rname][rp] += 1
                            added += 1
                rpos += n
                qpos += n
            elif op in "IS":
                qpos += n
            elif op in "DN":
                rpos += n
            elif op in "HP":
                pass

        if added > 0:
            reads_per_oligo[rname] += 1
            aligned_bases_per_oligo[rname] += added
            used_records += 1
        else:
            skipped["no_aligned_bases"] += 1

passed_records = []
summary_rows = []
window_lengths = Counter()

for name in sorted(refs, key=oligo_sort_key):
    refseq = refs[name]
    cons_chars = []
    ambiguous = 0
    low_depth = 0

    for i in range(len(refseq)):
        depth = depths[name][i]
        total_vote = sum(votes[name][i].values())
        if depth < MIN_DEPTH or total_vote == 0:
            cons_chars.append("N")
            low_depth += 1
            continue

        best_base, best_vote = votes[name][i].most_common(1)[0]
        frac = best_vote / total_vote
        if frac < MIN_FRAC:
            cons_chars.append("N")
            ambiguous += 1
        else:
            cons_chars.append(best_base)

    cons = "".join(cons_chars)
    window = best_consensus_window(cons, depths[name], MIN_DEPTH, EXPECT_LEN)
    written = "NO"
    decoder_seq = ""
    start_1based = 0
    end_1based = 0
    min_decoder_depth = 0
    mean_decoder_depth = 0.0

    if window is not None:
        s, e = window
        # Optional fixed trimming, off by default and not needed when cutadapt has already trimmed.
        s2 = s + TRIM_LEFT
        e2 = e - TRIM_RIGHT if TRIM_RIGHT else e
        if 0 <= s2 < e2 <= len(cons):
            decoder_seq = cons[s2:e2]
            decoder_depths = depths[name][s2:e2]
            if decoder_seq and "N" not in decoder_seq:
                passed_records.append((name, decoder_seq, cons, reads_per_oligo[name], s2, e2))
                window_lengths[len(decoder_seq)] += 1
                written = "YES"
                start_1based = s2 + 1
                end_1based = e2
                min_decoder_depth = min(decoder_depths) if decoder_depths else 0
                mean_decoder_depth = sum(decoder_depths) / len(decoder_depths) if decoder_depths else 0.0

    mismatches_to_ref = sum(1 for a, b in zip(cons, refseq) if a != "N" and a != b)
    covered_positions = sum(1 for d in depths[name] if d >= MIN_DEPTH)
    summary_rows.append((
        name,
        reads_per_oligo[name],
        len(refseq),
        covered_positions,
        start_1based,
        end_1based,
        len(decoder_seq),
        min_decoder_depth,
        f"{mean_decoder_depth:.2f}",
        low_depth,
        ambiguous,
        mismatches_to_ref,
        written,
    ))

# Sort decoder input sequences by length from longest to shortest before DNA Fountain.
# Tie-breakers: higher mapped-read support first, then oligo numeric order for reproducibility.
# passed_records tuple: (name, decoder_seq, full_cons, nreads, start, end)
passed_records.sort(key=lambda x: (-len(x[1]), -x[3], oligo_sort_key(x[0])))

with open("bwa_files/decoder_input.seq", "w") as out:
    for _, decoder_seq, _, _, _, _ in passed_records:
        out.write(decoder_seq + "\n")

with open("bwa_files/consensus_decoder_sequences.fa", "w") as out:
    for name, decoder_seq, _, nreads, s, e in passed_records:
        out.write(f">{name} mapped_reads={nreads} covered_window={s+1}-{e}\n{decoder_seq}\n")

with open("bwa_files/consensus_full_reference_coordinates.fa", "w") as out:
    for name, _, full_cons, nreads, _, _ in passed_records:
        out.write(f">{name} mapped_reads={nreads}\n{full_cons}\n")

with open("bwa_files/consensus_summary.tsv", "w") as out:
    cols = [
        "oligo_id", "mapped_records_used", "reference_len", "covered_positions_min_depth",
        "decoder_start_1based", "decoder_end_1based", "decoder_seq_len",
        "min_decoder_depth", "mean_decoder_depth", "low_depth_positions_full_ref",
        "ambiguous_positions_full_ref", "mismatches_to_reference_full_ref", "written_to_decoder_input"
    ]
    out.write("\t".join(cols) + "\n")
    for row in summary_rows:
        out.write("\t".join(map(str, row)) + "\n")

print(f"Consensus settings: MIN_MAPQ={MIN_MAPQ}, MIN_DEPTH={MIN_DEPTH}, MIN_FRAC={MIN_FRAC}, EXPECT_LEN={EXPECT_LEN}, TRIM_LEFT={TRIM_LEFT}, TRIM_RIGHT={TRIM_RIGHT}")
print("Quality handling: unweighted consensus; FASTQ base qualities ignored. MAPQ cutoff defaults to 0 for simulator data.")
print("Window handling: no fixed primer/payload trimming by default; decoder_input.seq is built from the actual covered consensus window.")
print(f"Used SAM records: {used_records}")
print(f"Skipped records: {dict(skipped)}")
print(f"Wrote {len(passed_records)} consensus-derived sequences to bwa_files/decoder_input.seq")
print(f"Decoder sequence lengths written: {dict(window_lengths)}")
print("Decoder input ordering: sequences sorted by length, longest to shortest; ties by mapped-read count.")
print("Also wrote bwa_files/consensus_decoder_sequences.fa, bwa_files/consensus_full_reference_coordinates.fa, and bwa_files/consensus_summary.tsv")
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
