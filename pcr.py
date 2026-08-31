import argparse
import json
import os
import random
from collections import defaultdict

import Levenshtein
import numpy as np

from GC_content import gc_error_probability
from Error_module import Error_simulation



err_rates = {
    "1": {
        "raw_rate": 0.000043,
        "substitution": 0.99,
        "deletion": 0.01,
        "insertion": 0.0,
    },  # Taq
    "2": {
        "raw_rate": 0.0000024,
        "substitution": 1.0,
        "deletion": 0.0,
        "insertion": 0.0,
    },  # Pwo
    "3": {
        "raw_rate": 0.0000028,
        "substitution": 1.0,
        "deletion": 0.0,
        "insertion": 0.0,
    },  # Pfu
    "4": {
        "raw_rate": 0.0000026,
        "substitution": 0.84,
        "deletion": 0.08,
        "insertion": 0.08,
    },  # Phusion
    "5": {
        "raw_rate": 0.0,
        "substitution": 0.3334,
        "deletion": 0.3333,
        "insertion": 0.3333,
    },  # None / mutation-free polymerase control
}

POLYMERASE_NAMES = {
    "1": "Taq",
    "2": "Pwo",
    "3": "Pfu",
    "4": "Phusion",
    "5": "None",
}


mutation_attributes = {
    # Taq
    "1": {
        "deletion": {
            "position": {"homopolymer": 0, "random": 1},
            "pattern": {"G": 0.25, "C": 0.25, "A": 0.25, "T": 0.25},
        },
        "insertion": {
            "position": {"homopolymer": 0, "random": 1},
            "pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25},
        },
        "substitution": {
            "pattern": {
                "A": {"G": 0.97, "T": 0.01, "C": 0.02},
                "T": {"C": 0.97, "A": 0.01, "G": 0.02},
                "G": {"A": 1.0, "T": 0.0, "C": 0.0},
                "C": {"T": 1.0, "G": 0.0, "A": 0.0},
            }
        },
    },
    # Pwo
    "2": {
        "deletion": {
            "position": {"homopolymer": 0, "random": 1},
            "pattern": {"G": 0.25, "C": 0.25, "A": 0.25, "T": 0.25},
        },
        "insertion": {
            "position": {"homopolymer": 0, "random": 1},
            "pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25},
        },
        "substitution": {
            "pattern": {
                "A": {"G": 1.0, "T": 0.0, "C": 0.0},
                "T": {"C": 0.67, "A": 0.33, "G": 0.0},
                "G": {"A": 0.57, "T": 0.0, "C": 0.43},
                "C": {"T": 1.0, "G": 0.0, "A": 0.0},
            }
        },
    },
    # Pfu
    "3": {
        "deletion": {
            "position": {"homopolymer": 0, "random": 1},
            "pattern": {"G": 0.25, "C": 0.25, "A": 0.25, "T": 0.25},
        },
        "insertion": {
            "position": {"homopolymer": 0, "random": 1},
            "pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25},
        },
        "substitution": {
            "pattern": {
                "A": {"G": 0.75, "T": 0.25, "C": 0.0},
                "T": {"C": 0.75, "A": 0.25, "G": 0.0},
                "G": {"A": 1.0, "T": 0.0, "C": 0.0},
                "C": {"T": 1.0, "G": 0.0, "A": 0.0},
            }
        },
    },
    # Phusion
    "4": {
        "deletion": {
            "position": {"homopolymer": 0, "random": 1},
            "pattern": {"G": 0.25, "C": 0.25, "A": 0.25, "T": 0.25},
        },
        "insertion": {
            "position": {"homopolymer": 0, "random": 1},
            "pattern": {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25},
        },
        "substitution": {
            "pattern": {
                "A": {"G": 1.0, "T": 0.0, "C": 0.0},
                "T": {"C": 1.0, "A": 0.0, "G": 0.0},
                "G": {"A": 1.0, "T": 0.0, "C": 0.0},
                "C": {"T": 1.0, "G": 0.0, "A": 0.0},
            }
        },
    },
    # None
    "5": {
        "deletion": {
            "pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
            "position": {"homopolymer": 0.5, "random": 0.5},
        },
        "insertion": {
            "pattern": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
            "position": {"homopolymer": 0.5, "random": 0.5},
        },
        "substitution": {"pattern": {}},
    },
}


MUTATION_TYPES = ("insertion", "deletion", "substitution")
MUTATION_SCALE_PRESETS = {
    "0": 0.0,  # mutation-free polymerase control
    "1": 0.5,  # half reference error rate
    "2": 1.0,  # reference error rate
    "3": 2.0,  # double reference error rate / stress condition
}

POOL_DTYPE = np.dtype(
    [
        ("pid", np.int64),
        ("count", np.int64),
        ("seq", object),
    ]
)



def percentage_0_to_100(value):
    value = float(value)
    if not 0.0 <= value <= 100.0:
        raise argparse.ArgumentTypeError("value must be between 0 and 100")
    return value


def positive_percentage(value):
    value = float(value)
    if not 0.0 < value <= 100.0:
        raise argparse.ArgumentTypeError("value must satisfy 0 < value <= 100")
    return value


def nonnegative_scale(value):
    value = float(value)
    if not 0.0 <= value <= 10.0:
        raise argparse.ArgumentTypeError("custom mutation scale must be between 0 and 10")
    return value


def positive_int(value):
    value = int(value)
    if value <= 0:
        raise argparse.ArgumentTypeError("value must be > 0")
    return value


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "PCR simulator with stochastic amplification and polymerase-error "
            "propagation. Mutation abundance emerges from polymerase fidelity "
            "instead of a fixed 80/10/10 copy split."
        )
    )

    parser.add_argument(
        "--s",
        type=positive_percentage,
        default=100.0,
        help="PCR sampling percentage from the input pool (default: 100)",
    )
    parser.add_argument(
        "--n",
        type=positive_int,
        default=30,
        help="number of PCR cycles (default: 30)",
    )
    parser.add_argument(
        "--polymerase",
        choices=sorted(POLYMERASE_NAMES.keys()),
        default="1",
        help="1=Taq, 2=Pwo, 3=Pfu, 4=Phusion, 5=None (default: 1)",
    )
    parser.add_argument(
        "--mut",
        choices=sorted(MUTATION_SCALE_PRESETS.keys()),
        default="2",
        help=(
            "polymerase error-rate preset: 0=0x, 1=0.5x, "
            "2=1x reference, 3=2x reference (default: 2)"
        ),
    )
    parser.add_argument(
        "--c",
        type=nonnegative_scale,
        default=None,
        help=(
            "optional custom multiplier for the polymerase raw error rate "
            "(0-10); overrides --mut"
        ),
    )
    parser.add_argument(
        "--chimera_rate",
        type=percentage_0_to_100,
        default=5.0,
        help="percentage of final PCR molecules reassigned to chimeras (default: 0)",
    )
    parser.add_argument(
        "--variant_cap",
        type=positive_int,
        default=3,
        help=(
            "maximum representative mutant variants stored per parent oligo "
            "(default: 3; counts remain conserved)"
        ),
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="optional RNG seed for reproducible simulations",
    )
    parser.add_argument(
        "--keep_intermediate",
        action="store_true",
        help="keep PCR intermediate/debug files",
    )
    parser.add_argument(
        "--in_file",
        required=True,
        help="input pool CSV: parent_id,count,sequence",
    )
    parser.add_argument(
        "--out_file",
        default="pcr_output.txt",
        help="output filename written under ./pcr unless an absolute path is supplied",
    )
    return parser.parse_args()



def safe_sum(values):
    return sum(int(value) for value in values)


def resolve_output_path(pcr_dir, requested_path):
    if os.path.isabs(requested_path):
        return requested_path
    return os.path.join(pcr_dir, requested_path)


def read_pool(path):
    """Read parent_id,count,sequence"""
    records = []

    with open(path, "r") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip("\n\r")
            if not line or "," not in line:
                continue

            parts = line.split(",", maxsplit=2)
            if len(parts) != 3:
                continue

            try:
                pid = int(parts[0].strip())
                count = int(parts[1].strip())
            except ValueError:
                continue

            seq = parts[2].strip()
            if count > 0 and seq:
                records.append((pid, count, seq))

    if not records:
        raise ValueError(f"No valid PCR input records found in: {path}")

    return np.array(records, dtype=POOL_DTYPE)


def write_pool(path, data):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as handle:
        handle.write("parent_id,count,sequence\n")
        for row in data:
            handle.write(f"{int(row['pid'])},{int(row['count'])},{row['seq']}\n")


def aggregate_records(records):
    """Merge duplicate (parent_id, sequence) records while conserving counts."""
    merged = defaultdict(int)

    for pid, count, seq in records:
        count = int(count)
        if count <= 0:
            continue
        merged[(int(pid), str(seq))] += count

    if not merged:
        return np.empty(0, dtype=POOL_DTYPE)

    result = np.empty(len(merged), dtype=POOL_DTYPE)
    for i, ((pid, seq), count) in enumerate(merged.items()):
        result[i] = (pid, count, seq)
    return result


# -----------------------------------------------------------------------------
# PCR efficiency model
# -----------------------------------------------------------------------------
def sequence_efficiency(initial_efficiency, sequence, primer_f, primer_r, pf_len, pr_len):
    """
    Compute initial PCR amplification efficiency from global GC behavior and
    primer similarity.
    """
    eff = float(initial_efficiency)

    gc_error = gc_error_probability(sequence)
    if gc_error <= 0.0:
        pass
    elif gc_error <= 0.10:
        eff *= 0.90
    elif gc_error <= 0.25:
        eff *= 0.75
    elif gc_error <= 0.40:
        eff *= 0.60
    elif gc_error <= 0.60:
        eff *= 0.40
    else:
        eff *= 0.20

    seq_pf = sequence[:pf_len]
    seq_pr = sequence[-pr_len:]

    sim_pf = Levenshtein.ratio(primer_f, seq_pf)
    sim_pr = Levenshtein.ratio(primer_r, seq_pr)
    primer_similarity = min(sim_pf, sim_pr)

    if primer_similarity < 0.80:
        eff *= 0.50
    elif primer_similarity < 0.90:
        eff *= 0.60
    elif primer_similarity < 0.95:
        eff *= 0.70

    return float(np.clip(eff, 0.0, 1.0))


def cycle_efficiency(current_eff, cycle_index, c1, c2):
    """
    Stage 1: constant exponential efficiency.
    Stage 2: efficiency decreases by 10% before each cycle.
    Stage 3: amplification factor (0.4 + eff), therefore the equivalent new-copy probability is max(0, eff - 0.6).
    """
    if cycle_index < c1:
        return current_eff, current_eff

    if cycle_index < c2:
        current_eff = current_eff * 0.90
        return current_eff, current_eff

    plateau_eff = np.maximum(current_eff - 0.60, 0.0)
    return current_eff, plateau_eff


# -----------------------------------------------------------------------------
# Polymerase error model
# -----------------------------------------------------------------------------
def scaled_error_model(polymerase_key, mutation_scale):
    model = dict(err_rates[polymerase_key])
    model["raw_rate"] = min(max(model["raw_rate"] * mutation_scale, 0.0), 1.0)
    return model


def per_new_copy_error_probability(sequence_length, error_model):
    """
    Probability that a newly synthesized copy receives >=1 polymerase error.
    """
    if sequence_length <= 0 or error_model["raw_rate"] <= 0.0:
        return 0.0

    p_no_error = 1.0
    for mutation_type in MUTATION_TYPES:
        per_base_type_rate = (
            error_model["raw_rate"] * error_model[mutation_type]
        )
        per_base_type_rate = min(max(per_base_type_rate, 0.0), 1.0)
        p_no_error *= (1.0 - per_base_type_rate) ** sequence_length

    return float(np.clip(1.0 - p_no_error, 0.0, 1.0))


def choose_mutation_type(error_model):
    mutation_types = []
    weights = []

    for mutation_type in MUTATION_TYPES:
        weight = float(error_model[mutation_type])
        if weight > 0.0:
            mutation_types.append(mutation_type)
            weights.append(weight)

    if not mutation_types:
        return None

    total = sum(weights)
    normalized = [weight / total for weight in weights]
    return random.choices(mutation_types, weights=normalized, k=1)[0]


def force_one_error_module_mutation(sequence, attributes, error_model, max_attempts=25):
    """
    Generate one representative mutant using Error_module's own mutation
    handlers and mutation attributes.

    The probability that a molecule belongs to the mutant population has already
    been decided during PCR amplification. Therefore this function deliberately
    conditions on 'an error occurred' and samples the error TYPE from the
    configured insertion/deletion/substitution proportions.
    """
    for _ in range(max_attempts):
        mutation_type = choose_mutation_type(error_model)
        if mutation_type is None:
            return None

        simulator = Error_simulation(
            sequence,
            "pcr",
            attribute=attributes,
            error_rate=error_model,
        )
        simulator.reset_visited()

        before = simulator.seq
        handler = {
            "insertion": simulator.insertion,
            "deletion": simulator.deletion,
            "substitution": simulator.substitution,
        }[mutation_type]
        handler()

        if simulator.seq != before:
            return simulator.seq

    return None


def build_representative_pool(
    pids,
    sequences,
    clean_counts,
    total_counts,
    attributes,
    error_model,
    variant_cap,
):
    """
    Compress the PCR population into original sequences plus a small number of
    representative mutant variants per parent. Counts remain exactly conserved.
    """
    records = []
    total_failed_variant_copies = 0

    for pid, seq, clean_count, total_count in zip(
        pids,
        sequences,
        clean_counts,
        total_counts,
    ):
        clean_count = int(clean_count)
        total_count = int(total_count)
        mutant_count = total_count - clean_count

        if clean_count > 0:
            records.append((int(pid), clean_count, seq))

        if mutant_count <= 0:
            continue

        number_of_variants = min(
            random.randint(1, variant_cap),
            mutant_count,
        )

        variants = []
        variant_attempts = 0
        max_variant_attempts = max(10, number_of_variants * 10)

        while len(variants) < number_of_variants and variant_attempts < max_variant_attempts:
            mutant_seq = force_one_error_module_mutation(
                seq,
                attributes,
                error_model,
            )
            variant_attempts += 1

            if mutant_seq is None or mutant_seq == seq:
                continue
            if mutant_seq not in variants:
                variants.append(mutant_seq)

        if not variants:
            # Conservative fallback: if no representative sequence can be generated, return those copies to the original sequence
            records.append((int(pid), mutant_count, seq))
            total_failed_variant_copies += mutant_count
            continue

        portions = np.random.multinomial(
            mutant_count,
            [1.0 / len(variants)] * len(variants),
        )

        for variant_seq, portion in zip(variants, portions):
            if portion > 0:
                records.append((int(pid), int(portion), variant_seq))

    return aggregate_records(records), total_failed_variant_copies


# -----------------------------------------------------------------------------
# Chimera model
# -----------------------------------------------------------------------------
def build_chimeras(pre_chimera_data, source_pids, source_sequences, source_counts, chimera_fraction):
    """
    Reassign a requested fraction of final PCR molecules to chimeric sequences.
    Molecule counts are conserved exactly.
    """
    if chimera_fraction <= 0.0 or len(pre_chimera_data) == 0:
        return pre_chimera_data.copy(), np.empty(0, dtype=POOL_DTYPE)

    reduced = pre_chimera_data.copy()
    removed_counts = np.random.binomial(reduced["count"], chimera_fraction)
    reduced["count"] -= removed_counts
    chimera_total = safe_sum(removed_counts)

    if chimera_total <= 0:
        return reduced[reduced["count"] > 0], np.empty(0, dtype=POOL_DTYPE)

    groups = defaultdict(list)
    for index, (seq, count) in enumerate(zip(source_sequences, source_counts)):
        if int(count) > 0 and len(seq) >= 2:
            groups[len(seq)].append(index)

    eligible_groups = {
        length: indices
        for length, indices in groups.items()
        if len(indices) >= 2
    }

    if not eligible_groups:
        reduced["count"] += removed_counts
        print("WARNING: chimera_rate > 0 but fewer than two compatible parent templates exist; no chimeras created.")
        return reduced, np.empty(0, dtype=POOL_DTYPE)

    lengths = list(eligible_groups.keys())
    group_weights = [
        safe_sum(source_counts[eligible_groups[length]])
        for length in lengths
    ]
    group_weights = np.asarray(group_weights, dtype=np.float64)
    group_weights /= group_weights.sum()

    number_of_chimera_variants = min(
        random.randint(10, 30),
        chimera_total,
    )

    portions = np.random.multinomial(
        chimera_total,
        [1.0 / number_of_chimera_variants] * number_of_chimera_variants,
    )

    chimera_records = []

    for portion in portions:
        if portion <= 0:
            continue

        chosen_length = int(np.random.choice(lengths, p=group_weights))
        candidate_indices = np.asarray(eligible_groups[chosen_length], dtype=np.int64)

        candidate_weights = np.asarray(
            [int(source_counts[index]) for index in candidate_indices],
            dtype=np.float64,
        )
        candidate_weights /= candidate_weights.sum()

        parent_indices = np.random.choice(
            candidate_indices,
            size=2,
            replace=False,
            p=candidate_weights,
        )

        i1, i2 = int(parent_indices[0]), int(parent_indices[1])
        seq1 = source_sequences[i1]
        seq2 = source_sequences[i2]

        breakpoint = random.randrange(1, chosen_length)
        chimera_seq = seq1[:breakpoint] + seq2[breakpoint:]

        chimera_records.append(
            (
                int(source_pids[i1]),
                int(portion),
                chimera_seq,
            )
        )

    chimera_data = aggregate_records(chimera_records)
    reduced = reduced[reduced["count"] > 0]

    assert safe_sum(reduced["count"]) + safe_sum(chimera_data["count"]) == safe_sum(pre_chimera_data["count"])

    return reduced, chimera_data



def main():
    args = parse_args()

    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    base_dir = os.getcwd()
    pcr_dir = os.path.join(base_dir, "pcr")
    os.makedirs(pcr_dir, exist_ok=True)

    metadata_path = os.path.join(base_dir, "metadata.json")
    with open(metadata_path, "r") as handle:
        metadata = json.load(handle)

    primer_f = metadata["primer_F"]
    primer_r = metadata["primer_R"]
    orig_length = int(metadata["orig_length"])
    pf_length = int(metadata["pF_length"])
    pr_length = int(metadata["pR_length"])

    sampling_fraction = float(args.s) / 100.0
    num_cycles = int(args.n)
    polymerase_key = args.polymerase
    polymerase_name = POLYMERASE_NAMES[polymerase_key]

    if args.c is not None:
        mutation_scale = float(args.c)
        mutation_scale_source = "custom --c"
    else:
        mutation_scale = MUTATION_SCALE_PRESETS[args.mut]
        mutation_scale_source = f"preset --mut {args.mut}"

    error_model = scaled_error_model(polymerase_key, mutation_scale)
    attributes = mutation_attributes[polymerase_key]
    chimera_fraction = float(args.chimera_rate) / 100.0

    input_path = args.in_file
    output_path = resolve_output_path(pcr_dir, args.out_file)

    intermediate_paths = {
        "dropouts": os.path.join(pcr_dir, "pcr_dropouts.txt"),
        "filtered": os.path.join(pcr_dir, "pcr_filtered.txt"),
        "sampled": os.path.join(pcr_dir, "pcr_sampled.txt"),
        "complete": os.path.join(pcr_dir, "pcr_complete.txt"),
    }

    print("\n--- PCR CONFIGURATION ---")
    print(f"Input:                         {input_path}")
    print(f"Output:                        {output_path}")
    print(f"Sampling:                      {args.s:.4g}%")
    print(f"PCR cycles:                    {num_cycles}")
    print(f"Polymerase:                    {polymerase_name} ({polymerase_key})")
    print(f"Reference raw error rate:      {err_rates[polymerase_key]['raw_rate']:.8g} / base / new copy")
    print(f"Mutation scale:                {mutation_scale:.4g}x ({mutation_scale_source})")
    print(f"Effective raw error rate:      {error_model['raw_rate']:.8g} / base / new copy")
    print(f"Chimera rate:                  {args.chimera_rate:.4g}%")
    print(f"Representative variant cap:    {args.variant_cap}")
    print(f"Seed:                          {args.seed}")
    print("-------------------------\n")

    # ------------------------------------------------------------------
    # Load and pre-filter input pool
    # ------------------------------------------------------------------
    np_data = read_pool(input_path)
    sequences = np_data["seq"]
    number_of_records = len(np_data)

    similarity_pf = np.fromiter(
        (Levenshtein.ratio(primer_f, seq[:pf_length]) for seq in sequences),
        dtype=np.float64,
        count=number_of_records,
    )
    similarity_pr = np.fromiter(
        (Levenshtein.ratio(primer_r, seq[-pr_length:]) for seq in sequences),
        dtype=np.float64,
        count=number_of_records,
    )
    length_difference = np.fromiter(
        (abs(len(seq) - orig_length) for seq in sequences),
        dtype=np.int64,
        count=number_of_records,
    )

    primer_dropout_mask = (similarity_pf < 0.75) & (similarity_pr < 0.75)
    length_dropout_mask = length_difference > 20
    dropout_mask = primer_dropout_mask | length_dropout_mask

    dropouts = np_data[dropout_mask]
    filtered = np_data[~dropout_mask]

    print("--- PCR PRE-FILTER DEBUG ---")
    print(f"Input oligos:                  {len(np_data)}")
    print(f"Primer dropouts only:          {np.count_nonzero(primer_dropout_mask & ~length_dropout_mask)}")
    print(f"Length dropouts only:          {np.count_nonzero(length_dropout_mask & ~primer_dropout_mask)}")
    print(f"Primer and length dropouts:    {np.count_nonzero(primer_dropout_mask & length_dropout_mask)}")
    print(f"Total dropouts:                {np.count_nonzero(dropout_mask)}")
    print(f"Surviving filtered oligos:     {len(filtered)}")
    print("----------------------------\n")

    write_pool(intermediate_paths["dropouts"], dropouts)

    if len(filtered) == 0:
        raise ValueError("PCR FAILED: every input oligo was removed by the pre-filter")

    write_pool(intermediate_paths["filtered"], filtered)

    # ------------------------------------------------------------------
    # PCR sampling
    # ------------------------------------------------------------------
    sampled = filtered.copy()
    sampled["count"] = np.random.binomial(
        sampled["count"],
        sampling_fraction,
    ).astype(np.int64)

    sampled = sampled[sampled["count"] > 0]
    if len(sampled) == 0:
        raise ValueError("PCR FAILED: sampling produced zero molecules")

    write_pool(intermediate_paths["sampled"], sampled)

    starting_molecules = safe_sum(sampled["count"])

    # ------------------------------------------------------------------
    # Per-oligo efficiency and per-new-copy polymerase error probability
    # ------------------------------------------------------------------
    initial_efficiency = 0.90
    base_efficiencies = np.fromiter(
        (
            sequence_efficiency(
                initial_efficiency,
                seq,
                primer_f,
                primer_r,
                pf_length,
                pr_length,
            )
            for seq in sampled["seq"]
        ),
        dtype=np.float64,
        count=len(sampled),
    )

    p_error_per_new_copy = np.fromiter(
        (
            per_new_copy_error_probability(len(seq), error_model)
            for seq in sampled["seq"]
        ),
        dtype=np.float64,
        count=len(sampled),
    )

    # ------------------------------------------------------------------
    # Stochastic PCR amplification with mutation inheritance
    # ------------------------------------------------------------------
    total_counts = sampled["count"].astype(np.int64, copy=True)
    clean_counts = sampled["count"].astype(np.int64, copy=True)
    current_eff = base_efficiencies.copy()

    c1 = num_cycles // 3
    c2 = c1 + num_cycles // 3
    max_int64 = np.iinfo(np.int64).max

    for cycle in range(num_cycles):
        current_eff, this_cycle_eff = cycle_efficiency(
            current_eff,
            cycle,
            c1,
            c2,
        )
        this_cycle_eff = np.clip(this_cycle_eff, 0.0, 1.0)

        mutant_counts_before = total_counts - clean_counts

        # Clean parents and mutant parents amplify independently.
        new_from_clean_parents = np.random.binomial(
            clean_counts,
            this_cycle_eff,
        ).astype(np.int64)

        new_from_mutant_parents = np.random.binomial(
            mutant_counts_before,
            this_cycle_eff,
        ).astype(np.int64)

        # Among copies synthesized from clean templates, some acquire their
        # first polymerase error. Mutant-parent descendants remain mutant.
        new_clean_copies = np.random.binomial(
            new_from_clean_parents,
            1.0 - p_error_per_new_copy,
        ).astype(np.int64)

        new_total = new_from_clean_parents + new_from_mutant_parents

        if np.any(total_counts > (max_int64 - new_total)):
            raise OverflowError(
                "PCR molecule counts exceeded np.int64 capacity. Reduce input "
                "copy counts/cycles or add a count-scaling strategy."
            )

        total_counts += new_total
        clean_counts += new_clean_copies

        assert np.all(clean_counts <= total_counts), (
            "Internal PCR accounting error: clean molecule count exceeded total count"
        )

    mutant_counts = total_counts - clean_counts

    # ------------------------------------------------------------------
    # Amplification summary file
    # ------------------------------------------------------------------
    with open(intermediate_paths["complete"], "w") as handle:
        handle.write(
            "parent_id,sampled_count,pcr_count,clean_count,mutant_count,"
            "final_internal_efficiency,p_error_per_new_copy,sequence\n"
        )
        for pid, sampled_count, total_count, clean_count, mutant_count, eff, p_err, seq in zip(
            sampled["pid"],
            sampled["count"],
            total_counts,
            clean_counts,
            mutant_counts,
            current_eff,
            p_error_per_new_copy,
            sampled["seq"],
        ):
            handle.write(
                f"{int(pid)},{int(sampled_count)},{int(total_count)},"
                f"{int(clean_count)},{int(mutant_count)},"
                f"{float(eff):.8f},{float(p_err):.10f},{seq}\n"
            )


    pre_chimera_data, failed_variant_copies = build_representative_pool(
        sampled["pid"],
        sampled["seq"],
        clean_counts,
        total_counts,
        attributes,
        error_model,
        args.variant_cap,
    )

    pcr_total = safe_sum(total_counts)
    clean_total = safe_sum(clean_counts)
    mutant_total = safe_sum(mutant_counts)
    represented_pre_chimera_total = safe_sum(pre_chimera_data["count"])

    assert clean_total + mutant_total == pcr_total
    assert represented_pre_chimera_total == pcr_total


    non_chimera_data, chimera_data = build_chimeras(
        pre_chimera_data,
        sampled["pid"],
        sampled["seq"],
        total_counts,
        chimera_fraction,
    )

    if len(chimera_data) > 0:
        final_data = aggregate_records(
            list(non_chimera_data) + list(chimera_data)
        )
    else:
        final_data = non_chimera_data

    final_total = safe_sum(final_data["count"])
    chimera_total = safe_sum(chimera_data["count"]) if len(chimera_data) else 0

    assert final_total == pcr_total, (
        "Final output does not conserve the post-PCR molecule count"
    )

    write_pool(output_path, final_data)


    weighted_mean_p_error = np.average(
        p_error_per_new_copy,
        weights=sampled["count"].astype(np.float64),
    )

    print("\n--- PCR / MUTATION AUDIT ---")
    print(f"Molecules entering PCR:        {starting_molecules}")
    print(f"Molecules after PCR:           {pcr_total}")
    print(f"Clean molecules:               {clean_total}")
    print(f"Mutant molecules:              {mutant_total}")
    print(f"Mutant fraction:               {(mutant_total / pcr_total) if pcr_total else 0.0:.6%}")
    print(f"Weighted P(error/new copy):    {weighted_mean_p_error:.6%}")
    print(f"P(error/new copy) min/max:     {p_error_per_new_copy.min():.6%} / {p_error_per_new_copy.max():.6%}")
    print(f"Representative PCR records:    {len(pre_chimera_data)}")
    print(f"Variant fallback copies:       {failed_variant_copies}")
    print(f"Chimera molecules:             {chimera_total}")
    print(f"Final output molecules:        {final_total}")
    print(f"Final output records:          {len(final_data)}")
    print("----------------------------\n")

    if not args.keep_intermediate:
        for path in intermediate_paths.values():
            if os.path.exists(path):
                os.remove(path)

    print("pcr.py run completed")


if __name__ == "__main__":
    main()