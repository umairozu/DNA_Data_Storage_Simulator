
import argparse
import math
import os
from collections import defaultdict

import numpy as np

from Arrhenius_decay import Arrhenius_decay


SEC_PER_WEEK = 7 * 24 * 3600

STORAGE_SCALE_PRESETS = {
    "0": 0.0,
    "1": 0.5,
    "2": 1.0,
    "3": 2.0,
}

POOL_DTYPE = np.dtype([
    ("pid", np.int64),
    ("count", np.int64),
    ("seq", object),
])

def pH_value(value):
    value = float(value)
    if not 0.0 <= value <= 14.0:
        raise argparse.ArgumentTypeError("pH must be between 0 and 14")
    return value


def rh_value(value):
    value = float(value)
    if not 0.0 <= value <= 100.0:
        raise argparse.ArgumentTypeError("relative humidity must be between 0 and 100")
    return value


def nonnegative_float(value):
    value = float(value)
    if value < 0.0:
        raise argparse.ArgumentTypeError("value must be >= 0")
    return value


def nonnegative_scale(value):
    value = float(value)
    if not 0.0 <= value <= 10.0:
        raise argparse.ArgumentTypeError("custom storage scale must be between 0 and 10")
    return value


def positive_int(value):
    value = int(value)
    if value <= 0:
        raise argparse.ArgumentTypeError("value must be > 0")
    return value


def parse_args():
    parser = argparse.ArgumentParser(
        description="DNA-tape storage simulator with validated Arrhenius population decay"
    )
    parser.add_argument("--temp", type=float, required=True, help="storage temperature in Celsius")
    parser.add_argument(
        "--ph", type=pH_value, default=7.0,
        help="pH; affects only optional --damage_model depurination (default: 7.0)",
    )
    parser.add_argument(
        "--rh", type=rh_value, default=50.0,
        help=("relative humidity percentage. The cassette-tape Arrhenius model was "
              "calibrated at 50%% RH; other values are logged but not corrected (default: 50)")
    )
    parser.add_argument(
        "--week", type=nonnegative_float, required=True,
        help="storage duration in weeks; zero is allowed for a baseline control",
    )
    parser.add_argument(
        "--encap", choices=["0", "1"], default="0",
        help="0=D-DNA  (decapsulated), 1=E-DNA  (ZIF encapsulated)",
    )
    parser.add_argument(
        "--mut", choices=sorted(STORAGE_SCALE_PRESETS), default="2",
        help="storage-rate scale: 0=0x, 1=0.5x, 2=1x reference, 3=2x (default: 2)",
    )
    parser.add_argument(
        "--c", type=nonnegative_scale, default=None,
        help="custom storage-rate multiplier 0-10; overrides --mut",
    )
    parser.add_argument(
        "--damage_model", choices=["none", "depurination"], default="depurination",
        help="sequence-damage model. 'depurination' enables an optional mechanistic stress model (default: none)"
    )
    parser.add_argument(
        "--damage_fate", choices=["variants", "dropout"], default="variants",
        help="optional depurination fate: representative deletion variants or dropout",
    )
    parser.add_argument(
        "--encap_damage_factor", default="bulk_ratio",
        help=("for --damage_model depurination with E-DNA: 'bulk_ratio', '1', '0', or a "
              "numeric 0-1 factor multiplying aqueous depurination.")
    )
    parser.add_argument(
        "--variant_cap", type=positive_int, default=8,
        help="max Monte-Carlo representative damaged variants per parent (default: 8)",
    )
    parser.add_argument(
        "--survival_mode", choices=["binomial", "expected"], default="binomial",
        help="realize bulk survival stochastically or by rounded expectation (default: binomial)",
    )
    parser.add_argument("--seed", type=int, default=None, help="optional RNG seed")
    parser.add_argument(
        "--arrhenius_xlsx", default="RawData.xlsx",
        help="cassette-tape RawData.xlsx path (default: RawData.xlsx)",
    )
    parser.add_argument("--in_file", required=True, help="input CSV: parent_id,count,sequence")
    parser.add_argument(
        "--out_file", default="storage_output.txt",
        help="output filename under ./storage unless absolute",
    )
    return parser.parse_args()


def resolve_existing_file(path):
    if os.path.isabs(path) and os.path.isfile(path):
        return path
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [os.path.abspath(path), os.path.join(script_dir, path)]
    for candidate in candidates:
        if os.path.isfile(candidate):
            return candidate
    raise FileNotFoundError("Could not find file. Checked: " + ", ".join(candidates))


def resolve_output_path(storage_dir, path):
    return path if os.path.isabs(path) else os.path.join(storage_dir, path)


def read_pool(path):
    records = []
    with open(path, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or "," not in line:
                continue
            parts = line.split(",", 2)
            if len(parts) != 3:
                continue
            try:
                pid = int(parts[0].strip())
                count = int(float(parts[1].strip()))
            except ValueError:
                continue
            sequence = parts[2].strip().upper()
            if count > 0 and sequence:
                records.append((pid, count, sequence))
    if not records:
        raise ValueError(f"No valid storage input records in {path}")
    return np.array(records, dtype=POOL_DTYPE)


def write_pool(path, data):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as handle:
        handle.write("parent_id,count,sequence\n")
        for row in data:
            handle.write(f"{int(row['pid'])},{int(row['count'])},{row['seq']}\n")


def aggregate_records(records):
    merged = defaultdict(int)
    for pid, count, sequence in records:
        count = int(count)
        if count > 0:
            merged[(int(pid), str(sequence))] += count
    result = np.empty(len(merged), dtype=POOL_DTYPE)
    for index, ((pid, sequence), count) in enumerate(merged.items()):
        result[index] = (pid, count, sequence)
    return result


def safe_sum(values):
    return sum(int(v) for v in values)


def aqueous_depurination_rate(pH, temp_C):
    """
    Legacy per-purine depurination formula
    """
    T = 273.15 + float(temp_C)
    if T <= 0.0:
        raise ValueError("Temperature must be above absolute zero")
    if pH < 2.5:
        lgk = 14.6 - 0.707 * pH - (5.63e3 / T)
    else:
        lgk = 16.5 - 0.982 * pH - (5.85e3 / T)
    return 10.0 ** lgk


def parse_encap_damage_factor(value, model, temp_C):
    text = str(value).strip().lower()
    if text == "bulk_ratio":
        kd = model.k(temp_C, False)
        ke = model.k(temp_C, True)
        if kd <= 0.0:
            return 1.0, "bulk_ratio"
        return float(np.clip(ke / kd, 0.0, 1.0)), "bulk_ratio"
    try:
        factor = float(text)
    except ValueError as exc:
        raise ValueError("--encap_damage_factor must be bulk_ratio or a numeric value 0-1") from exc
    if not 0.0 <= factor <= 1.0:
        raise ValueError("--encap_damage_factor numeric value must be between 0 and 1")
    return factor, "numeric"


def p_damage_per_molecule(sequence, p_site):
    n_purines = sequence.count("A") + sequence.count("G")
    if n_purines == 0 or p_site <= 0.0:
        return 0.0
    return float(-math.expm1(n_purines * math.log1p(-p_site)))


def sample_positive_binomial(n, p, rng):
    """
    Draw K ~ Binomial(n, p) conditioned on K >= 1.

    This avoids inefficient rejection sampling when p is extremely small,
    which is important for mild storage conditions such as 5 C.
    """
    n = int(n)
    p = float(p)

    if n <= 0 or p <= 0.0:
        return 0
    if p >= 1.0:
        return n

    ks = np.arange(1, n + 1, dtype=np.int64)

    log_p = math.log(p)
    log_q = math.log1p(-p)

    log_weights = np.array([
        (
            math.lgamma(n + 1)
            - math.lgamma(int(k) + 1)
            - math.lgamma(n - int(k) + 1)
            + int(k) * log_p
            + (n - int(k)) * log_q
        )
        for k in ks
    ], dtype=np.float64)

    # Stable normalization.
    log_weights -= np.max(log_weights)
    probs = np.exp(log_weights)
    probs /= probs.sum()

    return int(rng.choice(ks, p=probs))


def make_depurination_variant(sequence, p_site, rng):
    """
    Generate one representative damaged molecule.

    A damaged molecule is already known to contain at least one lesion, so the
    lesion count is sampled from Binomial(n_purines, p_site) conditioned on >= 1.

    Distinct purine (A/G) positions are then selected without replacement and
    represented as one-base deletions.

    NOTE:
    Chemically, depurination creates AP sites. Representing each AP site as a
    deletion is a sequence-level simulator approximation for downstream testing.
    """
    purine_positions = [
        index for index, base in enumerate(sequence)
        if base in ("A", "G")
    ]

    n_purines = len(purine_positions)
    if n_purines == 0 or p_site <= 0.0:
        return None, 0

    lesion_count = sample_positive_binomial(
        n=n_purines,
        p=p_site,
        rng=rng,
    )

    lesion_count = max(1, min(int(lesion_count), n_purines))

    chosen_positions = rng.choice(
        np.asarray(purine_positions, dtype=np.int64),
        size=lesion_count,
        replace=False,
    )

    variant = sequence
    for pos in sorted((int(x) for x in chosen_positions), reverse=True):
        variant = variant[:pos] + variant[pos + 1:]

    if variant == sequence:
        return None, 0

    return variant, lesion_count


def build_damage_variants(
    data,
    intact_counts,
    damaged_counts,
    p_site,
    variant_cap,
    rng,
):

    records = []
    representative_variant_records = 0
    fallback_copies = 0
    damaged_parent_populations = 0

    representative_lesion_counts = []
    weighted_lesion_copies = 0

    for row, intact, damaged in zip(data, intact_counts, damaged_counts):
        pid = int(row["pid"])
        seq = str(row["seq"])
        intact = int(intact)
        damaged = int(damaged)

        if intact > 0:
            records.append((pid, intact, seq))

        if damaged <= 0:
            continue

        damaged_parent_populations += 1

        n_purines = seq.count("A") + seq.count("G")
        if n_purines <= 0 or p_site <= 0.0:
            records.append((pid, damaged, seq))
            fallback_copies += damaged
            continue

        n_representatives = min(int(variant_cap), damaged)
        if n_representatives <= 0:
            records.append((pid, damaged, seq))
            fallback_copies += damaged
            continue

        sampled = []
        for _ in range(n_representatives):
            variant, lesion_count = make_depurination_variant(
                sequence=seq,
                p_site=p_site,
                rng=rng,
            )
            if variant is not None and lesion_count > 0:
                sampled.append((variant, int(lesion_count)))

        if not sampled:
            records.append((pid, damaged, seq))
            fallback_copies += damaged
            continue


        portions = rng.multinomial(
            damaged,
            [1.0 / len(sampled)] * len(sampled),
        )

        for (variant, lesion_count), portion in zip(sampled, portions):
            portion = int(portion)
            if portion <= 0:
                continue

            records.append((pid, portion, variant))
            representative_variant_records += 1
            representative_lesion_counts.append(int(lesion_count))
            weighted_lesion_copies += portion * int(lesion_count)

    aggregated = aggregate_records(records)

    if representative_lesion_counts:
        lesion_min = min(representative_lesion_counts)
        lesion_max = max(representative_lesion_counts)
    else:
        lesion_min = 0
        lesion_max = 0

    represented_damaged_copies = safe_sum(damaged_counts) - fallback_copies
    weighted_mean_lesions = (
        weighted_lesion_copies / represented_damaged_copies
        if represented_damaged_copies > 0
        else 0.0
    )

    return (
        aggregated,
        representative_variant_records,
        fallback_copies,
        damaged_parent_populations,
        lesion_min,
        weighted_mean_lesions,
        lesion_max,
    )

def realize_survival(input_counts, survival_probability, mode, rng):
    p = float(np.clip(survival_probability, 0.0, 1.0))
    if mode == "expected":
        return np.rint(input_counts.astype(np.float64) * p).astype(np.int64)
    return rng.binomial(input_counts, p).astype(np.int64)


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    storage_dir = os.path.join(os.getcwd(), "storage")
    os.makedirs(storage_dir, exist_ok=True)

    input_path = args.in_file
    output_path = resolve_output_path(storage_dir, args.out_file)
    xlsx_path = resolve_existing_file(args.arrhenius_xlsx)

    temp_C = float(args.temp)
    pH = float(args.ph)
    weeks = float(args.week)
    rh = float(args.rh)
    encapsulated = args.encap == "1"
    tape_name = "E-DNA (ZIF-encapsulated)" if encapsulated else "D-DNA (decapsulated)"

    if args.c is not None:
        storage_scale = float(args.c)
        scale_source = "custom --c"
    else:
        storage_scale = STORAGE_SCALE_PRESETS[args.mut]
        scale_source = f"preset --mut {args.mut}"

    model = Arrhenius_decay.from_xlsx(xlsx_path)
    reference_k = float(model.k(temp_C, encapsulated))
    effective_k = reference_k * storage_scale
    t_seconds = weeks * SEC_PER_WEEK
    remaining_fraction = math.exp(-effective_k * t_seconds) if effective_k > 0 else 1.0
    remaining_fraction = float(np.clip(remaining_fraction, 0.0, 1.0))

    data = read_pool(input_path)
    input_counts = data["count"].astype(np.int64, copy=True)
    surviving_counts = realize_survival(
        input_counts, remaining_fraction, args.survival_mode, rng
    )

    # Optional, explicitly non-paper-calibrated sequence-lesion model.
    p_site = 0.0
    p_damage = np.zeros(len(data), dtype=np.float64)
    damaged_counts = np.zeros(len(data), dtype=np.int64)
    intact_counts = surviving_counts.copy()
    damage_factor = None
    damage_factor_source = None

    if args.damage_model == "depurination" and storage_scale > 0.0 and weeks > 0.0:
        dep_k = aqueous_depurination_rate(pH, temp_C)
        if encapsulated:
            damage_factor, damage_factor_source = parse_encap_damage_factor(
                args.encap_damage_factor, model, temp_C
            )
            dep_k *= damage_factor
        dep_k *= storage_scale
        p_site = float(-math.expm1(-dep_k * t_seconds))
        p_site = float(np.clip(p_site, 0.0, 1.0))

        p_damage = np.fromiter(
            (p_damage_per_molecule(seq, p_site) for seq in data["seq"]),
            dtype=np.float64,
            count=len(data),
        )
        damaged_counts = rng.binomial(surviving_counts, p_damage).astype(np.int64)
        intact_counts = surviving_counts - damaged_counts

    # Output population.
    representative_variant_records = 0
    damaged_parent_populations = 0
    fallback_copies = 0
    lesion_dropout_total = 0
    representative_lesion_min = 0
    representative_lesion_mean = 0.0
    representative_lesion_max = 0

    if args.damage_model == "depurination" and args.damage_fate == "variants":
        (
            final_data,
            representative_variant_records,
            fallback_copies,
            damaged_parent_populations,
            representative_lesion_min,
            representative_lesion_mean,
            representative_lesion_max,
        ) = build_damage_variants(
            data,
            intact_counts,
            damaged_counts,
            p_site,
            args.variant_cap,
            rng,
        )
        assert safe_sum(final_data["count"]) == safe_sum(surviving_counts)
    elif args.damage_model == "depurination" and args.damage_fate == "dropout":
        final_data = aggregate_records([
            (int(row["pid"]), int(count), str(row["seq"]))
            for row, count in zip(data, intact_counts)
            if int(count) > 0
        ])
        lesion_dropout_total = safe_sum(damaged_counts)
        assert safe_sum(final_data["count"]) + lesion_dropout_total == safe_sum(surviving_counts)
    else:
        final_data = aggregate_records([
            (int(row["pid"]), int(count), str(row["seq"]))
            for row, count in zip(data, surviving_counts)
            if int(count) > 0
        ])
        assert safe_sum(final_data["count"]) == safe_sum(surviving_counts)

    write_pool(output_path, final_data)

    input_total = safe_sum(input_counts)
    survived_total = safe_sum(surviving_counts)
    final_total = safe_sum(final_data["count"])
    bulk_lost = input_total - survived_total
    parent_dropouts = int(np.count_nonzero(surviving_counts == 0))
    damaged_total = safe_sum(damaged_counts)
    damaged_fraction = (damaged_total / survived_total) if survived_total else 0.0


    print("\n--- STORAGE CONFIGURATION ---")
    print(f"Input:                         {input_path}")
    print(f"Output:                        {output_path}")
    print(f"RawData workbook:              {xlsx_path}")
    print(f"Tape state:                    {tape_name}")
    print(f"Temperature:                   {temp_C:g} C")
    print(f"Relative humidity:             {rh:g}%")
    print(f"Storage duration:              {weeks:g} weeks")
    print(f"Storage-rate scale:            {storage_scale:g}x ({scale_source})")
    print(f"Reference k(T):                {reference_k:.8g} s^-1")
    print(f"Effective k(T):                {effective_k:.8g} s^-1")
    print(f"Bulk remaining fraction:       {remaining_fraction:.10%}")
    print(f"Survival realization:          {args.survival_mode}")
    print(f"Sequence-damage model:         {args.damage_model}")
    if args.damage_model == "depurination":
        print(f"pH (lesion model only):        {pH:g}")
        print(f"P(depurination/purine):        {p_site:.10%}")
        if encapsulated:
            print(f"Encap lesion factor:           {damage_factor} ({damage_factor_source})")
        print(f"Damage fate:                   {args.damage_fate}")
    if abs(rh - 50.0) > 1e-12:
        print("WARNING: empirical Arrhenius kinetics were measured at 50% RH;")
        print("         this script does not apply an unvalidated RH correction.")

    print("--- STORAGE POPULATION AUDIT ---")
    print(f"Input records:                  {len(data)}")
    print(f"Input molecules:                {input_total}")
    print(f"Surviving molecules:            {survived_total}")
    print(f"Bulk molecules lost:            {bulk_lost}")
    print(f"Bulk loss fraction:             {(bulk_lost / input_total if input_total else 0):.8%}")
    print(f"Parent oligo dropouts:          {parent_dropouts}")
    print(f"Damaged molecules:              {damaged_total}")
    if args.damage_model == "depurination":
        print(f"Damaged fraction of survivors:  {damaged_fraction:.8%}")
        print(f"Parents with damaged copies:     {damaged_parent_populations}")
        print(f"Representative variant records: {representative_variant_records}")
        print(f"Representative variant cap:     {args.variant_cap} per parent")
        if args.damage_fate == "variants":
            print(
                "Representative lesions/variant: "
                f"{representative_lesion_min}/"
                f"{representative_lesion_mean:.4f}/"
                f"{representative_lesion_max} (min/weighted-mean/max)"
            )
        print(f"Variant fallback copies:        {fallback_copies}")
        print(f"Lesion-induced dropout copies:  {lesion_dropout_total}")
    print(f"Final output molecules:         {final_total}")
    print(f"Final output records:           {len(final_data)}")
    print("--------------------------------\n")

    print("storage.py run completed")


if __name__ == "__main__":
    main()