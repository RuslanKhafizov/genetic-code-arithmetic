#!/usr/bin/env python3
"""Post hoc sensitivity checks for genetic-code arithmetic (standard library).

Compare S2, S4sense and the mean positive distance S4pos under relabel20
and fixed_atg_met19. The six main statistics in mc_significance.py are
unchanged. Repeated group entries retain their original multiplicities.
"""
import argparse
from collections import Counter
import csv
import math
from pathlib import Path
import random

MODULUS = 37
SERVICE = {"ATG", "TAA", "TAG", "TGA"}
MODELS = ("relabel20", "fixed_atg_met19")
STATISTICS = ("S2", "S4sense", "S4pos")
HEADLINES = {
    "ALL: {C, G, A, T}", "Keto: {G, T}", "Amino: {A, C}",
    "Strong: {C, G}", "Weak: {A, T}", "Purine: {A, G}",
    "Pyrimidine: {C, T}", "Octet I: {C, G, A, T}",
    "Octet II: {C, G, A, T}",
}
FIELDS = ("Null_Model", "Statistic", "Tail", "Observed", "Null_Mean",
          "Null_SD", "Tail_Count", "P_MC", "CI_Low", "CI_High",
          "N_Trials", "Seed", "Observed_Positive_Count", "All_Zero_Trials")


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def load_problem(data_dir):
    amino_rows = read_csv(data_dir / "amino_acids_nucleons.csv")
    names = sorted(row["Amino_Acid"] for row in amino_rows)
    if len(names) != 20 or len(set(names)) != 20:
        raise ValueError("Expected 20 distinct amino acids")
    profiles = {row["Amino_Acid"]: (int(row["Protons"]), int(row["Neutrons"]))
                for row in amino_rows}
    index = {name: i for i, name in enumerate(names)}
    codon_to_product = {}
    for row in read_csv(data_dir / "genetic_code_codons.csv"):
        for codon in row["Codons"].split(";"):
            codon = codon.strip()
            if codon in codon_to_product:
                raise ValueError("Duplicate codon: " + codon)
            codon_to_product[codon] = row["Product"]
    codons = {a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"}
    if set(codon_to_product) != codons:
        raise ValueError("Expected exactly the 64 DNA codons")
    if {c for c, product in codon_to_product.items() if product == "STOP"} != SERVICE - {"ATG"}:
        raise ValueError("Unexpected stop codons")
    if codon_to_product["ATG"] != "Methionine":
        raise ValueError("The reference ATG block must encode methionine")
    if set(codon_to_product.values()) - {"STOP"} != set(names):
        raise ValueError("Codon products and amino-acid profiles disagree")

    # Identical weight rows are evaluated once, with their full multiplicity.
    # This is an arithmetic optimisation, not a deduplicated statistic.
    all_rows, headline_rows = Counter(), Counter()
    seen_headlines = Counter()
    group_rows = read_csv(data_dir / "codon_groups.csv")
    if len(group_rows) != 33:
        raise ValueError("Expected 33 groups")
    for group in group_rows:
        group_codons = [c.strip() for c in group["Codon_List"].split(";")]
        if len(group_codons) != len(set(group_codons)) or not set(group_codons) <= codons:
            raise ValueError("Invalid codon list in " + group["Group_Name"])
        counts = Counter(index[codon_to_product[c]] for c in group_codons if c not in SERVICE)
        row = tuple(sorted(counts.items()))
        all_rows[row] += 1
        if group["Group_Name"] in HEADLINES:
            headline_rows[row] += 1
            seen_headlines[group["Group_Name"]] += 1
    if seen_headlines != Counter({name: 1 for name in HEADLINES}):
        raise ValueError("Expected each of the nine headline groups exactly once")
    compiled = [(row, multiplicity, headline_rows[row]) for row, multiplicity in all_rows.items()]
    return names, [profiles[name] for name in names], index["Methionine"], compiled


def score(protons, neutrons, compiled):
    """Return S2, the integer sum of 132 distances, and the positive count."""
    s2 = distance_sum = positive_count = 0
    for row, multiplicity, headline_multiplicity in compiled:
        p = sum(weight * protons[i] for i, weight in row)
        n = sum(weight * neutrons[i] for i, weight in row)
        distances = []
        for value in (p, n, p + n, p - n):
            r = value % MODULUS
            distances.append(min(r, MODULUS - r))
        zeros = distances.count(0)
        s2 += headline_multiplicity * zeros
        distance_sum += multiplicity * sum(distances)
        positive_count += multiplicity * (4 - zeros)
    return s2, distance_sum, positive_count


def values_and_hits(raw, observed):
    count, total, positive = raw
    observed_count, observed_total, observed_positive = observed
    # Explicit convention: all distances zero => S4pos=0. Never drop a trial.
    pos_value = total / positive if positive else 0.0
    pos_hit = (positive == 0 or
               (observed_positive > 0 and total * observed_positive <= observed_total * positive))
    return ((count, total / 132, pos_value),
            (count >= observed_count, total <= observed_total, pos_hit))


def permutation(rng, model, fixed_index):
    result = list(range(20))
    if model == "relabel20":
        rng.shuffle(result)
    elif model == "fixed_atg_met19":
        active = [i for i in result if i != fixed_index]
        shuffled = active[:]
        rng.shuffle(shuffled)
        for target, source in zip(active, shuffled):
            result[target] = source
    else:
        raise ValueError("Unknown model: " + model)
    return result


def wilson_interval(hits, trials):
    z = 1.959963984540054
    rate = hits / trials
    denom = 1 + z * z / trials
    centre = (rate + z * z / (2 * trials)) / denom
    half = z * math.sqrt(rate * (1 - rate) / trials + z * z / (4 * trials * trials)) / denom
    return (0.0 if hits == 0 else max(0.0, centre - half),
            1.0 if hits == trials else min(1.0, centre + half))


def run_model(profiles, fixed_index, compiled, model, trials, seed):
    original_p = [p for p, _ in profiles]
    original_n = [n for _, n in profiles]
    observed = score(original_p, original_n, compiled)
    if observed != (13, 636, 104):
        raise AssertionError("Unexpected reference scores: " + str(observed))
    observed_values, _ = values_and_hits(observed, observed)
    sums, squares, hits = [0.0] * 3, [0.0] * 3, [0] * 3
    all_zero = 0
    rng = random.Random(seed)
    step = max(1, trials // 10)
    for trial in range(1, trials + 1):
        order = permutation(rng, model, fixed_index)
        raw = score([original_p[i] for i in order], [original_n[i] for i in order], compiled)
        vals, hit = values_and_hits(raw, observed)
        all_zero += raw[2] == 0
        for j in range(3):
            sums[j] += vals[j]
            squares[j] += vals[j] * vals[j]
            hits[j] += hit[j]
        if trial % step == 0 or trial == trials:
            print(f"{model}: {trial:,}/{trials:,}", flush=True)
    rows = []
    for j, statistic in enumerate(STATISTICS):
        lower, upper = wilson_interval(hits[j], trials)
        variance = max(0.0, (squares[j] - sums[j] * sums[j] / trials) / (trials - 1))
        rows.append(dict(zip(FIELDS, (
            model, statistic, ">=" if j == 0 else "<=", observed_values[j],
            sums[j] / trials, math.sqrt(variance), hits[j], (hits[j] + 1) / (trials + 1),
            lower, upper, trials, seed, observed[2], all_zero))))
    return rows


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=here)
    parser.add_argument("--output", type=Path, default=here / "mc_sensitivity.csv")
    parser.add_argument("--trials", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    if args.trials < 2:
        parser.error("--trials must be at least 2")
    _, profiles, fixed_index, compiled = load_problem(args.data_dir)
    rows = []
    for model in MODELS:
        rows.extend(run_model(profiles, fixed_index, compiled, model, args.trials, args.seed))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print("Model                 Statistic  Observed    Null mean    Tail       P_MC")
    for row in rows:
        print(f"{row['Null_Model']:21s} {row['Statistic']:9s} {row['Observed']:10.6f} "
              f"{row['Null_Mean']:11.6f} {row['Tail_Count']:8d} {row['P_MC']:12.6g}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
