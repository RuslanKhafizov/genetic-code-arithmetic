#!/usr/bin/env python3
"""Simple Monte Carlo test for the genetic-code arithmetic statistics.

One trial uniformly relabels the 20 amino-acid blocks with the 20 real
amino-acid (P, N) profiles.  The code architecture stays fixed.

S1--S3 are the original sense-pool count statistics.  S4sense measures the
mean distance to the 37-lattice without assigning service codons.  S4a counts
the relabelled product at ATG and zero at stops.  S4b rederives the minimal
Key 1 for every random code by the same balance and divisibility rules.
"""

import argparse
import csv
import math
import random
from collections import defaultdict
from pathlib import Path


MODULUS = 37
STOPS = {"TAA", "TAG", "TGA"}
SERVICE = STOPS | {"ATG"}


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def dot(row, values):
    return sum(weight * values[index] for index, weight in row)


def distance_to_37(value):
    residue = value % MODULUS
    return min(residue, MODULUS - residue)


def wilson_interval(hits, trials):
    """95% interval for the Monte Carlo tail probability."""
    z = 1.959963984540054
    rate = hits / trials
    denominator = 1.0 + z * z / trials
    centre = (rate + z * z / (2.0 * trials)) / denominator
    half_width = (
        z
        * math.sqrt(rate * (1.0 - rate) / trials + z * z / (4.0 * trials * trials))
        / denominator
    )
    lower = 0.0 if hits == 0 else max(0.0, centre - half_width)
    upper = 1.0 if hits == trials else min(1.0, centre + half_width)
    return lower, upper


def parse_args():
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=script_dir)
    parser.add_argument("--output", type=Path, default=script_dir / "mc_significance.csv")
    parser.add_argument("--trials", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.trials <= 0:
        raise ValueError("--trials must be positive")

    amino_rows = read_csv(args.data_dir / "amino_acids_nucleons.csv")
    codon_rows = read_csv(args.data_dir / "genetic_code_codons.csv")
    group_rows = read_csv(args.data_dir / "codon_groups.csv")

    profiles = {
        row["Amino_Acid"]: (int(row["Protons"]), int(row["Neutrons"]))
        for row in amino_rows
    }
    amino_acids = sorted(profiles)
    if len(amino_acids) != 20:
        raise ValueError(f"Expected 20 amino acids, found {len(amino_acids)}")

    block_index = {amino_acid: index for index, amino_acid in enumerate(amino_acids)}
    proton_values = [profiles[amino_acid][0] for amino_acid in amino_acids]
    neutron_values = [profiles[amino_acid][1] for amino_acid in amino_acids]
    met_block = block_index["Methionine"]

    codon_to_block = {}
    for row in codon_rows:
        for codon in row["Codons"].split(";"):
            codon_to_block[codon.strip()] = row["Product"]
    if len(codon_to_block) != 64:
        raise ValueError(f"Expected 64 codons, found {len(codon_to_block)}")

    def codons_of(group_row):
        return [codon.strip() for codon in group_row["Codon_List"].split(";")]

    def weight_row(codons):
        counts = [0] * len(amino_acids)
        for codon in codons:
            if codon in SERVICE:
                continue
            counts[block_index[codon_to_block[codon]]] += 1
        return [(index, count) for index, count in enumerate(counts) if count]

    if len(group_rows) != 33:
        raise ValueError(f"Expected 33 codon groups, found {len(group_rows)}")

    groups = []
    rows_by_name = defaultdict(list)
    for group_row in group_rows:
        codons = codons_of(group_row)
        groups.append(
            (
                weight_row(codons),
                int("ATG" in codons),
                sum(codon in STOPS for codon in codons),
            )
        )
        rows_by_name[group_row["Group_Name"]].append(group_row)

    def unique_weight(name):
        matches = rows_by_name[name]
        if len(matches) != 1:
            raise ValueError(f"Expected one group named {name!r}, found {len(matches)}")
        return weight_row(codons_of(matches[0]))

    headline_names = [
        "ALL: {C, G, A, T}",
        "Keto: {G, T}",
        "Amino: {A, C}",
        "Strong: {C, G}",
        "Weak: {A, T}",
        "Purine: {A, G}",
        "Pyrimidine: {C, T}",
        "Octet I: {C, G, A, T}",
        "Octet II: {C, G, A, T}",
    ]
    headline_rows = [unique_weight(name) for name in headline_names]
    all_sense = headline_rows[0]
    keto_sense = headline_rows[1]
    amino_sense = headline_rows[2]
    octet_i = headline_rows[7]

    independent_rows = []
    for group_row in group_rows:
        codons = codons_of(group_row)
        if not (set(codons) & SERVICE):
            independent_rows.append(weight_row(codons))
    if len(independent_rows) != 17:
        raise ValueError(
            f"Expected 17 service-free groups for S3, found {len(independent_rows)}"
        )

    def divisibility_count(rows, protons, neutrons):
        count = 0
        for row in rows:
            p_total = dot(row, protons)
            n_total = dot(row, neutrons)
            for value in (p_total, n_total, p_total + n_total, p_total - n_total):
                count += value % MODULUS == 0
        return count

    def key_pair(total, difference):
        """Return the minimal nonnegative (start, stop) pair for Key 1."""
        start = ((-total - 3 * difference) * pow(4, -1, MODULUS)) % MODULUS
        while start + difference < 0:
            start += MODULUS
        return start, start + difference

    def lattice_scores(protons, neutrons):
        p_difference = dot(keto_sense, protons) - dot(amino_sense, protons)
        n_difference = dot(keto_sense, neutrons) - dot(amino_sense, neutrons)
        p_start, p_stop = key_pair(dot(all_sense, protons), p_difference)
        n_start, n_stop = key_pair(dot(all_sense, neutrons), n_difference)

        sense_distance = 0
        key0_distance = 0
        key1_distance = 0
        for row, has_atg, stop_count in groups:
            base_p = dot(row, protons)
            base_n = dot(row, neutrons)

            for value in (base_p, base_n, base_p + base_n, base_p - base_n):
                sense_distance += distance_to_37(value)

            key0_p = base_p + has_atg * protons[met_block]
            key0_n = base_n + has_atg * neutrons[met_block]
            for value in (key0_p, key0_n, key0_p + key0_n, key0_p - key0_n):
                key0_distance += distance_to_37(value)

            key1_p = base_p + has_atg * p_start + stop_count * p_stop
            key1_n = base_n + has_atg * n_start + stop_count * n_stop
            for value in (key1_p, key1_n, key1_p + key1_n, key1_p - key1_n):
                key1_distance += distance_to_37(value)

        denominator = len(groups) * 4
        return (
            sense_distance / denominator,
            key0_distance / denominator,
            key1_distance / denominator,
            (p_start, p_stop, n_start, n_stop),
        )

    observed_s4sense, observed_s4a, observed_s4b, observed_key = lattice_scores(
        proton_values, neutron_values
    )
    observed = {
        "S1": int(
            dot(all_sense, neutron_values) % MODULUS == 0
            and (dot(octet_i, proton_values) + dot(octet_i, neutron_values))
            % MODULUS
            == 0
        ),
        "S2": divisibility_count(headline_rows, proton_values, neutron_values),
        "S3": divisibility_count(independent_rows, proton_values, neutron_values),
        "S4sense": observed_s4sense,
        "S4a": observed_s4a,
        "S4b": observed_s4b,
    }

    expected_key = (1, 37, 0, 37)
    if observed_key != expected_key:
        raise AssertionError(f"Expected standard Key 1 {expected_key}, found {observed_key}")
    if observed["S1"] != 1 or observed["S2"] != 13 or observed["S3"] != 14:
        raise AssertionError(f"Unexpected observed count statistics: {observed}")

    descriptions = {
        "S1": "both 37 anchors are divisible",
        "S2": "divisible count over 9 headline sense-pool groups (higher is stronger)",
        "S3": "divisible count over 17 service-free groups (higher is stronger)",
        "S4sense": "mean distance over 33 sense-pool groups (lower is stronger)",
        "S4a": "mean distance under Key 0 with the relabelled ATG product (lower is stronger)",
        "S4b": "mean distance after rederiving Key 1 for each code (lower is stronger)",
    }
    tail_direction = {
        "S1": "event",
        "S2": ">=",
        "S3": ">=",
        "S4sense": "<=",
        "S4a": "<=",
        "S4b": "<=",
    }

    tail_hits = {name: 0 for name in observed}
    null_sums = {name: 0.0 for name in observed}
    rng = random.Random(args.seed)
    identity = list(range(len(amino_acids)))
    progress_step = max(1, args.trials // 10)

    for trial in range(1, args.trials + 1):
        permutation = identity[:]
        rng.shuffle(permutation)
        protons = [proton_values[index] for index in permutation]
        neutrons = [neutron_values[index] for index in permutation]

        values = {}
        values["S1"] = int(
            dot(all_sense, neutrons) % MODULUS == 0
            and (dot(octet_i, protons) + dot(octet_i, neutrons)) % MODULUS == 0
        )
        values["S2"] = divisibility_count(headline_rows, protons, neutrons)
        values["S3"] = divisibility_count(independent_rows, protons, neutrons)
        values["S4sense"], values["S4a"], values["S4b"], _ = lattice_scores(
            protons, neutrons
        )

        for name, value in values.items():
            null_sums[name] += value
        tail_hits["S1"] += values["S1"] == 1
        tail_hits["S2"] += values["S2"] >= observed["S2"]
        tail_hits["S3"] += values["S3"] >= observed["S3"]
        tail_hits["S4sense"] += values["S4sense"] <= observed["S4sense"]
        tail_hits["S4a"] += values["S4a"] <= observed["S4a"]
        tail_hits["S4b"] += values["S4b"] <= observed["S4b"]

        if trial % progress_step == 0 or trial == args.trials:
            print(f"Completed {trial:,} / {args.trials:,} trials", flush=True)

    rows = []
    for name in observed:
        hits = tail_hits[name]
        p_value = (hits + 1) / (args.trials + 1)
        ci_low, ci_high = wilson_interval(hits, args.trials)
        rows.append(
            {
                "Statistic": name,
                "Description": descriptions[name],
                "Tail": tail_direction[name],
                "Observed": observed[name],
                "Null_Mean": null_sums[name] / args.trials,
                "Tail_Count": hits,
                "P_MC": p_value,
                "CI_Low": ci_low,
                "CI_High": ci_high,
                "N_Trials": args.trials,
                "Seed": args.seed,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Standard Key 1 (Pstart, Pstop, Nstart, Nstop): {observed_key}")
    print("Statistic  Observed  Null mean  Tail count  P_MC       95% MC interval")
    for row in rows:
        print(
            f"{row['Statistic']:9s} "
            f"{float(row['Observed']):8.4f}  "
            f"{row['Null_Mean']:9.4f}  "
            f"{row['Tail_Count']:10d}  "
            f"{row['P_MC']:.6g}  "
            f"[{row['CI_Low']:.6g}, {row['CI_High']:.6g}]"
        )
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
