#!/usr/bin/env python3
"""Independent check of mc_significance.py.

This script rebuilds the 33 groups from codon rules instead of reading
codon_groups.csv, derives Key 1 by direct search, and uses a different seed.
"""

import argparse
import csv
import math
import random
from collections import defaultdict
from pathlib import Path


MODULUS = 37
BASES = "ACGT"
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


def parse_args():
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=script_dir)
    parser.add_argument("--reference", type=Path, default=script_dir / "mc_significance.csv")
    parser.add_argument("--trials", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=12345)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.trials <= 0:
        raise ValueError("--trials must be positive")

    amino_rows = read_csv(args.data_dir / "amino_acids_nucleons.csv")
    codon_rows = read_csv(args.data_dir / "genetic_code_codons.csv")
    reference_rows = {
        row["Statistic"]: row for row in read_csv(args.reference)
    }

    profiles = {
        row["Amino_Acid"]: (int(row["Protons"]), int(row["Neutrons"]))
        for row in amino_rows
    }
    amino_acids = sorted(profiles)
    block_index = {amino_acid: index for index, amino_acid in enumerate(amino_acids)}
    proton_values = [profiles[amino_acid][0] for amino_acid in amino_acids]
    neutron_values = [profiles[amino_acid][1] for amino_acid in amino_acids]
    met_block = block_index["Methionine"]

    codon_to_block = {}
    for row in codon_rows:
        for codon in row["Codons"].split(";"):
            codon_to_block[codon.strip()] = row["Product"]

    all_codons = [a + b + c for a in BASES for b in BASES for c in BASES]
    boxes = defaultdict(list)
    for codon in all_codons:
        boxes[codon[:2]].append(codon)
    octet_i = {
        codon
        for codons in boxes.values()
        if len({codon_to_block[codon] for codon in codons}) == 1
        and codon_to_block[codons[0]] != "Stop"
        for codon in codons
    }

    def select(letters, section):
        codons = [codon for codon in all_codons if codon[2] in set(letters)]
        if section == "I":
            return [codon for codon in codons if codon in octet_i]
        if section == "II":
            return [codon for codon in codons if codon not in octet_i]
        return codons

    axes = ["GT", "CG", "AC", "AT", "AG", "CT"]
    single_bases = ["C", "G", "T", "A"]
    group_codons = []
    for section in ("all", "I", "II"):
        group_codons.extend(select(axis, section) for axis in axes)
        if section == "all":
            group_codons.append(list(all_codons))
        elif section == "I":
            group_codons.append([codon for codon in all_codons if codon in octet_i])
        else:
            group_codons.append([codon for codon in all_codons if codon not in octet_i])
        group_codons.extend(select(base, section) for base in single_bases)
    if len(group_codons) != 33:
        raise AssertionError(f"Expected 33 rebuilt groups, found {len(group_codons)}")

    def weight_row(codons):
        counts = [0] * len(amino_acids)
        for codon in codons:
            if codon in SERVICE:
                continue
            counts[block_index[codon_to_block[codon]]] += 1
        return [(index, count) for index, count in enumerate(counts) if count]

    groups = [
        (
            weight_row(codons),
            int("ATG" in codons),
            sum(codon in STOPS for codon in codons),
        )
        for codons in group_codons
    ]

    sense_codons = [codon for codon in all_codons if codon not in SERVICE]
    all_sense = weight_row(sense_codons)
    keto_sense = weight_row([codon for codon in sense_codons if codon[2] in "GT"])
    amino_sense = weight_row([codon for codon in sense_codons if codon[2] in "AC"])
    headline_rows = [
        all_sense,
        keto_sense,
        amino_sense,
        weight_row([codon for codon in sense_codons if codon[2] in "CG"]),
        weight_row([codon for codon in sense_codons if codon[2] in "AT"]),
        weight_row([codon for codon in sense_codons if codon[2] in "AG"]),
        weight_row([codon for codon in sense_codons if codon[2] in "CT"]),
        weight_row(sorted(octet_i)),
        weight_row([codon for codon in sense_codons if codon not in octet_i]),
    ]
    octet_i_row = headline_rows[7]
    independent_rows = [
        weight_row(codons) for codons in group_codons if not (set(codons) & SERVICE)
    ]
    if len(independent_rows) != 17:
        raise AssertionError(
            f"Expected 17 service-free groups, found {len(independent_rows)}"
        )

    def divisibility_count(rows, protons, neutrons):
        count = 0
        for row in rows:
            p_total = dot(row, protons)
            n_total = dot(row, neutrons)
            count += p_total % MODULUS == 0
            count += n_total % MODULUS == 0
            count += (p_total + n_total) % MODULUS == 0
            count += (p_total - n_total) % MODULUS == 0
        return count

    def key_pair_by_search(total, difference):
        """Independent Key 1 derivation without using a modular inverse."""
        first_possible_start = max(0, -difference)
        for start in range(first_possible_start, first_possible_start + MODULUS):
            stop = start + difference
            if (total + start + 3 * stop) % MODULUS == 0:
                return start, stop
        raise AssertionError("No Key 1 pair found")

    def lattice_scores(protons, neutrons):
        p_difference = dot(keto_sense, protons) - dot(amino_sense, protons)
        n_difference = dot(keto_sense, neutrons) - dot(amino_sense, neutrons)
        p_start, p_stop = key_pair_by_search(dot(all_sense, protons), p_difference)
        n_start, n_stop = key_pair_by_search(dot(all_sense, neutrons), n_difference)

        distances = [0, 0, 0]
        for row, has_atg, stop_count in groups:
            base_p = dot(row, protons)
            base_n = dot(row, neutrons)
            for value in (base_p, base_n, base_p + base_n, base_p - base_n):
                distances[0] += distance_to_37(value)

            key0_p = base_p + has_atg * protons[met_block]
            key0_n = base_n + has_atg * neutrons[met_block]
            for value in (key0_p, key0_n, key0_p + key0_n, key0_p - key0_n):
                distances[1] += distance_to_37(value)

            key1_p = base_p + has_atg * p_start + stop_count * p_stop
            key1_n = base_n + has_atg * n_start + stop_count * n_stop
            for value in (key1_p, key1_n, key1_p + key1_n, key1_p - key1_n):
                distances[2] += distance_to_37(value)

        denominator = len(groups) * 4
        return (
            distances[0] / denominator,
            distances[1] / denominator,
            distances[2] / denominator,
            (p_start, p_stop, n_start, n_stop),
        )

    observed_s4sense, observed_s4a, observed_s4b, observed_key = lattice_scores(
        proton_values, neutron_values
    )
    observed = {
        "S1": int(
            dot(all_sense, neutron_values) % MODULUS == 0
            and (dot(octet_i_row, proton_values) + dot(octet_i_row, neutron_values))
            % MODULUS
            == 0
        ),
        "S2": divisibility_count(headline_rows, proton_values, neutron_values),
        "S3": divisibility_count(independent_rows, proton_values, neutron_values),
        "S4sense": observed_s4sense,
        "S4a": observed_s4a,
        "S4b": observed_s4b,
    }
    if observed_key != (1, 37, 0, 37):
        raise AssertionError(f"Unexpected standard Key 1: {observed_key}")

    expected_names = set(observed)
    if set(reference_rows) != expected_names:
        raise AssertionError(
            f"Reference statistics are {sorted(reference_rows)}, expected {sorted(expected_names)}"
        )
    for name, value in observed.items():
        reference_value = float(reference_rows[name]["Observed"])
        if not math.isclose(float(value), reference_value, rel_tol=0.0, abs_tol=1e-12):
            raise AssertionError(
                f"Observed {name} differs: verify={value}, reference={reference_value}"
            )

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

        values = {
            "S1": int(
                dot(all_sense, neutrons) % MODULUS == 0
                and (dot(octet_i_row, protons) + dot(octet_i_row, neutrons))
                % MODULUS
                == 0
            ),
            "S2": divisibility_count(headline_rows, protons, neutrons),
            "S3": divisibility_count(independent_rows, protons, neutrons),
        }
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
            print(f"Completed {trial:,} / {args.trials:,} verification trials", flush=True)

    print("Statistic  Reference p  Verify p    Difference in MC standard errors")
    failures = []
    for name in observed:
        reference = reference_rows[name]
        reference_p = float(reference["P_MC"])
        reference_trials = int(reference["N_Trials"])
        verify_p = (tail_hits[name] + 1) / (args.trials + 1)
        standard_error = math.sqrt(
            reference_p * (1.0 - reference_p) / reference_trials
            + verify_p * (1.0 - verify_p) / args.trials
        )
        z_score = abs(reference_p - verify_p) / standard_error if standard_error else 0.0
        print(f"{name:9s}  {reference_p:11.6g}  {verify_p:9.6g}  {z_score:8.3f}")
        if z_score > 5.0:
            failures.append((name, z_score))

    if failures:
        raise AssertionError(f"Monte Carlo results disagree: {failures}")
    print("Independent verification passed.")


if __name__ == "__main__":
    main()
