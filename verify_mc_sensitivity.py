#!/usr/bin/env python3
"""Independent standard-library check of mc_sensitivity.csv.

Rebuild groups from codon rules, evaluate eight position/section bins,
and sample with a different seed. No import from mc_sensitivity.py and
no reading of codon_groups.csv. Optional --output saves the audit table.
"""
import argparse
import csv
from fractions import Fraction
import itertools
import math
from pathlib import Path
import random

SERVICE = {"ATG", "TAA", "TAG", "TGA"}
MODELS = ("relabel20", "fixed_atg_met19")
STATS = ("S2", "S4sense", "S4pos")


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def build_problem(folder):
    amino = read_rows(folder / "amino_acids_nucleons.csv")
    names = sorted(row["Amino_Acid"] for row in amino)
    if len(set(names)) != 20 or len(names) != 20:
        raise ValueError("Expected 20 unique profiles")
    lookup = {row["Amino_Acid"]: (int(row["Protons"]), int(row["Neutrons"])) for row in amino}
    assignment = {}
    for row in read_rows(folder / "genetic_code_codons.csv"):
        for raw_codon in row["Codons"].split(";"):
            codon = raw_codon.strip()
            if codon in assignment:
                raise ValueError("Duplicate codon")
            assignment[codon] = row["Product"]
    codons = ["".join(chars) for chars in itertools.product("ACGT", repeat=3)]
    if set(assignment) != set(codons) or assignment["ATG"] != "Methionine":
        raise ValueError("Invalid standard-code assignments")
    if {c for c in codons if assignment[c] == "STOP"} != SERVICE - {"ATG"}:
        raise ValueError("Invalid stop assignments")
    if set(assignment.values()) - {"STOP"} != set(names):
        raise ValueError("Profiles and codon assignments disagree")
    families = [a + b for a in "ACGT" for b in "ACGT"]
    octet = {family for family in families
             if len({assignment[family + c] for c in "ACGT"}) == 1
             and assignment[family + "A"] != "STOP"}
    if len(octet) != 8:
        raise ValueError("Expected eight fourfold families")
    bins = []
    bin_tags = []
    for section in ("I", "II"):
        for letter in "ACGT":
            selected = [c for c in codons if c[2] == letter
                        and ((c[:2] in octet) == (section == "I")) and c not in SERVICE]
            bins.append(tuple(names.index(assignment[c]) for c in selected))
            bin_tags.append((section, letter))
    groups = []
    for section in ("all", "I", "II"):
        for letters in ("ACGT", "GT", "CG", "AC", "AT", "AG", "CT", "C", "G", "T", "A"):
            members = tuple(i for i, (s, letter) in enumerate(bin_tags)
                            if (section == "all" or s == section) and letter in letters)
            headline = letters == "ACGT" or (section == "all" and len(letters) == 2)
            groups.append((members, headline))
    assert len(groups) == 33 and sum(h for _, h in groups) == 9
    return names, [lookup[name] for name in names], names.index("Methionine"), bins, groups


def evaluate(p, n, bins, groups):
    bin_p = [sum(p[i] for i in members) for members in bins]
    bin_n = [sum(n[i] for i in members) for members in bins]
    s2 = total = nonzero = 0
    for members, headline in groups:
        gp = sum(bin_p[i] for i in members)
        gn = sum(bin_n[i] for i in members)
        for quantity in (gp + gn, gp, gn, gp - gn):
            # Nearest multiple, with the balanced residue in [-18, 18].
            distance = abs((quantity + 18) % 37 - 18)
            total += distance
            nonzero += distance != 0
            if headline and distance == 0:
                s2 += 1
    return s2, total, nonzero


def statistics_and_tails(raw, observed):
    count, total, positive = raw
    reference_count, reference_total, reference_positive = observed
    mean_positive = Fraction(total, positive) if positive else Fraction(0)
    reference_positive_mean = (Fraction(reference_total, reference_positive)
                               if reference_positive else Fraction(0))
    return ((count, total / 132, float(mean_positive)),
            (count >= reference_count, total <= reference_total,
             mean_positive <= reference_positive_mean))


def sample_assignment(rng, model, fixed):
    targets = list(range(20))
    if model == "fixed_atg_met19":
        targets.remove(fixed)
    elif model != "relabel20":
        raise ValueError("Unknown model")
    labels = targets.copy()
    rng.shuffle(labels)
    mapping = dict(zip(targets, labels))
    return [mapping.get(i, i) for i in range(20)]


def simulate(profiles, fixed, bins, groups, model, trials, seed):
    p0 = [p for p, _ in profiles]
    n0 = [n for _, n in profiles]
    observed = evaluate(p0, n0, bins, groups)
    if observed != (13, 636, 104):
        raise AssertionError("Unexpected observed statistics")
    obs_values, _ = statistics_and_tails(observed, observed)
    means, m2, hits = [0.0] * 3, [0.0] * 3, [0] * 3
    empty = 0
    rng = random.Random(seed)
    step = max(1, trials // 10)
    for trial in range(1, trials + 1):
        order = sample_assignment(rng, model, fixed)
        raw = evaluate([p0[i] for i in order], [n0[i] for i in order], bins, groups)
        values, tail = statistics_and_tails(raw, observed)
        empty += raw[2] == 0
        for j in range(3):
            delta = values[j] - means[j]
            means[j] += delta / trial
            m2[j] += delta * (values[j] - means[j])
            hits[j] += tail[j]
        if trial % step == 0 or trial == trials:
            print(f"Verify {model}: {trial:,}/{trials:,}", flush=True)
    return [{"Statistic": name, "Observed": obs_values[j], "Null_Mean": means[j],
             "Null_SD": math.sqrt(max(0.0, m2[j] / (trials - 1))), "Tail_Count": hits[j],
             "P_MC": (hits[j] + 1) / (trials + 1), "All_Zero_Trials": empty}
            for j, name in enumerate(STATS)]


def standard_errors_apart(first, second, se):
    if se == 0:
        return 0.0 if first == second else math.inf
    return abs(first - second) / se


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=here)
    parser.add_argument("--reference", type=Path, default=here / "mc_sensitivity.csv")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--trials", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=12345)
    args = parser.parse_args()
    if args.trials < 2:
        parser.error("--trials must be at least 2")
    reference_rows = read_rows(args.reference)
    reference = {(r["Null_Model"], r["Statistic"]): r for r in reference_rows}
    expected = set(itertools.product(MODELS, STATS))
    if len(reference_rows) != 6 or set(reference) != expected:
        raise ValueError("Reference must contain all six model/statistic combinations")
    if any(int(r["N_Trials"]) < 2 for r in reference_rows):
        raise ValueError("Invalid reference trial count")
    if any(int(r["Seed"]) == args.seed for r in reference_rows):
        parser.error("Use a verification seed different from the reference seed")
    _, profiles, fixed, bins, groups = build_problem(args.data_dir)
    output, failures = [], []
    for model in MODELS:
        results = simulate(profiles, fixed, bins, groups, model, args.trials, args.seed)
        for result in results:
            name = result["Statistic"]
            ref = reference[(model, name)]
            if not math.isclose(result["Observed"], float(ref["Observed"]), rel_tol=0, abs_tol=1e-12):
                raise AssertionError("Observed value mismatch: " + model + "/" + name)
            if ref["Tail"] != (">=" if name == "S2" else "<="):
                raise AssertionError("Reference tail direction mismatch")
            count, nref = int(ref["Tail_Count"]), int(ref["N_Trials"])
            if not 0 <= count <= nref or not math.isclose(float(ref["P_MC"]), (count + 1) / (nref + 1), rel_tol=1e-12):
                raise AssertionError("Invalid reference p or tail count")
            pref, pver = float(ref["P_MC"]), result["P_MC"]
            tail_se = math.sqrt(pref * (1 - pref) / nref + pver * (1 - pver) / args.trials)
            mean_se = math.sqrt(float(ref["Null_SD"]) ** 2 / nref + result["Null_SD"] ** 2 / args.trials)
            pz = standard_errors_apart(pref, pver, tail_se)
            mz = standard_errors_apart(float(ref["Null_Mean"]), result["Null_Mean"], mean_se)
            passed = pz <= 5 and mz <= 5
            if not passed:
                failures.append(model + "/" + name)
            row = {"Null_Model": model, **result, "N_Trials": args.trials, "Seed": args.seed,
                   "Reference_P_MC": pref, "P_MC_SE_Difference": pz,
                   "Reference_Null_Mean": float(ref["Null_Mean"]),
                   "Mean_MC_SE_Difference": mz, "Passed": passed}
            output.append(row)
            print(f"{model}/{name}: tail {pz:.3f} MC SE, mean {mz:.3f} MC SE; passed={passed}", flush=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(output[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(output)
    if failures:
        raise AssertionError("Independent checks failed: " + ", ".join(failures))
    print("Independent sensitivity verification passed.")


if __name__ == "__main__":
    main()
