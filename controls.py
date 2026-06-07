#!/usr/bin/env python3
"""
controls.py — Supplementary specificity controls for the standard genetic code.

This script is NOT part of the core reproducible pipeline (reproduce.py). It is
a standalone, deterministic, exact-integer companion that consumes the same
input datasets and writes two supplementary control datasets:

  * position_axis_analysis.csv
        The three chemical axes (Keto/Amino, Strong/Weak, Purine/Pyrimidine)
        applied independently to the FIRST, SECOND and THIRD codon position,
        under the sense pool and under Key 1. Establishes that the divisibility
        structure modulo 37 is specific to the third codon position.

  * prime_divisibility_scan.csv
        For every prime p up to 97 and for three value sets -- the
        parametrization-independent groups (no service codons, so 37 cannot be
        imposed), the groups under Key 0 (37 not imposed by construction), and
        the groups under Key 1 (37 imposed by the derivation) -- the number of
        group quantities divisible by p, counted both over all instances and
        over distinct nonzero values. Establishes that, after the trivial
        moduli 2 (proton parity) and 5 (decimal round-number scale), 37 is the
        only modulus dividing a non-chance share of the quantities, and that
        this holds in the value sets where 37 is not imposed, so the result is
        not an artifact of the Key 1 derivation.

Output rules match reproduce.py:
  - Encoding UTF-8, LF line endings, trailing newline.
  - Quoting: csv.QUOTE_NONNUMERIC (string fields quoted, numbers bare).
  - One deterministic sort rule per file. No randomization, no float numerics.
"""

import csv
import os
from math import gcd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

BASES = "ACGT"
ALL_CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
SERVICE_CODONS = {"TAA", "TAG", "TGA", "ATG"}

# The six two-nucleotide subsets, grouped into the three complementary axes.
AXES = [
    ("Keto",       "Amino",      frozenset("GT"), frozenset("AC")),
    ("Strong",     "Weak",       frozenset("CG"), frozenset("AT")),
    ("Purine",     "Pyrimidine", frozenset("AG"), frozenset("CT")),
]
# Flat list (name, members) preserving a deterministic order.
AXIS_GROUPS = []
for n1, n2, s1, s2 in AXES:
    AXIS_GROUPS.append((n1, s1))
    AXIS_GROUPS.append((n2, s2))

PRIME_LIMIT = 97


def read_csv(filename):
    with open(os.path.join(SCRIPT_DIR, filename), newline="",
              encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(filename, header, rows):
    path = os.path.join(SCRIPT_DIR, filename)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC, lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow(r)


# ── Read inputs ────────────────────────────────────────────────────────────

aa_rows = read_csv("amino_acids_nucleons.csv")
codon_rows = read_csv("genetic_code_codons.csv")
group_rows = read_csv("codon_groups.csv")
key_rows = read_csv("key_parameters.csv")

AA_BY_NAME = {
    r["Amino_Acid"]: {"P": int(r["Protons"]), "N": int(r["Neutrons"])}
    for r in aa_rows
}

CODON_TO_AA = {}
for r in codon_rows:
    for c in r["Codons"].split(";"):
        CODON_TO_AA[c.strip()] = r["Product"]

KEY_PARAMS = {}
for r in key_rows:
    KEY_PARAMS.setdefault(r["Key_ID"], {})[r["Codon"]] = {
        "P": int(r["Final_Protons"]), "N": int(r["Final_Neutrons"]),
    }


# ── Nucleon counting ───────────────────────────────────────────────────────

def codon_pn(codon, mode):
    """Return (P, N) for a codon under a counting mode.

    mode == "sense" : service codons are excluded (contribute nothing).
    mode == a key id : service codons take their Key parameters.
    Sense (non-service) codons always take their amino acid's (P, N).
    """
    if codon in SERVICE_CODONS:
        if mode == "sense":
            return (0, 0)
        kp = KEY_PARAMS[mode][codon]
        return (kp["P"], kp["N"])
    aa = AA_BY_NAME[CODON_TO_AA[codon]]
    return (aa["P"], aa["N"])


def group_quantities(codons, mode):
    """Return (T, P, N, Delta) summed over a codon list under a counting mode,
    together with the number of contributing (counted) codons."""
    p = n = counted = 0
    for c in codons:
        if mode == "sense" and c in SERVICE_CODONS:
            continue
        cp, cn = codon_pn(c, mode)
        p += cp
        n += cn
        counted += 1
    return (p + n, p, n, p - n, counted)


# ── Control 1: chemical axes at each codon position ────────────────────────

def build_position_analysis():
    header = ["Position", "Axis", "Members", "Pool", "Group_Size",
              "Total_Nucleons", "Protons", "Neutrons", "Delta_P_N",
              "Total_mod37", "Protons_mod37", "Neutrons_mod37", "Delta_mod37"]
    rows = []
    for pos in (1, 2, 3):
        for pool in ("sense", "key1"):
            for axis_name, members in AXIS_GROUPS:
                codons = [c for c in ALL_CODONS if c[pos - 1] in members]
                T, P, N, D, counted = group_quantities(codons, pool)
                members_str = "{" + ", ".join(sorted(members)) + "}"
                rows.append([
                    pos, axis_name, members_str, pool, counted,
                    T, P, N, D,
                    T % 37, P % 37, N % 37, D % 37,
                ])
    return header, rows


def verify_position3_against_registry():
    # Canonical cross-check: the axis groups this script builds at the THIRD
    # codon position must reproduce the full-code two-nucleotide groups of the
    # published registry (codon_groups.csv). Positions 1 and 2 have no registry
    # counterpart -- the registry defines third-position groups only -- and are
    # therefore built directly from the codon set. Failing loudly here prevents
    # the generated third-position groups from silently diverging from the
    # registry the rest of the analysis is based on.
    registry = {}
    for gr in group_rows:
        name = gr["Group_Name"]
        for axis_name, _ in AXIS_GROUPS:
            if name.startswith(axis_name + ":"):
                registry[axis_name] = frozenset(
                    c.strip() for c in gr["Codon_List"].split(";"))
    for axis_name, members in AXIS_GROUPS:
        built = frozenset(c for c in ALL_CODONS if c[2] in members)
        if built != registry.get(axis_name):
            raise AssertionError(
                f"third-position group '{axis_name}' built by controls.py does "
                f"not match codon_groups.csv")


# ── Control 2: prime divisibility scan over Key 1 group quantities ─────────

def sieve(limit):
    flags = [True] * (limit + 1)
    flags[0] = flags[1] = False
    for i in range(2, int(limit ** 0.5) + 1):
        if flags[i]:
            for j in range(i * i, limit + 1, i):
                flags[j] = False
    return [i for i in range(2, limit + 1) if flags[i]]


def reduced(num, den):
    if den == 0:
        return "0/0"
    g = gcd(num, den) or 1
    return f"{num // g}/{den // g}"


def build_modulus_scan():
    # Recompute the 33 group quantities under both keys from codon_groups.csv,
    # using the same counting rule as reproduce.py. A group is
    # parametrization-independent -- in the literal sense -- when its four
    # nucleon quantities are identical under Key 0 and Key 1; equivalently
    # (Section 3.1) these are exactly Octet I with its sub-groups and the
    # pyrimidine-restricted groups, the only groups containing no service codon.
    # 37 cannot be imposed on this value set by any parametrization.
    key0_values, key1_values, indep_values = [], [], []
    for gr in group_rows:
        codons = [c.strip() for c in gr["Codon_List"].split(";")]
        q0 = group_quantities(codons, "key0")[:4]
        q1 = group_quantities(codons, "key1")[:4]
        key0_values.extend(q0)
        key1_values.extend(q1)
        if q0 == q1:                       # value invariant under the key
            indep_values.extend(q1)

    value_sets = [
        ("param_independent", indep_values),
        ("key0", key0_values),
        ("key1", key1_values),
    ]

    # Counts are pure integers; the chance expectation of divisible distinct
    # values, Total_Distinct / Prime, is stored as an exact reduced fraction
    # (no float numerics). The excess factor over chance is read off as
    # Divisible_Distinct relative to Expected_Distinct.
    header = ["Value_Set", "Prime", "Divisible_All", "Total_All",
              "Divisible_Distinct", "Total_Distinct", "Expected_Distinct"]
    primes = sieve(PRIME_LIMIT)
    rows = []
    for set_name, values in value_sets:
        distinct = sorted({v for v in values if v != 0})
        nd = len(distinct)
        for p in primes:
            d_all = sum(1 for v in values if v != 0 and v % p == 0)
            d_dist = sum(1 for v in distinct if v % p == 0)
            rows.append([set_name, p, d_all, len(values),
                         d_dist, nd, reduced(nd, p)])
    return header, rows


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    verify_position3_against_registry()
    h1, r1 = build_position_analysis()
    write_csv("position_axis_analysis.csv", h1, r1)
    h2, r2 = build_modulus_scan()
    write_csv("prime_divisibility_scan.csv", h2, r2)
    print("Done. 2 control files written.")


if __name__ == "__main__":
    main()
