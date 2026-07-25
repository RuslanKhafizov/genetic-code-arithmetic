#!/usr/bin/env python3
"""
controls.py — Supplementary specificity controls for the standard genetic code.

This script is NOT part of the core reproducible pipeline (reproduce.py). It is
a standalone, deterministic, exact-integer companion that consumes the same
input datasets and writes four supplementary control datasets:

  * position_axis_analysis.csv
        The three chemical axes (Keto/Amino, Strong/Weak, Purine/Pyrimidine)
        applied separately to the FIRST, SECOND and THIRD codon position,
        under the sense pool, Key 0, and Key 1. Establishes that the divisibility
        structure modulo 37 is specific to the third codon position.

  * prime_divisibility_scan.csv
        For every prime p up to 97 and for three value sets -- the
        parametrization-independent groups (no service codons, so 37 cannot be
        imposed), the groups under Key 0 (37 not imposed by construction), and
        the groups under Key 1 (37 imposed by the derivation) -- the number of
        group quantities divisible by p, counted both over all instances and
        over distinct nonzero values. Establishes that, after the small primes
        2 (proton parity) and 5 (the round-number regularity), 37 is the only
        modulus dividing a non-chance share of the quantities, and that this
        holds in the value sets where 37 is not imposed, so the result is not
        an artifact of the Key 1 derivation. The scan ranges over primes, not
        all integers, because divisibility by a composite is the conjunction of
        divisibilities by its prime-power factors (no new modulus), and because
        a higher prime power only re-scores the SAME quantities against a
        smaller chance rate 1/m and so inflates the excess factor without
        adding information. 37 is distinguished from 2 and 5 on base-independent
        grounds: it is a larger prime (rarer per hit under the null), has no
        parity-type structural cause, and occurs at a single power -- no
        quantity is divisible by 37^2, so its excess does not grow with the
        modulus, unlike 5, whose multiples here are almost all multiples of 25.

  * position_asymmetry_ncbi.csv
        For each of the 27 NCBI translation tables, the Keto/Amino
        functional-group axis ({G,T} vs {A,C}) sense-pool asymmetry for T, P, N
        and Delta at each of the three codon positions, with residues modulo 37
        under two parallel sense-pool models. Model 1 excludes every codon
        listed as a stop for that table and the ATG initiator. Model 2 returns
        context-dependent stop codons to the pool with their amino-acid
        assignments. Under Model 1, Tables 1, 11 and 28 reach all three
        positions; under Model 2, only Tables 1 and 11 do. In neither model
        does a table with an analyzed codon-to-amino-acid pool nonequivalent
        to the standard code reach all three positions.

  * representation_sensitivity.csv
        The load-bearing sense-pool quantities recomputed under alternative
        molecular representations of the encoded amino acid: the neutral free
        amino acid (baseline), the peptide residue (free minus one water), and
        the zwitterion, each a fixed per-residue nucleon offset from the free
        amino acid. Establishes which regularities are properties of the
        representation and which are not: every relational quantity (the
        Trp/Ile and the Keto/Amino, Strong/Weak axis differences) is invariant
        across representations, whereas the absolute lattice alignments
        (N(All sense) = 97*37, the deficit, the axis half-pools) hold only for
        the free amino acid, since an offset (dP, dN) preserves neutron
        divisibility by 37 only when dN is itself a multiple of 37.

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

# NCBI "ncbieaa" string codon order: base1, base2, base3 each cycling T,C,A,G.
NCBI_ORDER = "TCAG"
NCBI_CODONS = [NCBI_ORDER[i // 16] + NCBI_ORDER[(i // 4) % 4] + NCBI_ORDER[i % 4]
               for i in range(64)]

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
ncbi_rows = read_csv("ncbi_genetic_code_registry.csv")

AA_BY_NAME = {
    r["Amino_Acid"]: {"P": int(r["Protons"]), "N": int(r["Neutrons"])}
    for r in aa_rows
}

CODON_TO_AA = {}
for r in codon_rows:
    for c in r["Codons"].split(";"):
        CODON_TO_AA[c.strip()] = r["Product"]

# One-letter amino acid code -> (P, N), for decoding NCBI ncbieaa strings.
ONE_LETTER_PN = {}
for r in codon_rows:
    if r["Product"] in AA_BY_NAME:
        ONE_LETTER_PN[r["One_Letter"]] = (AA_BY_NAME[r["Product"]]["P"],
                                          AA_BY_NAME[r["Product"]]["N"])

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
        for pool in ("sense", "key0", "key1"):
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


def odd_core(value):
    """The odd part of a positive integer (all factors of two removed).

    The nested group hierarchy (Octet I and its 16- and 8-codon subgroups)
    mechanically produces, for any quantity V, the multiples 2V and 4V as
    well. A prime that divides V then scores up to three times over the same
    underlying fact. Reducing each divisible value to its odd core collapses
    those power-of-two duplicates. The resulting count records distinct odd
    cores after this structural deduplication; it does not establish
    statistical independence among the remaining quantities."""
    while value % 2 == 0:
        value //= 2
    return value


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
              "Divisible_Distinct", "Total_Distinct", "Expected_Distinct",
              "Divisible_Distinct_OddCores"]
    # The scan ranges over primes only: divisibility by a composite is the
    # conjunction of its prime-power factors (no new modulus), and a higher
    # prime power merely re-scores the same values against a smaller chance
    # rate 1/m, inflating the excess factor without adding information. The
    # prime is therefore the irreducible unit of comparison.
    primes = sieve(PRIME_LIMIT)
    rows = []
    for set_name, values in value_sets:
        distinct = sorted({v for v in values if v != 0})
        nd = len(distinct)
        for p in primes:
            d_all = sum(1 for v in values if v != 0 and v % p == 0)
            divisible = [v for v in distinct if v % p == 0]
            d_dist = len(divisible)
            d_cores = len({odd_core(v) for v in divisible})
            rows.append([set_name, p, d_all, len(values),
                         d_dist, nd, reduced(nd, p), d_cores])
    return header, rows


# ── Main ───────────────────────────────────────────────────────────────────

def build_position_asymmetry_ncbi():
    # Control 3. For each of the 27 NCBI translation tables, the Keto/Amino
    # functional-group axis ({G,T} vs {A,C}) sense-pool asymmetry
    # (Keto minus Amino) for T, P, N and Delta at each of the three codon
    # positions, with residues modulo 37, under two parallel pool models.
    #
    # Model 1 excludes ATG and every codon listed in the table's Stop_Codons
    # field, including context-dependent stops. Model 2 excludes ATG and only
    # positions marked '*' in Amino_Acids, thereby returning context-dependent
    # stops to the pool with their amino-acid assignments. The two models differ
    # only for Tables 27, 28 and 31 in the fixed registry snapshot.
    header = ["Transl_Table", "Code_Name", "Pool_Model", "Position",
              "Keto_minus_Amino_T", "Keto_minus_Amino_P",
              "Keto_minus_Amino_N", "Keto_minus_Amino_Delta",
              "T_mod37", "P_mod37", "N_mod37", "Delta_mod37",
              "Position_Has_Div37", "Positions_With_Div37"]
    rows = []
    full_signature = {"Model_1": set(), "Model_2": set()}
    for tr in ncbi_rows:
        table = int(tr["Transl_Table"])
        name = tr["Code_Name"]
        aa_string = tr["Amino_Acids"]
        listed_stops = {
            codon.strip()
            for codon in tr["Stop_Codons"].split(",")
            if codon.strip()
        }
        table_results = {}
        for model in ("Model_1", "Model_2"):
            sense = {}
            for i, codon in enumerate(NCBI_CODONS):
                aa = aa_string[i]
                if model == "Model_1":
                    is_stop = codon in listed_stops
                else:
                    is_stop = aa == "*"
                if is_stop or codon == "ATG":
                    continue
                sense[codon] = ONE_LETTER_PN[aa]

            per_position = []
            for pos in (1, 2, 3):
                pk = nk = pa = na = 0
                for codon, (p, n) in sense.items():
                    base = codon[pos - 1]
                    if base in "GT":        # Keto
                        pk += p
                        nk += n
                    elif base in "AC":      # Amino
                        pa += p
                        na += n
                T = (pk + nk) - (pa + na)
                P = pk - pa
                N = nk - na
                D = (pk - nk) - (pa - na)
                has = int(any(v % 37 == 0 for v in (T, P, N, D)))
                per_position.append((pos, T, P, N, D, has))

            table_results[model] = per_position
            positions_with = sum(h for *_, h in per_position)
            if positions_with == 3:
                full_signature[model].add(table)
            for pos, T, P, N, D, has in per_position:
                rows.append([table, name, model, pos, T, P, N, D,
                             T % 37, P % 37, N % 37, D % 37,
                             has, positions_with])

        if (table not in {27, 28, 31}
                and table_results["Model_1"] != table_results["Model_2"]):
            raise AssertionError(
                f"position asymmetry: Models 1 and 2 differ for Table {table}"
            )

    expected = {
        "Model_1": {1, 11, 28},
        "Model_2": {1, 11},
    }
    if full_signature != expected:
        raise AssertionError(
            "position asymmetry: unexpected full-signature tables: "
            f"{full_signature!r}"
        )

    rows.sort(key=lambda r: (r[0], r[2], r[3]))
    return header, rows


# -- Control 4: representation (counting-convention) sensitivity ------------

# Alternative molecular representations of the amino acid a codon encodes, each
# a fixed per-residue nucleon offset from the neutral free amino acid (the
# baseline used throughout). The free amino acid is the canonical molecular
# identity of the encoded species; a peptide residue is that molecule minus one
# water lost on peptide-bond formation (H2O = 10 protons, 8 neutrons); the
# zwitterion has the same atoms as the neutral form. All are UNIFORM offsets, so
# every difference between codon sets of equal cardinality is invariant across
# them, while an absolute sum's divisibility by 37 survives an offset (dP, dN)
# only when 60*dN is a multiple of 37, i.e. dN is a multiple of 37 (gcd(60,37)=1)
# -- which among chemically meaningful representations selects the free amino
# acid alone.
REPRESENTATIONS = [
    ("free",            0,   0),   # neutral free amino acid (baseline)
    ("peptide_residue", -10, -8),  # free minus one H2O (10 protons, 8 neutrons)
    ("zwitterion",      0,   0),   # same atoms as the neutral free form
]

def _octet1_codons():
    """Octet I: the eight four-fold-degenerate boxes (32 codons), built from the
    degeneracy rule rather than a group name, matching reproduce.py's Octet I."""
    box = {}
    for c in ALL_CODONS:
        box.setdefault(c[:2], []).append(c)
    return [c for _, cs in box.items()
            if len({CODON_TO_AA[x] for x in cs}) == 1
            and CODON_TO_AA[cs[0]] != "Stop"
            for c in cs]

def _shifted_sense(codons, dP, dN):
    """(T, P, N, Delta) over the SENSE codons of a list, each amino acid shifted
    by a per-residue offset (dP, dN); service codons excluded."""
    p = n = 0
    for c in codons:
        if c in SERVICE_CODONS:
            continue
        aa = AA_BY_NAME[CODON_TO_AA[c]]
        p += aa["P"] + dP
        n += aa["N"] + dN
    return (p + n, p, n, p - n)

def _representation_anchors(dP, dN):
    """Every anchor quantity under a representation offset, as (name, class, value)."""
    sense = [c for c in ALL_CODONS if c not in SERVICE_CODONS]
    third = lambda members: [c for c in sense if c[2] in members]
    keto, amino = third(frozenset("GT")), third(frozenset("AC"))
    strong, weak = third(frozenset("CG")), third(frozenset("AT"))
    Tall, Pall, Nall, _ = _shifted_sense(sense, dP, dN)
    TocI = _shifted_sense(_octet1_codons(), dP, dN)[0]
    _, Pk, Nk, _ = _shifted_sense(keto,   dP, dN)
    _, Pa, Na, _ = _shifted_sense(amino,  dP, dN)
    _, Ps, Ns, _ = _shifted_sense(strong, dP, dN)
    _, Pw, Nw, _ = _shifted_sense(weak,   dP, dN)
    trp, ile = AA_BY_NAME[CODON_TO_AA["TGG"]], AA_BY_NAME[CODON_TO_AA["ATT"]]
    return [
        ("N_sense",                     "absolute",   Nall),
        ("P_sense",                     "absolute",   Pall),
        ("deficit_T_OctetI_minus_N_sense", "absolute", TocI - Nall),
        ("N_Keto",                      "absolute",   Nk),
        ("N_Amino",                     "absolute",   Na),
        ("N_Strong",                    "absolute",   Ns),
        ("N_Weak",                      "absolute",   Nw),
        ("dN_Trp_Ile",                  "relational", (trp["N"] + dN) - (ile["N"] + dN)),
        ("dP_Trp_Ile",                  "relational", (trp["P"] + dP) - (ile["P"] + dP)),
        ("KetoAmino_dN",                "relational", Nk - Na),
        ("KetoAmino_dP",                "relational", Pk - Pa),
        ("StrongWeak_dN",               "relational", Ns - Nw),
        ("StrongWeak_dP",               "relational", Ps - Pw),
    ]

def build_representation_sensitivity():
    header = ["Representation", "Delta_P", "Delta_N", "Anchor", "Class",
              "Value", "Value_mod37", "Matches_Free"]
    free_val = {name: val for name, _, val in _representation_anchors(0, 0)}
    rows = []
    for rep_name, dP, dN in REPRESENTATIONS:
        for name, cls, val in _representation_anchors(dP, dN):
            rows.append([rep_name, dP, dN, name, cls,
                         val, val % 37, val == free_val[name]])
    return header, rows

def verify_representation_free_against_preprint():
    # The free-amino-acid row must reproduce the preprint's sense-pool values.
    want = {"N_sense": 3589, "P_sense": 4180,
            "deficit_T_OctetI_minus_N_sense": 111,
            "N_Keto": 1813, "N_Amino": 1776,
            "dN_Trp_Ile": 37, "dP_Trp_Ile": 36}
    got = {name: val for name, _, val in _representation_anchors(0, 0)}
    for k, v in want.items():
        if got.get(k) != v:
            raise AssertionError(
                f"representation control: free '{k}' = {got.get(k)}, "
                f"expected preprint value {v}")


def main():
    verify_position3_against_registry()
    h1, r1 = build_position_analysis()
    write_csv("position_axis_analysis.csv", h1, r1)
    h2, r2 = build_modulus_scan()
    write_csv("prime_divisibility_scan.csv", h2, r2)
    h3, r3 = build_position_asymmetry_ncbi()
    write_csv("position_asymmetry_ncbi.csv", h3, r3)
    verify_representation_free_against_preprint()
    h4, r4 = build_representation_sensitivity()
    write_csv("representation_sensitivity.csv", h4, r4)
    print("Done. 4 control files written.")


if __name__ == "__main__":
    main()
