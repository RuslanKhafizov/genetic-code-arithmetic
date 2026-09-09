#!/usr/bin/env python3
"""
reproduce.py — Generates all 18 derived CSV files from 5 input CSV files.

Unified output rules for every file:
  - Encoding UTF-8, LF line endings, trailing newline.
  - Quoting: csv.QUOTE_NONNUMERIC (string fields quoted, numbers bare).
  - One deterministic sort rule per file, applied to every row.

Group-level outputs are generated for pool60, Key 0, and Key 1. The
pool60 dataset excludes ATG, TAA, TAG, and TGA from every group; it does
not model them as zero-valued codons.
"""

import csv
import os
from itertools import combinations

ACTIVE_KEYS = ["key0", "key1"]
POOL60_ID = "pool60"
OUTPUT_DATASETS = [POOL60_ID] + ACTIVE_KEYS
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


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


def gcd(a, b):
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def reduced_fraction(p, q):
    if q == 0:
        return "0/0"
    g = gcd(p, q)
    return f"{p // g}/{q // g}"


# ── Read inputs ────────────────────────────────────────────────────────────

aa_rows = read_csv("amino_acids_nucleons.csv")
codon_rows = read_csv("genetic_code_codons.csv")
group_rows = read_csv("codon_groups.csv")
key_rows = read_csv("key_parameters.csv")
ncbi_rows = read_csv("ncbi_genetic_code_registry.csv")

AA_BY_NAME, AA_BY_CODE1, AA_BY_THREE, AA_ORDER = {}, {}, {}, []
for idx, r in enumerate(aa_rows):
    name = r["Amino_Acid"]
    AA_BY_NAME[name] = {
        "P": int(r["Protons"]), "N": int(r["Neutrons"]),
        "T": int(r["Total_Nucleons"]), "D": int(r["Delta_P_N"]),
        "idx": idx,
    }
    AA_ORDER.append(name)

for r in codon_rows:
    if r["Product"] == "STOP":
        continue
    d = AA_BY_NAME[r["Product"]]
    d["three"] = r["Three_Letter"]
    AA_BY_THREE[r["Three_Letter"]] = d
    AA_BY_CODE1[r["One_Letter"]] = d

CODON_TO_AA = {}
for r in codon_rows:
    for c in r["Codons"].split(";"):
        codon = c.strip()
        if codon in CODON_TO_AA:
            raise ValueError("Codon {} is listed more than once".format(codon))
        CODON_TO_AA[codon] = r["Product"]

KEY_PARAMS = {}
for r in key_rows:
    KEY_PARAMS.setdefault(r["Key_ID"], {})[r["Codon"]] = {
        "P": int(r["Final_Protons"]), "N": int(r["Final_Neutrons"]),
    }

SERVICE_CODONS = {"TAA", "TAG", "TGA", "ATG"}


def group_codons(group_row):
    return [c.strip() for c in group_row["Codon_List"].split(";")]


def validate_inputs():
    expected_codons = {
        a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"
    }
    actual_codons = set(CODON_TO_AA)
    if actual_codons != expected_codons:
        missing = sorted(expected_codons - actual_codons)
        extra = sorted(actual_codons - expected_codons)
        raise ValueError(
            "Standard code must contain all 64 codons; missing={}, extra={}".format(
                missing, extra
            )
        )

    if len(group_rows) != 33:
        raise ValueError("codon_groups.csv must contain exactly 33 groups")

    seen_group_ids = set()
    for gr in group_rows:
        group_id = int(gr["Group_ID"])
        if group_id in seen_group_ids:
            raise ValueError("Duplicate Group_ID {}".format(group_id))
        seen_group_ids.add(group_id)

        codons = group_codons(gr)
        if len(codons) != len(set(codons)):
            raise ValueError("Group {} contains duplicate codons".format(group_id))
        if len(codons) != int(gr["Group_Size"]):
            raise ValueError(
                "Group {} size does not match Codon_List".format(group_id)
            )
        unknown = sorted(set(codons) - actual_codons)
        if unknown:
            raise ValueError(
                "Group {} contains unknown codons: {}".format(group_id, unknown)
            )

    expected_group_ids = set(range(1, 34))
    if seen_group_ids != expected_group_ids:
        raise ValueError("Group_ID values must be the integers 1 through 33")

    for key_id in ACTIVE_KEYS:
        missing = sorted(SERVICE_CODONS - set(KEY_PARAMS.get(key_id, {})))
        if missing:
            raise ValueError(
                "{} is missing service-codon parameters: {}".format(
                    key_id, missing
                )
            )

    pool60_codons = actual_codons - SERVICE_CODONS
    if len(pool60_codons) != 60:
        raise ValueError(
            "The service-codon-excluded pool must contain exactly 60 codons"
        )


validate_inputs()


# ── 1. Pairwise differences ────────────────────────────────────────────────

def write_differences(value_key, diff_col, filename):
    buckets = {}
    for a_name, b_name in combinations(AA_ORDER, 2):
        a, b = AA_BY_NAME[a_name], AA_BY_NAME[b_name]
        d = abs(a[value_key] - b[value_key])
        t1, t2 = sorted([a["three"], b["three"]])
        buckets.setdefault(d, []).append(f"{t1}-{t2}")
    max_d = max(buckets)
    rows = []
    for d in range(max_d + 1):
        pairs = sorted(buckets.get(d, []))
        if pairs:
            rows.append([d, len(pairs), "; ".join(pairs)])
        else:
            rows.append([d, 0, "No pairs found"])
    write_csv(filename, [diff_col, "Number_Of_Pairs", "Amino_Acid_Pairs"],
              rows)


# ── 2. Nucleon data per dataset ────────────────────────────────────────────

def compute_nucleon_data(dataset_id):
    exclude_service = dataset_id == POOL60_ID
    kp = None if exclude_service else KEY_PARAMS[dataset_id]
    results = []
    for gr in group_rows:
        tp = tn = 0
        codons = group_codons(gr)
        if exclude_service:
            codons = [c for c in codons if c not in SERVICE_CODONS]
        for codon in codons:
            if codon in SERVICE_CODONS:
                tp += kp[codon]["P"]
                tn += kp[codon]["N"]
            else:
                aa = AA_BY_NAME[CODON_TO_AA[codon]]
                tp += aa["P"]
                tn += aa["N"]
        results.append({
            "Group_Name": gr["Group_Name"],
            "Group_ID": int(gr["Group_ID"]),
            "Twin_Group_ID": gr["Twin_Group_ID"],
            "Code_Section": gr["Code_Section"],
            "Group_Size": len(codons),
            "Total_Nucleons": tp + tn, "Protons": tp, "Neutrons": tn,
            "Delta_P_N": tp - tn,
        })
    results.sort(key=lambda d: d["Group_ID"])
    return results


def twin_field(d):
    t = d["Twin_Group_ID"]
    return int(t) if t else ""


def write_nucleon_data(dataset_id, data):
    header = ["Group_Name", "Group_ID", "Twin_Group_ID", "Code_Section",
              "Group_Size", "Total_Nucleons", "Protons", "Neutrons",
              "Delta_P_N"]
    rows = [[d["Group_Name"], d["Group_ID"], twin_field(d),
             d["Code_Section"], d["Group_Size"], d["Total_Nucleons"],
             d["Protons"], d["Neutrons"], d["Delta_P_N"]] for d in data]
    write_csv(f"{dataset_id}_nucleon-data.csv", header, rows)


# ── 3. Divisibility by 37 ──────────────────────────────────────────────────

def div37_expr(v):
    v = int(v)
    if v < 37:
        return str(v)
    q, r = divmod(v, 37)
    return f"{q}*37" if r == 0 else f"{q}*37+{r}"


def write_divisibility(dataset_id, data):
    header = ["Group_Name", "Group_ID", "Twin_Group_ID", "Code_Section",
              "Group_Size", "Total_mod37", "Protons_mod37",
              "Neutrons_mod37", "Delta_mod37"]
    rows = [[d["Group_Name"], d["Group_ID"], twin_field(d),
             d["Code_Section"], d["Group_Size"],
             div37_expr(d["Total_Nucleons"]), div37_expr(d["Protons"]),
             div37_expr(d["Neutrons"]), div37_expr(d["Delta_P_N"])]
            for d in data]
    write_csv(f"{dataset_id}_divisibility-37.csv", header, rows)


# ── 4. Equalities ──────────────────────────────────────────────────────────

EQ_PARAMS = ["Total_Nucleons", "Protons", "Neutrons", "Delta_P_N"]
MATCH_ORDER = {"DIFF_PARAM": 0, "CROSS_PARAM": 1, "SAME_PARAM": 2}


def write_equalities(dataset_id, data):
    header = ["Nucleon_Value", "Match_Type", "Source_Section", "Source_ID",
              "Source_Name", "Source_Parameter", "Target_Parameter",
              "Target_Name", "Target_ID", "Target_Section"]
    n = len(data)
    eq_rows = []
    for i in range(n):
        for j in range(i, n):
            src, tgt = data[i], data[j]
            same_size = src["Group_Size"] == tgt["Group_Size"]
            for pi, sp in enumerate(EQ_PARAMS):
                for pj, tp in enumerate(EQ_PARAMS):
                    if i == j and pi >= pj:
                        continue
                    if src[sp] != tgt[tp]:
                        continue
                    if not same_size:
                        mt = "DIFF_PARAM"
                    elif sp != tp:
                        mt = "CROSS_PARAM"
                    else:
                        mt = "SAME_PARAM"
                    eq_rows.append([
                        src[sp], mt, src["Code_Section"], src["Group_ID"],
                        src["Group_Name"], sp, tp, tgt["Group_Name"],
                        tgt["Group_ID"], tgt["Code_Section"],
                    ])
    eq_rows.sort(key=lambda r: (r[0], MATCH_ORDER[r[1]], r[3], r[8],
                                EQ_PARAMS.index(r[5]),
                                EQ_PARAMS.index(r[6])))
    write_csv(f"{dataset_id}_equalities.csv", header, eq_rows)
    return len(eq_rows)


# ── 5. Ratios ──────────────────────────────────────────────────────────────

def write_ratios(dataset_id, data):
    header = ["Group_Name", "Group_ID", "Twin_Group_ID", "Code_Section",
              "Group_Size",
              "P_T_Rational", "N_T_Rational", "Delta_T_Rational",
              "P_N_Rational", "Delta_P_Rational", "Delta_N_Rational",
              "P_T_Ratio", "N_T_Ratio", "Delta_T_Ratio",
              "P_N_Ratio", "Delta_P_Ratio", "Delta_N_Ratio"]
    rows = []
    for d in data:
        p, n, t, dl = (d["Protons"], d["Neutrons"],
                       d["Total_Nucleons"], d["Delta_P_N"])
        pairs = [(p, t), (n, t), (dl, t), (p, n), (dl, p), (dl, n)]
        rats = [reduced_fraction(a, b) for a, b in pairs]
        decs = [f"{(a / b):.4f}" if b != 0 else "0.0000"
                for a, b in pairs]
        rows.append([d["Group_Name"], d["Group_ID"], twin_field(d),
                     d["Code_Section"], d["Group_Size"]] + rats + decs)
    write_csv(f"{dataset_id}_ratios.csv", header, rows)


# ── 6. Deficit models analysis ─────────────────────────────────────────────

NCBI_CODONS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
OCTET1_STRUCT = {"GC", "CG", "GG", "CC", "AC", "GT", "CT", "TC"}


def _per_table_metrics(nr):
    """Compute the per-codon metrics for a single NCBI table row.

    Returns a dict with everything needed by both
    write_deficit_models and write_keto_amino_balance.
    All sums are over the 64 codons; service codons (those in the
    excl set) contribute 0 to the corresponding "_pure" tally.
    Stop codons (aa_str letter == '*') contribute 0 to N_all/P_all
    as well.
    """
    aa_str = nr["Amino_Acids"]
    stops = set(c.strip() for c in nr["Stop_Codons"].split(","))
    # Unconditional stops are those the registry marks with '*' in the
    # amino-acid string; every other declared stop is context-dependent.
    # Model 2 returns the context-dependent ones to the pool and keeps the
    # unconditional ones as service positions. Deriving the set this way
    # rather than from the per-table Stop_Context_Dependent flag also
    # covers a table that would declare both kinds at once; on the present
    # registry the two definitions coincide, so published outputs are
    # unchanged.
    unconditional = {c for k, c in enumerate(NCBI_CODONS) if aa_str[k] == "*"}
    assert unconditional <= stops, (nr["Transl_Table"], sorted(unconditional))
    excl_m1 = stops | {"ATG"}
    excl_m2 = unconditional | {"ATG"}

    n_all = p_all = 0
    n_pure_m1 = p_pure_m1 = 0
    n_pure_m2 = p_pure_m2 = 0
    n_keto_m1 = n_amino_m1 = p_keto_m1 = p_amino_m1 = 0
    n_keto_m2 = n_amino_m2 = p_keto_m2 = p_amino_m2 = 0

    for k, codon in enumerate(NCBI_CODONS):
        L = aa_str[k]
        if L == "*":
            n = p = 0
        else:
            aa = AA_BY_CODE1[L]
            n, p = aa["N"], aa["P"]
        n_all += n
        p_all += p
        in_m1 = codon not in excl_m1
        in_m2 = codon not in excl_m2
        if in_m1:
            n_pure_m1 += n
            p_pure_m1 += p
        if in_m2:
            n_pure_m2 += n
            p_pure_m2 += p
        # Keto = third position is G or T; Amino = third is A or C.
        third = codon[2]
        if third in "GT":
            if in_m1:
                n_keto_m1 += n
                p_keto_m1 += p
            if in_m2:
                n_keto_m2 += n
                p_keto_m2 += p
        else:  # third in "AC"
            if in_m1:
                n_amino_m1 += n
                p_amino_m1 += p
            if in_m2:
                n_amino_m2 += n
                p_amino_m2 += p

    # Octet I: 8 structural families (codon[:2] in OCTET1_STRUCT),
    # T = sum(P+N).
    oct_struct_t = 0
    for k, codon in enumerate(NCBI_CODONS):
        if codon[:2] in OCTET1_STRUCT:
            L = aa_str[k]
            if L != "*":
                aa = AA_BY_CODE1[L]
                oct_struct_t += aa["P"] + aa["N"]

    # Functional Octet I: every 4-fold-degenerate family in this table.
    families = {}
    for k, codon in enumerate(NCBI_CODONS):
        families.setdefault(codon[:2], []).append(aa_str[k])
    oct_func_t = 0
    for letters in families.values():
        if "*" not in letters and len(set(letters)) == 1:
            aa = AA_BY_CODE1[letters[0]]
            oct_func_t += 4 * (aa["P"] + aa["N"])

    # Service stop positions carrying the shared stop parameter, per model.
    # Model 2 keeps only unconditional stops; ATG is never a stop.
    service_m1 = set(stops)
    service_m2 = set(unconditional)

    return {
        "tt": int(nr["Transl_Table"]),
        "name": nr["Code_Name"],
        "excl_m1": excl_m1, "excl_m2": excl_m2,
        "service_m1": service_m1, "service_m2": service_m2,
        "n_all": n_all, "p_all": p_all,
        "n_pure_m1": n_pure_m1, "n_pure_m2": n_pure_m2,
        "p_pure_m1": p_pure_m1, "p_pure_m2": p_pure_m2,
        "oct_struct_t": oct_struct_t, "oct_func_t": oct_func_t,
        "n_keto_m1": n_keto_m1, "n_amino_m1": n_amino_m1,
        "p_keto_m1": p_keto_m1, "p_amino_m1": p_amino_m1,
        "n_keto_m2": n_keto_m2, "n_amino_m2": n_amino_m2,
        "p_keto_m2": p_keto_m2, "p_amino_m2": p_amino_m2,
    }


def write_deficit_models(metrics_list):
    header = [
        "Transl_Table", "Code_Name",
        "Pure_Codons_Count", "Pure_Codons_Count_Alt",
        "N_All", "N_Pure", "N_Pure_Alt",
        "Octet1_Struct_T", "Octet1_Func_T",
        "Deficit_vs_3700", "Deficit_vs_Struct", "Deficit_vs_Func",
        "Deficit_vs_3700_Alt", "Deficit_vs_Struct_Alt",
        "Deficit_vs_Func_Alt",
        # Test 1: proton residues mod 37 of the sense pool.
        "P_All", "P_Pure", "P_Pure_Alt",
        "P_Pure_mod37", "P_Pure_Alt_mod37",
    ]
    rows = []
    for m in metrics_list:
        rows.append([
            m["tt"], m["name"],
            64 - len(m["excl_m1"]), 64 - len(m["excl_m2"]),
            m["n_all"], m["n_pure_m1"], m["n_pure_m2"],
            m["oct_struct_t"], m["oct_func_t"],
            3700 - m["n_pure_m1"],
            m["oct_struct_t"] - m["n_pure_m1"],
            m["oct_func_t"] - m["n_pure_m1"],
            3700 - m["n_pure_m2"],
            m["oct_struct_t"] - m["n_pure_m2"],
            m["oct_func_t"] - m["n_pure_m2"],
            m["p_all"], m["p_pure_m1"], m["p_pure_m2"],
            m["p_pure_m1"] % 37, m["p_pure_m2"] % 37,
        ])
    rows.sort(key=lambda r: r[0])
    write_csv("deficit_models_analysis.csv", header, rows)


def write_keto_amino_balance(metrics_list):
    """Test 2: Keto vs Amino balance of the sense pool.

    Keto = codons with G or T in the third position;
    Amino = codons with A or C in the third position.
    Computed only over the sense pool (excl set removed) under each
    of the two models.
    """
    header = [
        "Transl_Table", "Code_Name",
        "N_Keto_Pure", "N_Amino_Pure", "N_Diff_Pure",
        "P_Keto_Pure", "P_Amino_Pure", "P_Diff_Pure",
        "N_Keto_Pure_Alt", "N_Amino_Pure_Alt", "N_Diff_Pure_Alt",
        "P_Keto_Pure_Alt", "P_Amino_Pure_Alt", "P_Diff_Pure_Alt",
    ]
    rows = []
    for m in metrics_list:
        rows.append([
            m["tt"], m["name"],
            m["n_keto_m1"], m["n_amino_m1"],
            m["n_keto_m1"] - m["n_amino_m1"],
            m["p_keto_m1"], m["p_amino_m1"],
            m["p_keto_m1"] - m["p_amino_m1"],
            m["n_keto_m2"], m["n_amino_m2"],
            m["n_keto_m2"] - m["n_amino_m2"],
            m["p_keto_m2"], m["p_amino_m2"],
            m["p_keto_m2"] - m["p_amino_m2"],
        ])
    rows.sort(key=lambda r: r[0])
    write_csv("keto_amino_balance_models.csv", header, rows)


# ---------------------------------------------------------------------------
# Minimal Keto/Amino compensation, with and without global divisibility by 37.
#
# For each quantity Q in {P, N} independently:
#   Q_pool   -- service-excluded pool total for the model in question;
#   dQ       -- Q_pool(Keto) - Q_pool(Amino);
#   s_k, s_a -- stop positions in the Keto and Amino branches;
#   S        -- s_k + s_a;
#   y        -- the shared Q value carried by every stop codon;
#   x        -- the separate Q value carried by ATG (ATG lies in Keto).
#
# Balance:    dQ + s_k*y + x - s_a*y = 0   =>   x = (s_a - s_k)*y - dQ
# Objective:  C_Q = S*y + x = 2*s_a*y - dQ   (total added contribution)
#
# Two procedures:
#   Bal -- balance, non-negativity, integrality, minimal C_Q;
#   Mod -- the same plus (Q_pool + C_Q) % 37 == 0.
#
# P and N are solved independently. No equality between the P and N stop
# parameters is imposed, and no unit offset on ATG is assumed.
#
# Because C_Q = 2*s_a*y - dQ is strictly increasing in y whenever s_a > 0,
# minimising the total contribution selects the smallest admissible y even
# when s_a <= s_k. When s_a == 0 the contribution does not depend on y at
# all: the total is fixed while the split between stops and ATG is not,
# which is reported as "total_only".
# ---------------------------------------------------------------------------

BALANCE_MODULUS = 37


def _egcd(a, b):
    """Extended Euclid: return (g, u, v) with a*u + b*v == g == gcd(a, b)."""
    old_r, r = a, b
    old_u, u = 1, 0
    old_v, v = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_u, u = u, old_u - q * u
        old_v, v = v, old_v - q * v
    return old_r, old_u, old_v


def _admissible_y(dq, s_k, s_a):
    """Integer bounds (lo, hi) for y >= 0 with x = (s_a - s_k)*y - dq >= 0.

    hi is None when y is unbounded above. Returns None when no y qualifies.
    """
    d = s_a - s_k
    if d > 0:
        lo = -((-dq) // d)          # ceil(dq / d)
        return (max(0, lo), None)
    if d == 0:
        return (0, None) if dq <= 0 else None
    hi = dq // d                    # floor(dq / d); d < 0 flips the inequality
    return (0, hi) if hi >= 0 else None


def _solve_component(dq, q_pool, s_k, s_a, require_divisible):
    """Minimal compensation for one quantity.

    Returns a dict with:
      status  -- "unique", "total_only" or "none";
      y, x    -- the parameters when they are uniquely determined, else None
                 (y is always None when there is no stop position at all);
      added   -- C_Q, defined whenever status != "none";
      total   -- q_pool + C_Q, defined whenever status != "none".
    """
    empty = {"status": "none", "y": None, "x": None,
             "added": None, "total": None}
    stops = s_k + s_a

    if stops == 0:
        # No stop position exists: y is absent and x is forced by the balance.
        x = -dq
        if x < 0:
            return empty
        added = x
        if require_divisible and (q_pool + added) % BALANCE_MODULUS != 0:
            return empty
        return {"status": "unique", "y": None, "x": x,
                "added": added, "total": q_pool + added}

    span = _admissible_y(dq, s_k, s_a)
    if span is None:
        return empty
    lo, hi = span

    if s_a == 0:
        # C_Q does not depend on y: the total is fixed, the split is not.
        added = -dq
        if require_divisible and (q_pool + added) % BALANCE_MODULUS != 0:
            return empty
        count = hi - lo + 1          # hi is finite here because s_k > 0
        if count == 1:
            y = lo
            return {"status": "unique", "y": y, "x": (s_a - s_k) * y - dq,
                    "added": added, "total": q_pool + added}
        return {"status": "total_only", "y": None, "x": None,
                "added": added, "total": q_pool + added}

    # s_a > 0: the objective is strictly increasing in y.
    if not require_divisible:
        y = lo
    else:
        a = (2 * s_a) % BALANCE_MODULUS
        b = (dq - q_pool) % BALANCE_MODULUS
        g, u, _ = _egcd(a, BALANCE_MODULUS)
        if b % g != 0:
            return empty
        step = BALANCE_MODULUS // g
        y0 = (u * (b // g)) % step
        # smallest y >= lo congruent to y0 modulo step
        y = y0 + ((lo - y0 + step - 1) // step) * step
    if hi is not None and y > hi:
        return empty
    x = (s_a - s_k) * y - dq
    added = 2 * s_a * y - dq
    return {"status": "unique", "y": y, "x": x,
            "added": added, "total": q_pool + added}


def _fmt(v):
    """CSV cell: empty string for an undefined value."""
    return "" if v is None else v


def _flag(v):
    return "" if v is None else ("True" if v else "False")


def write_service_codon_balance_comparison(metrics_list):
    """Compare minimal Keto/Amino compensation with and without divisibility.

    One row per (translation table, sense-pool model). Model 1 excludes ATG
    and every declared stop position, context-dependent ones included;
    Model 2 returns context-dependent positions to the pool and keeps only
    unconditional stops as service positions. In the present registry all
    three context-dependent tables (27, 28, 31) declare no unconditional
    stop, so their Model 2 rows carry no stop parameter at all.
    """
    header = ["Transl_Table", "Code_Name", "Model",
              "Stop_Count", "Stops_Keto", "Stops_Amino", "Octet1_Struct_T"]
    for pref in ("Bal", "Mod"):
        header += [f"{pref}_Stop_P", f"{pref}_Stop_N",
                   f"{pref}_Start_P", f"{pref}_Start_N",
                   f"{pref}_Added_P", f"{pref}_Added_N",
                   f"{pref}_P_All", f"{pref}_N_All",
                   f"{pref}_Status_P", f"{pref}_Status_N",
                   f"{pref}_N_Diff_Octet1_Struct",
                   f"{pref}_N_All_eq_Octet1_Struct"]
    header += ["Bal_P_Divisible37", "Bal_N_Divisible37",
               "Cost_P", "Cost_N", "Mod_Stop_Params_Equal"]

    rows = []
    for m in metrics_list:
        for model in (1, 2):
            sfx = "m1" if model == 1 else "m2"
            service = sorted(m[f"service_{sfx}"])
            s_k = sum(1 for c in service if c[2] in "GT")
            s_a = len(service) - s_k
            p_pool = m[f"p_pure_{sfx}"]
            n_pool = m[f"n_pure_{sfx}"]
            dp = m[f"p_keto_{sfx}"] - m[f"p_amino_{sfx}"]
            dn = m[f"n_keto_{sfx}"] - m[f"n_amino_{sfx}"]
            oct1 = m["oct_struct_t"]

            sol = {}
            for pref, need in (("Bal", False), ("Mod", True)):
                sol[pref] = {
                    "P": _solve_component(dp, p_pool, s_k, s_a, need),
                    "N": _solve_component(dn, n_pool, s_k, s_a, need),
                }

            row = [m["tt"], m["name"], model,
                   len(service), s_k, s_a, oct1]
            for pref in ("Bal", "Mod"):
                sp, sn = sol[pref]["P"], sol[pref]["N"]
                n_diff = None if sn["total"] is None else sn["total"] - oct1
                row += [_fmt(sp["y"]), _fmt(sn["y"]),
                        _fmt(sp["x"]), _fmt(sn["x"]),
                        _fmt(sp["added"]), _fmt(sn["added"]),
                        _fmt(sp["total"]), _fmt(sn["total"]),
                        sp["status"], sn["status"],
                        _fmt(n_diff),
                        _flag(None if n_diff is None else n_diff == 0)]

            bal_p, bal_n = sol["Bal"]["P"], sol["Bal"]["N"]
            mod_p, mod_n = sol["Mod"]["P"], sol["Mod"]["N"]
            div_p = (None if bal_p["total"] is None
                     else bal_p["total"] % BALANCE_MODULUS == 0)
            div_n = (None if bal_n["total"] is None
                     else bal_n["total"] % BALANCE_MODULUS == 0)
            cost_p = (None if (bal_p["added"] is None or mod_p["added"] is None)
                      else mod_p["added"] - bal_p["added"])
            cost_n = (None if (bal_n["added"] is None or mod_n["added"] is None)
                      else mod_n["added"] - bal_n["added"])
            equal = None
            if service and mod_p["y"] is not None and mod_n["y"] is not None:
                equal = mod_p["y"] == mod_n["y"]
            row += [_flag(div_p), _flag(div_n),
                    _fmt(cost_p), _fmt(cost_n), _flag(equal)]
            rows.append(row)

    rows.sort(key=lambda r: (r[0], r[2]))
    write_csv("service_codon_balance_comparison.csv", header, rows)


def main():
    write_differences("P", "Proton_Difference",
                      "amino_acids_proton_differences.csv")
    write_differences("N", "Neutron_Difference",
                      "amino_acids_neutron_differences.csv")
    write_differences("T", "Nucleon_Difference",
                      "amino_acids_nucleon_differences.csv")
    counts = {}
    for dataset_id in OUTPUT_DATASETS:
        data = compute_nucleon_data(dataset_id)
        write_nucleon_data(dataset_id, data)
        write_divisibility(dataset_id, data)
        counts[dataset_id] = write_equalities(dataset_id, data)
        write_ratios(dataset_id, data)
    metrics = [_per_table_metrics(nr) for nr in ncbi_rows]
    write_deficit_models(metrics)
    write_keto_amino_balance(metrics)
    write_service_codon_balance_comparison(metrics)
    print("Done. 18 files written.")
    return counts


if __name__ == "__main__":
    main()
