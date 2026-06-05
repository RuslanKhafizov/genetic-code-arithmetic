#!/usr/bin/env python3
"""
reproduce.py — Generates all 13 derived CSV files from 5 input CSV files.

Unified output rules for every file:
  - Encoding UTF-8, LF line endings, trailing newline.
  - Quoting: csv.QUOTE_NONNUMERIC (string fields quoted, numbers bare).
  - One deterministic sort rule per file, applied to every row.
"""

import csv
import os
from itertools import combinations

ACTIVE_KEYS = ["key0", "key1"]
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
        CODON_TO_AA[c.strip()] = r["Product"]

KEY_PARAMS = {}
for r in key_rows:
    KEY_PARAMS.setdefault(r["Key_ID"], {})[r["Codon"]] = {
        "P": int(r["Final_Protons"]), "N": int(r["Final_Neutrons"]),
    }

SERVICE_CODONS = {"TAA", "TAG", "TGA", "ATG"}


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


# ── 2. Nucleon data per key ────────────────────────────────────────────────

def compute_nucleon_data(key_id):
    kp = KEY_PARAMS[key_id]
    results = []
    for gr in group_rows:
        tp = tn = 0
        for codon in gr["Codon_List"].split(";"):
            codon = codon.strip()
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
            "Group_Size": int(gr["Group_Size"]),
            "Total_Nucleons": tp + tn, "Protons": tp, "Neutrons": tn,
            "Delta_P_N": tp - tn,
        })
    results.sort(key=lambda d: d["Group_ID"])
    return results


def twin_field(d):
    t = d["Twin_Group_ID"]
    return int(t) if t else ""


def write_nucleon_data(key_id, data):
    header = ["Group_Name", "Group_ID", "Twin_Group_ID", "Code_Section",
              "Group_Size", "Total_Nucleons", "Protons", "Neutrons",
              "Delta_P_N"]
    rows = [[d["Group_Name"], d["Group_ID"], twin_field(d),
             d["Code_Section"], d["Group_Size"], d["Total_Nucleons"],
             d["Protons"], d["Neutrons"], d["Delta_P_N"]] for d in data]
    write_csv(f"{key_id}_nucleon-data.csv", header, rows)


# ── 3. Divisibility by 37 ──────────────────────────────────────────────────

def div37_expr(v):
    v = int(v)
    if v < 37:
        return str(v)
    q, r = divmod(v, 37)
    return f"{q}*37" if r == 0 else f"{q}*37+{r}"


def write_divisibility(key_id, data):
    header = ["Group_Name", "Group_ID", "Twin_Group_ID", "Code_Section",
              "Group_Size", "Total_mod37", "Protons_mod37",
              "Neutrons_mod37", "Delta_mod37"]
    rows = [[d["Group_Name"], d["Group_ID"], twin_field(d),
             d["Code_Section"], d["Group_Size"],
             div37_expr(d["Total_Nucleons"]), div37_expr(d["Protons"]),
             div37_expr(d["Neutrons"]), div37_expr(d["Delta_P_N"])]
            for d in data]
    write_csv(f"{key_id}_divisibility-37.csv", header, rows)


# ── 4. Equalities ──────────────────────────────────────────────────────────

EQ_PARAMS = ["Total_Nucleons", "Protons", "Neutrons", "Delta_P_N"]
MATCH_ORDER = {"DIFF_PARAM": 0, "CROSS_PARAM": 1, "SAME_PARAM": 2}


def write_equalities(key_id, data):
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
                    if sp == tp:
                        mt = "SAME_PARAM"
                    elif not same_size:
                        mt = "DIFF_PARAM"
                    else:
                        mt = "CROSS_PARAM"
                    eq_rows.append([
                        src[sp], mt, src["Code_Section"], src["Group_ID"],
                        src["Group_Name"], sp, tp, tgt["Group_Name"],
                        tgt["Group_ID"], tgt["Code_Section"],
                    ])
    eq_rows.sort(key=lambda r: (r[0], MATCH_ORDER[r[1]], r[3], r[8],
                                EQ_PARAMS.index(r[5]),
                                EQ_PARAMS.index(r[6])))
    write_csv(f"{key_id}_equalities.csv", header, eq_rows)
    return len(eq_rows)


# ── 5. Ratios ──────────────────────────────────────────────────────────────

def write_ratios(key_id, data):
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
    write_csv(f"{key_id}_ratios.csv", header, rows)


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
    ctx = nr["Stop_Context_Dependent"].strip() == "True"
    excl_m1 = stops | {"ATG"}
    excl_m2 = {"ATG"} if ctx else excl_m1

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

    return {
        "tt": int(nr["Transl_Table"]),
        "name": nr["Code_Name"],
        "excl_m1": excl_m1, "excl_m2": excl_m2,
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


def main():
    write_differences("P", "Proton_Difference",
                      "amino_acids_proton_differences.csv")
    write_differences("N", "Neutron_Difference",
                      "amino_acids_neutron_differences.csv")
    write_differences("T", "Nucleon_Difference",
                      "amino_acids_nucleon_differences.csv")
    counts = {}
    for k in ACTIVE_KEYS:
        data = compute_nucleon_data(k)
        write_nucleon_data(k, data)
        write_divisibility(k, data)
        counts[k] = write_equalities(k, data)
        write_ratios(k, data)
    metrics = [_per_table_metrics(nr) for nr in ncbi_rows]
    write_deficit_models(metrics)
    write_keto_amino_balance(metrics)
    print("Done. 13 files written.")
    return counts


if __name__ == "__main__":
    main()
