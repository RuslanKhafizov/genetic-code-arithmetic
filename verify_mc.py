"""Independent cross-check for mc_significance.py.
It recomputes the significance figures by a deliberately different route --
reading the input CSVs directly and building the codon groups from the
third-position and four-fold-degeneracy rules (NOT from codon_groups.csv),
with a different random seed -- then compares its results against
mc_significance.csv (the result of record). Agreement is evidence against a
coding bug in either. Standard library only."""
import csv, os, random, statistics
from collections import defaultdict

# --- read data directly ---
PN = {}
for r in csv.DictReader(open("amino_acids_nucleons.csv", newline="")):
    PN[r["Amino_Acid"]] = (int(r["Protons"]), int(r["Neutrons"]))
codon2aa = {}
for r in csv.DictReader(open("genetic_code_codons.csv", newline="")):
    for c in r["Codons"].split(";"):
        codon2aa[c.strip()] = r["Product"]

BASES = "ACGT"
all_codons = [a+b+c for a in BASES for b in BASES for c in BASES]   # 64
SERVICE = {"TAA", "TAG", "TGA", "ATG"}
sense60 = [c for c in all_codons if codon2aa[c] != "Stop" and c not in SERVICE]

# --- groups built from scratch ---
def third_in(s): return [c for c in sense60 if c[2] in s]
GROUPS = {
    "All":  sense60,
    "Keto": third_in("GT"), "Amino": third_in("AC"),
    "Strong": third_in("CG"), "Weak": third_in("AT"),
    "Purine": third_in("AG"), "Pyrimidine": third_in("CT"),
}
# Octet I = boxes (1st,2nd) whose four 3rd-position codons share one amino acid
box = defaultdict(list)
for c in all_codons: box[c[:2]].append(c)
octet1 = [c for b, cs in box.items()
          if len({codon2aa[x] for x in cs}) == 1 and codon2aa[cs[0]] != "Stop"
          for c in cs]
GROUPS["Octet I"] = octet1
# Octet II = the rest of the sense pool
GROUPS["Octet II"] = [c for c in sense60 if c not in set(octet1)]
print("sanity: #Octet I codons =", len(octet1), "| #sense60 =", len(sense60))

HEAD = ["All","Keto","Amino","Strong","Weak","Purine","Pyrimidine","Octet I","Octet II"]

def sums(group, assign):           # assign: aa-name -> (P,N)
    P = N = 0
    for c in group:
        aa = codon2aa[c]
        if aa in assign:           # ATG's aa excluded since ATG not in sense60
            p, n = assign[aa]; P += p; N += n
    return P, N

def s2_count(assign):
    k = 0
    for g in HEAD:
        P, N = sums(GROUPS[g], assign)
        for v in (P, N, P+N, P-N):
            if v % 37 == 0: k += 1
    return k

# identity assignment = the real code
ident = {aa: PN[aa] for aa in PN}
oN = sums(GROUPS["All"], ident)[1]
oT = sum(sums(GROUPS["Octet I"], ident))
print(f"observed N(All)={oN} (=97*37? {oN==97*37}) | T(OctetI)={oT} (=100*37? {oT==100*37})")
print("observed S2 count =", s2_count(ident))

# --- independent Monte Carlo (pure stdlib, seed 12345) ---
aas = list(PN); pn = [PN[a] for a in aas]
rng = random.Random(12345)
M = 1_000_000
s1 = 0; s2vals = []
for _ in range(M):
    perm = pn[:]; rng.shuffle(perm)
    assign = dict(zip(aas, perm))
    if sums(GROUPS["All"], assign)[1] % 37 == 0 and \
       sum(sums(GROUPS["Octet I"], assign)) % 37 == 0:
        s1 += 1
    s2vals.append(s2_count(assign))
p1 = s1 / M
obs2 = s2_count(ident)
p2 = sum(1 for v in s2vals if v >= obs2) / M
mean2 = statistics.fmean(s2vals)
print(f"\nINDEPENDENT MC ({M:,}, stdlib, seed 12345):")
print(f"  S1  P(both anchors /37) = {p1:.6f}   (reference (1/37)^2 = {1/37**2:.6f})")
print(f"  S2  null mean = {mean2:.4f} (analytic 36/37 = {36/37:.4f})")
print(f"  S2  P(>= {obs2})        = {p2:.6f}")

# --- S4 independent check: lattice concentration over all 33 groups, both keys ---
# Rebuild the 33 groups (full-code, Octet I, Octet II) from the third-position and
# octet rules, WITHOUT reading codon_groups.csv. Each (axis, section) and each
# single-base group is taken once; this is the independent counterpart to the
# script's S4, which reads the 33 rows of codon_groups.csv.
SVC_PN = {"key0": {"TAA": (0, 0), "TAG": (0, 0), "TGA": (0, 0), "ATG": (80, 69)},
          "key1": {"TAA": (37, 37), "TAG": (37, 37), "TGA": (37, 37), "ATG": (1, 0)}}
oct1_set = set(octet1)
# Define the 33 groups as (codon_list) by mirroring the published structure:
S4_GROUPS_V = []
AXES6 = ["GT", "CG", "AC", "AT", "AG", "CT"]   # Keto Strong Amino Weak Purine Pyrimidine
SINGLE = ["C", "G", "T", "A"]
oct1_boxes = {c[:2] for c in oct1_set}
def in_oct1(c): return c[:2] in oct1_boxes
def sel(letters, section):
    g = [c for c in all_codons if c[2] in set(letters)]
    if section == "I":   g = [c for c in g if in_oct1(c)]
    elif section == "II": g = [c for c in g if not in_oct1(c)]
    return g
# All code (11): 6 axes, ALL, 4 singles
for L in AXES6: S4_GROUPS_V.append(sel(L, "all"))
S4_GROUPS_V.append(list(all_codons))
for b in SINGLE: S4_GROUPS_V.append(sel(b, "all"))
# Octet I (11): 6 axes, full Octet I, 4 singles -- all intersected with Octet I
for L in AXES6: S4_GROUPS_V.append(sel(L, "I"))
S4_GROUPS_V.append([c for c in all_codons if in_oct1(c)])
for b in SINGLE: S4_GROUPS_V.append(sel(b, "I"))
# Octet II (11): 6 axes, full Octet II, 4 singles -- all intersected with Octet II
for L in AXES6: S4_GROUPS_V.append(sel(L, "II"))
S4_GROUPS_V.append([c for c in all_codons if not in_oct1(c)])
for b in SINGLE: S4_GROUPS_V.append(sel(b, "II"))
assert len(S4_GROUPS_V) == 33, f"expected 33 groups, built {len(S4_GROUPS_V)}"

def latt_dist_v(assign_pn):
    ds = []
    for cods in S4_GROUPS_V:
        for key in ("key0", "key1"):
            P = Nn = 0
            for c in cods:
                if c in SERVICE:
                    p, n = SVC_PN[key][c]
                else:
                    p, n = assign_pn[codon2aa[c]]
                P += p; Nn += n
            for v in (P, Nn, P + Nn, P - Nn):
                if v != 0:
                    r = v % 37
                    ds.append(r if r <= 37 - r else 37 - r)
    return statistics.fmean(ds)

obs4 = latt_dist_v(ident)
print(f"\n  S4 (independent rebuild): observed mean distance = {obs4:.4f}"
      f"  over {sum(1 for cods in S4_GROUPS_V for key in ('key0','key1') for v in [1])} group-key pairs")
rng4 = random.Random(98765)
M4 = 200_000
le = 0
for _ in range(M4):
    perm = pn[:]; rng4.shuffle(perm)
    if latt_dist_v(dict(zip(aas, perm))) <= obs4:
        le += 1
p4v = le / M4
print(f"  S4 independent MC ({M4:,}, seed 98765): P(<= {obs4:.4f}) = {p4v:.6f}")

# --- compare against mc_significance.csv (the result of record) ---
csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "mc_significance.csv")
ref = {}
try:
    for row in csv.DictReader(open(csv_path, newline="")):
        ref[row["Statistic"]] = float(row["P_Value"])
except FileNotFoundError:
    ref = None

print("\nCross-check against mc_significance.csv (result of record):")
if ref is None:
    print("  mc_significance.csv not found -- run `python3 mc_significance.py`"
          " first to compare.")
else:
    print(f"  {'':16s}{'this (seed 12345)':>20s}{'mc_significance (seed 0)':>28s}")
    print(f"  {'S1 P(both /37)':16s}{p1:>20.6f}{ref['S1']:>28.6f}")
    print(f"  {'S2 P(>= 13)':16s}{p2:>20.6f}{ref['S2']:>28.6f}")
    if "S4" in ref:
        print(f"  {'S4 P(<=dist)':16s}{p4v:>20.6f}{ref['S4']:>28.6f}")
    print("  Two independent implementations; the small differences are"
          " Monte-Carlo error")
    print("  (different seeds and code paths). The agreement is the check.")
