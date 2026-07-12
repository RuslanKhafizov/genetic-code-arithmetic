#!/usr/bin/env python3
"""
mc_significance.py -- Simple, pre-specified Monte Carlo significance test.

NULL MODEL (one model, stated once): the "random genetic code" of the
Freeland-Hurst type (Freeland & Hurst 1998), restricted to relabelling. The
architecture of the standard code is held fixed -- the 64 codons, which codons
are synonyms (the 20 amino-acid blocks with their sizes), the three stop
codons, the Rumer Octet I/II partition and the third-position chemical axes --
and the ONLY thing randomized is which amino acid (carrying its fixed
proton/neutron counts P,N) occupies each of the 20 blocks. A trial is thus a
uniformly random permutation of the 20 amino acids over the 20 blocks. Nothing
is "constructed": the amino-acid P,N values are the real ones, 37 is taken a
priori from the prior literature, and the test statistics are fixed in advance.

SENSE POOL: the 60 codons excluding the three stops and the initiator ATG,
matching the preprint's definition (so N(All sense)=3589=97*37 for the
standard code).

SERVICE CODONS: S1-S3 use the 60-codon sense pool, so the amino acid a
permutation places at the ATG block does not contribute. The lattice rows S4a
and S4b hold the service-codon contributions fixed at the Key 0 and Key 1
values in every trial -- they are fixed-key conditional diagnostics, not full
relabellings at the service positions.

PRE-SPECIFIED STATISTICS (decided before running; all reported, none selected
post hoc):
  S1  two independent 37-anchors: N(All sense) and T(Octet I) BOTH divisible
      by 37.
  S2  37-saturation of the headline partition: among the four nucleon
      quantities (T,P,N,Delta) of the nine headline sense-pool groups
      {All, Keto, Amino, Strong, Weak, Purine, Pyrimidine, Octet I, Octet II},
      the number divisible by 37. p = P(random code reaches >= the standard
      code's count).
  S3  same count, but over the parametrization-independent groups only (those
      containing no service codon, on which 37 can never be imposed) -- the
      most conservative, fully non-circular variant.

REPRODUCIBILITY AND DEPENDENCIES: standard library only (no third-party
packages). The test is randomized rather than exact, but seeded (default seed
0, 1,000,000 trials), so the figures it reports are reproducible run to run. A
full run of 1,000,000 trials takes about two minutes in pure Python (the S4
lattice statistic dominates the runtime). It writes a summary table
mc_significance.csv. An independent
re-implementation, verify_mc.py, reproduces these figures by a different route.
"""
import csv, os, random, statistics, collections

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def read_csv(filename):
    with open(os.path.join(SCRIPT_DIR, filename), newline="",
              encoding="utf-8") as f:
        return list(csv.DictReader(f))

# --- data ---
aa_rows = read_csv("amino_acids_nucleons.csv")
codon_rows = read_csv("genetic_code_codons.csv")
group_rows = read_csv("codon_groups.csv")

PN = {r["Amino_Acid"]: (int(r["Protons"]), int(r["Neutrons"])) for r in aa_rows}
aas = sorted(PN)
idx = {a: i for i, a in enumerate(aas)}
P = [PN[a][0] for a in aas]
N = [PN[a][1] for a in aas]

codon2aa = {}
for r in codon_rows:
    for c in r["Codons"].split(";"):
        codon2aa[c.strip()] = r["Product"]

SERVICE = {"TAA", "TAG", "TGA", "ATG"}
cluster_of = {c: idx[a] for c, a in codon2aa.items()
              if a in idx and c not in SERVICE}     # sense-pool-60 codons only

groups = {g["Group_Name"]: [c.strip() for c in g["Codon_List"].split(";")]
          for g in group_rows}

def wrow_nz(codons):
    """Sparse weight row: [(block_index, count), ...] over sense-pool codons."""
    v = [0] * 20
    for c in codons:
        if c in cluster_of:
            v[cluster_of[c]] += 1
    return [(i, w) for i, w in enumerate(v) if w]

ALL  = wrow_nz(groups["ALL: {C, G, A, T}"])
OCT1 = wrow_nz(groups["Octet I: {C, G, A, T}"])

HEAD = ["ALL: {C, G, A, T}", "Keto: {G, T}", "Amino: {A, C}", "Strong: {C, G}",
        "Weak: {A, T}", "Purine: {A, G}", "Pyrimidine: {C, T}",
        "Octet I: {C, G, A, T}", "Octet II: {C, G, A, T}"]
W_head = [wrow_nz(groups[g]) for g in HEAD]

# Parametrization-independent groups: those with no service codon. Build each
# weight row from the row's OWN Codon_List, iterating group_rows directly --
# NOT via the `groups` dict, whose keys (Group_Name) are non-unique across code
# sections ({C}, {G}, {T}, {A} and the two-letter axis names each recur at the
# All / Octet I / Octet II levels). A name lookup would silently collapse those
# to the last (Octet II) occurrence and compute S3 over the wrong groups.
W_indep = [wrow_nz([c.strip() for c in g["Codon_List"].split(";")])
           for g in group_rows
           if not (set(c.strip() for c in g["Codon_List"].split(";")) & SERVICE)]

# --- S4 setup: lattice concentration over all 33 groups under BOTH keys ---
# Each group contributes its sense-codon block counts (as a sparse weight row)
# plus a fixed service-codon contribution that depends on the key. Service
# codons (3 stops + ATG) carry (0,0) under Key 0's stops / (80,69) for ATG, and
# (37,37) per stop / (1,0) for ATG under Key 1.
SERVICE_PN = {
    "key0": {"TAA": (0, 0), "TAG": (0, 0), "TGA": (0, 0), "ATG": (80, 69)},
    "key1": {"TAA": (37, 37), "TAG": (37, 37), "TGA": (37, 37), "ATG": (1, 0)},
}
# Iterate group_rows directly: groups with the same name but different
# Code_Section (e.g. {C} over the full code vs within Octet I vs Octet II) have
# DIFFERENT codon lists, so keying by name alone would silently drop them. Each
# of the 33 rows is taken as a distinct group.
S4_GROUPS = []
for grow in group_rows:
    codons_g = [c.strip() for c in grow["Codon_List"].split(";")]
    nz = wrow_nz(codons_g)
    svc = {}
    for key in ("key0", "key1"):
        sp = sn = 0
        for c in codons_g:
            if c in SERVICE:
                cp, cn = SERVICE_PN[key][c]
                sp += cp; sn += cn
        svc[key] = (sp, sn)
    S4_GROUPS.append((nz, svc))

def lattice_mean_distance(Pv, Nv):
    """Mean distance to the nearest multiple of 37 over the nonzero T,P,N,Delta
    of all 33 groups under both keys (service contributions fixed per key)."""
    dsum = 0; dcount = 0
    for nz, svc in S4_GROUPS:
        base_p = dot(nz, Pv); base_n = dot(nz, Nv)
        for key in ("key0", "key1"):
            sp, sn = svc[key]
            p = base_p + sp; n = base_n + sn
            for v in (p, n, p + n, p - n):
                if v != 0:
                    r = v % 37
                    dsum += r if r <= 37 - r else 37 - r
                    dcount += 1
    return dsum / dcount

def lattice_mean_distance_key(Pv, Nv, which):
    """S4 restricted to a single key ('key0' or 'key1'): mean distance to the
    nearest multiple of 37 over the nonzero T,P,N,Delta of the 33 groups under
    that key alone. S4a = key0 (no divisibility by 37 imposed), S4b = key1."""
    dsum = 0; dcount = 0
    for nz, svc in S4_GROUPS:
        base_p = dot(nz, Pv); base_n = dot(nz, Nv)
        sp, sn = svc[which]
        p = base_p + sp; n = base_n + sn
        for v in (p, n, p + n, p - n):
            if v != 0:
                r = v % 37
                dsum += r if r <= 37 - r else 37 - r
                dcount += 1
    return dsum / dcount

def dot(nz, vec):
    s = 0
    for i, w in nz:
        s += w * vec[i]
    return s

def div37_count(rows, Pv, Nv):
    k = 0
    for nz in rows:
        p = dot(nz, Pv); n = dot(nz, Nv)
        if p % 37 == 0: k += 1
        if n % 37 == 0: k += 1
        if (p + n) % 37 == 0: k += 1
        if (p - n) % 37 == 0: k += 1
    return k

# --- observed (true assignment, perm = identity) ---
obs_Nall = dot(ALL, N); obs_TocI = dot(OCT1, P) + dot(OCT1, N)
obs_S2 = div37_count(W_head, P, N)
obs_S3 = div37_count(W_indep, P, N)
obs_S4 = lattice_mean_distance(P, N)
obs_S4a = lattice_mean_distance_key(P, N, "key0")
obs_S4b = lattice_mean_distance_key(P, N, "key1")
print("=== Standard code (observed) ===")
print(f"N(All sense) = {obs_Nall}  = {obs_Nall//37}*37 + {obs_Nall%37}")
print(f"T(Octet I)   = {obs_TocI}  = {obs_TocI//37}*37 + {obs_TocI%37}")
print(f"S1 both anchors divisible by 37: {obs_Nall%37==0 and obs_TocI%37==0}")
print(f"S2 headline /37 count    = {obs_S2}  of {len(W_head)*4}")
print(f"S3 param-indep /37 count = {obs_S3}  of {len(W_indep)*4}")
print(f"S4 mean distance to 37-lattice (33 groups x 2 keys) = {obs_S4:.4f}")
print(f"S4a mean distance, Key 0 only = {obs_S4a:.4f}")
print(f"S4b mean distance, Key 1 only = {obs_S4b:.4f}")

# --- Monte Carlo (standard library random) ---
M = 1_000_000
SEED = 0
rng = random.Random(SEED)
order = list(range(20))
s1 = 0
s2_sum = 0; s2_ge = 0; s2_hist = collections.Counter()
s3_sum = 0; s3_ge = 0
s4_sum = 0.0; s4_le = 0
s4a_sum = 0.0; s4a_le = 0
s4b_sum = 0.0; s4b_le = 0
for _ in range(M):
    perm = order[:]; rng.shuffle(perm)
    Pp = [P[i] for i in perm]; Np = [N[i] for i in perm]
    if dot(ALL, Np) % 37 == 0 and (dot(OCT1, Pp) + dot(OCT1, Np)) % 37 == 0:
        s1 += 1
    c2 = div37_count(W_head, Pp, Np); s2_sum += c2; s2_hist[c2] += 1
    if c2 >= obs_S2: s2_ge += 1
    c3 = div37_count(W_indep, Pp, Np); s3_sum += c3
    if c3 >= obs_S3: s3_ge += 1
    d4 = lattice_mean_distance(Pp, Np); s4_sum += d4
    if d4 <= obs_S4: s4_le += 1
    d4a = lattice_mean_distance_key(Pp, Np, "key0"); s4a_sum += d4a
    if d4a <= obs_S4a: s4a_le += 1
    d4b = lattice_mean_distance_key(Pp, Np, "key1"); s4b_sum += d4b
    if d4b <= obs_S4b: s4b_le += 1

p1 = s1 / M; p2 = s2_ge / M; p3 = s3_ge / M; p4 = s4_le / M
p4a = s4a_le / M; p4b = s4b_le / M
print(f"\n=== Monte Carlo, {M:,} random codes (stdlib random, seed {SEED}) ===")
print(f"S1  P(both anchors /37) = {p1:.6f}   "
      f"[independent reference (1/37)^2 = {1/37**2:.6f}]")
print(f"S2  observed {obs_S2}; null mean {s2_sum/M:.4f}; P(>= {obs_S2}) = {p2:.6f}")
print(f"S3  observed {obs_S3}; null mean {s3_sum/M:.4f}; P(>= {obs_S3}) = {p3:.6f}")
print(f"S4  observed {obs_S4:.4f}; null mean {s4_sum/M:.4f}; "
      f"P(<= {obs_S4:.4f}) = {p4:.6f}")
print(f"S4a observed {obs_S4a:.4f}; null mean {s4a_sum/M:.4f}; "
      f"P(<= {obs_S4a:.4f}) = {p4a:.6f}  [Key 0 only]")
print(f"S4b observed {obs_S4b:.4f}; null mean {s4b_sum/M:.4f}; "
      f"P(<= {obs_S4b:.4f}) = {p4b:.6f}  [Key 1 only]")

mx = max(s2_hist)
print("\nS2 null distribution (count : frequency):")
for k in range(mx + 1):
    print(f"  {k:>2} : {s2_hist.get(k, 0)}")

# --- summary table ---
with open(os.path.join(SCRIPT_DIR, "mc_significance.csv"),
          "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC, lineterminator="\n")
    w.writerow(["Statistic", "Description", "Observed", "Null_Mean",
                "P_Value", "N_Trials", "Seed"])
    w.writerow(["S1", "N(All sense) and T(Octet I) both divisible by 37",
                "both divisible", round(p1, 6), round(p1, 6), M, SEED])
    w.writerow(["S2",
                "count of T/P/N/Delta divisible by 37 over 9 headline "
                "sense-pool groups (of 36)",
                obs_S2, round(s2_sum / M, 4), round(p2, 6), M, SEED])
    w.writerow(["S3",
                "count of T/P/N/Delta divisible by 37 over "
                f"parametrization-independent groups (of {len(W_indep)*4})",
                obs_S3, round(s3_sum / M, 4), round(p3, 6), M, SEED])
    w.writerow(["S4",
                "mean distance to nearest multiple of 37 over T/P/N/Delta "
                "of all 33 groups under both keys (lower = more concentrated)",
                round(obs_S4, 4), round(s4_sum / M, 4), round(p4, 6), M, SEED])
    w.writerow(["S4a",
                "S4 restricted to Key 0 only (no divisibility by 37 imposed); "
                "shows the concentration is present before the derived Key 1 "
                "symmetrization (Key 0 is itself an assignment)",
                round(obs_S4a, 4), round(s4a_sum / M, 4), round(p4a, 6), M, SEED])
    w.writerow(["S4b",
                "S4 restricted to Key 1 only (the key derived under the "
                "balance and divisibility conditions); sharpens but does not "
                "create the concentration",
                round(obs_S4b, 4), round(s4b_sum / M, 4), round(p4b, 6), M, SEED])
print("\nWrote mc_significance.csv")
