"""Jeu synthétique conforme au schéma Deimos pour tester Pathfinder."""
import numpy as np
import pandas as pd
from itertools import combinations

rng = np.random.default_rng(7)
conds = ["Ctrl", "TreatA", "TreatB"]
reps = 4
N = 800
samples = [f"{c}_{i+1}" for c in conds for i in range(reps)]

# intensités log2 LFQ : effet condition sur un sous-ensemble de protéines
base = rng.normal(20, 2, size=(N, 1))
eff = np.zeros((N, len(conds)))
eff[:120, 1] = rng.normal(1.5, 0.4, 120)      # up TreatA
eff[120:240, 2] = rng.normal(-1.3, 0.4, 120)  # down TreatB
mat = np.zeros((N, len(samples)))
for j, s in enumerate(samples):
    ci = conds.index(s.rsplit("_", 1)[0])
    mat[:, j] = (base[:, 0] + eff[:, ci] + rng.normal(0, 0.6, N))

genes = [f"GENE{i:04d}" for i in range(N)]
accs = [f"P{i:05d}" for i in range(N)]
desc = [f"Protein {i} description" for i in range(N)]

comp = pd.DataFrame(mat, columns=samples)
comp.insert(0, "Protein.Group", accs)
comp.insert(1, "Genes", genes)
comp.insert(2, "First.Protein.Description", desc)

# stats par contraste
contrasts = ["_vs_".join(c) for c in combinations(conds, 2)]
sig_any = np.zeros(N, bool)
for c in contrasts:
    a, b = c.split("_vs_")
    ia, ib = conds.index(a), conds.index(b)
    diff = eff[:, ia] - eff[:, ib] + rng.normal(0, 0.1, N)
    p = np.clip(np.exp(-np.abs(diff) * 3) * rng.uniform(0.5, 1.5, N), 1e-6, 1)
    padj = np.clip(p * 1.5, 0, 1)
    comp[f"{c}_diff"] = diff
    comp[f"{c}_p.val"] = p
    comp[f"{c}_p.adj"] = padj
    comp[f"Robustness_Score_{c}"] = rng.integers(0, 101, N)
    comp[f"Pi_Score_{c}"] = np.abs(diff) * -np.log10(p)
    sig_any |= padj < 0.05
comp["significant"] = sig_any

# Z-score_Heatmap
z = (mat - mat.mean(axis=1, keepdims=True)) / (mat.std(axis=1, keepdims=True) + 1e-9)
zdf = pd.DataFrame(z, columns=samples)
zdf.insert(0, "Protein.Group", accs)
zdf.insert(1, "Genes", genes)
zdf.insert(2, "First.Protein.Description", desc)

with pd.ExcelWriter("rapportstatistique.xlsx") as xw:
    comp.to_excel(xw, sheet_name="Comparaison", index=False)
    zdf.to_excel(xw, sheet_name="Z-score_Heatmap", index=False)

# GO : une sheet par contraste, colonnes gprofiler + z_score
sources = ["GO:BP", "GO:CC", "GO:MF", "REAC", "KEGG"]
go_terms = [(f"GO:{1000000+i}", f"terme biologique {i}", sources[i % 5]) for i in range(60)]
with pd.ExcelWriter("rapportenrichissementGO.xlsx") as xw:
    for c in contrasts:
        abbr = c.replace("_vs_", "v")
        rows = []
        k = rng.integers(15, 30)
        chosen = rng.choice(len(go_terms), size=k, replace=False)
        for ti in chosen:
            tid, tname, src = go_terms[ti]
            zval = rng.normal(0, 1.5)
            inter = ",".join(rng.choice(accs[:300], size=rng.integers(3, 12), replace=False))
            rows.append({"source": src, "term_id": tid, "term_name": tname,
                         "p_value": np.clip(rng.exponential(0.02), 1e-6, 0.2),
                         "intersection": inter, "z_score": zval})
        pd.DataFrame(rows).to_excel(xw, sheet_name=f"GO_{abbr}", index=False)

print("Synthétique généré :", samples, "|", contrasts)
