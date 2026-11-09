"""Fake Deimos output folder, to test the biological reading (ex-InSight) on any design:
ProteomicAnalysis_Results workbook (Differential_Expression, Log2_Impute, _manifest,
GO_*, PCA_Scores) and matrices/ folder (samples.tsv, log2_imputed.tsv). Each element
can be left out to check that the figures switch themselves off cleanly."""
from __future__ import annotations

import os
import string

import numpy as np
import pandas as pd
from scipy import stats

L = string.ascii_uppercase


def _acc(i):
    return f"P{i % 10}{L[(i // 10) % 26]}{L[(i // 260) % 26]}{(i // 6760) % 10}{i % 7}"


def _bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def make_deimos_run(root, conditions, n_rep=3, n_prot=600, k_true=4, seed=0,
                    contrasts=None, go=True, matrices=True, pca=True, manifest=True,
                    terms_per_programme=2, cond_effects=None, gsea=False):
    rng = np.random.default_rng(seed)
    os.makedirs(root, exist_ok=True)
    conds = list(conditions)
    prot = [_acc(i) for i in range(n_prot)]
    genes = [f"GENE{i}" if rng.random() < 0.6 else "" for i in range(n_prot)]
    desc = [f"Protein family {i % 37} member" for i in range(n_prot)]
    cl = rng.choice(np.arange(k_true + 1), size=n_prot,
                    p=[0.6] + [0.4 / k_true] * k_true)
    eff = {c: rng.choice([-2.0, -1.2, 0, 1.2, 2.0], size=len(conds)) for c in range(1, k_true + 1)}
    eff[0] = np.zeros(len(conds))
    if cond_effects:                    # cluster -> effect per condition (factorial tests)
        eff.update({k: np.asarray(v, float) for k, v in cond_effects.items()})
    base = rng.uniform(8, 14, n_prot)
    labels, cond_of = [], []
    for c in conds:
        for r in range(1, n_rep + 1):
            labels.append(f"{c}_r{r}")
            cond_of.append(c)
    X = np.empty((n_prot, len(labels)))
    for j, c in enumerate(cond_of):
        X[:, j] = base + np.array([eff[k][conds.index(c)] for k in cl]) + rng.normal(0, 0.3, n_prot)
    A = pd.DataFrame(X, index=prot, columns=labels)

    contrasts = contrasts or [(c, conds[0]) for c in conds[1:]]
    de = pd.DataFrame({"name": prot, "Protein.Group": prot, "Genes": genes,
                       "First.Protein.Description": desc})
    cnames = []
    for num, den in contrasts:
        cn = f"{num}_vs_{den}"
        cnames.append(cn)
        a = A[[l for l, c in zip(labels, cond_of) if c == num]]
        b = A[[l for l, c in zip(labels, cond_of) if c == den]]
        t = stats.ttest_ind(a, b, axis=1)
        de[cn + "_diff"] = (a.mean(axis=1) - b.mean(axis=1)).values
        de[cn + "_p.val"] = t.pvalue
        de[cn + "_p.adj"] = _bh(t.pvalue)

    li = A.round(4).rename_axis("Protein.Group").reset_index()
    sheets = {"Differential_Expression": de, "Log2_Impute": li}
    man = [("schema", "name", "deimos-results"), ("stats", "volcano_rule",
                                                  "p.adj < 0.05 & ratio >= 2.0")]
    man += [("contrast", cn, f"{n}|{d}") for cn, (n, d) in zip(cnames, contrasts)]
    if go:
        for i, cn in enumerate(cnames):
            sig = (de[cn + "_p.adj"] < 0.05) & (de[cn + "_diff"].abs() >= 1)
            rows = []
            for k in range(1, k_true + 1):
                members = de.loc[sig.values & (cl == k), "name"].tolist()
                if len(members) < 4:
                    continue
                up = int((de.set_index("name").loc[members, cn + "_diff"] > 0).sum())
                z = (2 * up - len(members)) / np.sqrt(len(members))
                half = len(members) // 2
                for v in range(terms_per_programme):
                    # terms 0 and 1: the whole programme; beyond, halves -> several
                    # redundancy groups per programme
                    src = ("GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC")[v % 5]
                    pool = (members if v < 2 or v % 3 == 0 else
                            members[:half] if v % 3 == 1 else members[half:])
                    frac = 0.85 if v < 2 else rng.uniform(0.55, 0.95)
                    sub = sorted(rng.choice(pool, size=max(3, int(frac * len(pool))),
                                            replace=False))
                    nm = ("process" if v == 0 else "compartment" if v == 1 else f"function {v}")
                    rows.append(dict(source=src, term_id=f"GO:{k:03d}{v:04d}",
                                     term_name=f"programme {k} {nm}",
                                     p_value=10 ** -rng.uniform(3, 9), intersection_size=len(sub),
                                     intersection=",".join(sub), z_score=z))
            if rows:
                sh = f"GO_{i + 1}"
                sheets[sh] = pd.DataFrame(rows)
                man.append(("go_sheet", sh, cn))
    gsea_sheets = {}
    if gsea:                    # GSEA_* sheets in the Deimos layout (table from column P)
        for i, cn in enumerate(cnames):
            sig = (de[cn + "_p.adj"] < 0.05) & (de[cn + "_diff"].abs() >= 1)
            rows = []
            for k in range(1, k_true + 1):
                members = de.loc[sig.values & (cl == k), "name"].tolist()
                if len(members) < 6:
                    continue
                nes = float(np.sign(de.set_index("name").loc[members, cn + "_diff"].mean()) * 2.0)
                for v in range(3):          # 3 redundant pathways per programme
                    sub = sorted(rng.choice(members, size=max(4, int((0.9 - 0.1 * v) * len(members))),
                                            replace=False))
                    rows.append({"Name": "prerank",
                                 "Term": f"GO:BP__GO:{k:03d}{v + 50:04d} pathway {k}.{v}",
                                 "ES": nes / 3, "NES": nes - 0.1 * v, "NOM p-val": 0.0,
                                 "FDR q-val": 10 ** -rng.uniform(1.5, 3), "FWER p-val": 0.0,
                                 "Tag %": f"{len(sub)}/{len(members)}", "Gene %": "10%",
                                 "Lead_genes": ";".join(de.set_index("name").loc[sub, "Genes"]
                                                        .replace("", np.nan).dropna()),
                                 "Lead_proteins": ";".join(sub)})
            if rows:
                sh = f"GSEA_{cn.replace('_vs_', 'v').replace('_', '')}"[:31]
                gsea_sheets[sh] = pd.DataFrame(rows)
                man.append(("gsea_sheet", sh, cn))
    if pca:
        Xc = (A - A.mean(axis=1).values[:, None]).T.values
        U, S, _ = np.linalg.svd(Xc - Xc.mean(0), full_matrices=False)
        var = S ** 2 / (S ** 2).sum() * 100
        sheets["PCA_Scores"] = pd.DataFrame({
            "label": labels, "condition": cond_of, "PC1": U[:, 0] * S[0], "PC2": U[:, 1] * S[1],
            "PC1_var_pct": var[0], "PC2_var_pct": var[1]})
    if manifest:
        sheets["_manifest"] = pd.DataFrame(man, columns=["section", "key", "value"])
    xlsx = os.path.join(root, "ProteomicAnalysis_Results_test.xlsx")
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        for n, d in sheets.items():
            if n == "_manifest":
                continue
            d.to_excel(w, sheet_name=n, index=False)
        for n, d in gsea_sheets.items():
            d.to_excel(w, sheet_name=n, index=False, startcol=15)
            w.sheets[n].cell(row=1, column=1, value="GSEA pre-ranked (gseapy), 100 permutations, "
                                                    "FDR <= 0.25 · synthetic test sets")
        if "_manifest" in sheets:
            sheets["_manifest"].to_excel(w, sheet_name="_manifest", index=False)
    if matrices:
        md = os.path.join(root, "matrices")
        os.makedirs(md, exist_ok=True)
        pd.DataFrame({"label": labels, "condition": cond_of,
                      "replicate": [int(l.rsplit("_r", 1)[1]) for l in labels]}) \
            .to_csv(os.path.join(md, "samples.tsv"), sep="\t", index=False)
        A.rename_axis("protein").reset_index().to_csv(os.path.join(md, "log2_imputed.tsv"),
                                                      sep="\t", index=False)
    return root
