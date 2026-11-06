"""
Meta-GO axis - cross-contrast comparison at the GO/pathway level.

Builds a whole-project view by crossing every contrast on its GO/pathway terms.
Core object: a `term x contrast` matrix of signed z-scores, masked by significance
(NaN = not significant / not tested, never 0 - a 0 would distort clustering).

- Cross-directionality: sign(z) per contrast -> shared / specific / mixed pattern.
- Redundancy: terms clustered by Jaccard distance on their `intersection` accession
  sets (offline, no ontology download). Collapses terms hitting the same proteins
  and elects a representative.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import squareform

# GO / pathway terms are significant at q < 0.05: the threshold of the upstream
# enrichment (gProfiler g:SCS/FDR, STRING FDR) is the one Pathfinder keeps. It does
# not follow the protein rule (--pthr / --min-ratio): the terms are read, not retested.
GO_Q_THR = 0.05


def build_matrix(go_long: pd.DataFrame, p_thr: float = GO_Q_THR):
    """term x contrast matrix (signed z, masked by p_thr) + per-term metadata."""
    df = go_long.copy()
    df["sig"] = df["p_value"] < p_thr
    df["cell"] = np.where(df["sig"], df["z"], np.nan)

    # unique term key = term_id when available, else term_name+source
    df["key"] = np.where(
        df["term_id"].astype(str).str.len() > 0,
        df["term_id"].astype(str),
        df["term_name"].astype(str) + "|" + df["source"].astype(str),
    )

    matrix = df.pivot_table(index="key", columns="contrast", values="cell", aggfunc="mean")
    # keep contrasts that were read but hold no significant term: an all-NaN column is
    # honest ("tested, nothing passed"), whereas dropping it looks like a missing sheet.
    all_contrasts = list(dict.fromkeys(df["contrast"].dropna().tolist()))
    matrix = matrix.reindex(columns=all_contrasts)

    meta = (df.sort_values("p_value")
              .drop_duplicates("key")
              .set_index("key")[["term_name", "term_id", "source"]])
    # number of contrasts where the term is significant
    meta["n_sig"] = matrix.notna().sum(axis=1)
    meta = meta.loc[matrix.index]
    return matrix, meta, df


def directionality(matrix: pd.DataFrame, meta: pd.DataFrame, pairs=None) -> pd.DataFrame:
    """Per term: up / down contrasts and a compact pattern label.

    pairs (axis_proteins.comparable_pairs): 'opposite' = term reversed between two
    comparable contrasts (same comparison in different strata); 'mixed' = up in some
    contrasts and down in others that compare different things. Without pairs, any
    term both up and down is 'mixed' (legacy)."""
    from axis_proteins import _reversals
    recs = []
    for key, row in matrix.iterrows():
        up = [c for c in matrix.columns if pd.notna(row[c]) and row[c] > 0]
        dn = [c for c in matrix.columns if pd.notna(row[c]) and row[c] < 0]
        rev = (_reversals({c: float(row[c]) for c in up + dn}, pairs)
               if pairs and up and dn else [])
        if rev:
            pattern = "opposite"
        elif up and dn:
            pattern = "mixed"
        elif up:
            pattern = "shared up" if len(up) > 1 else "specific up"
        elif dn:
            pattern = "shared down" if len(dn) > 1 else "specific down"
        else:
            pattern = "—"
        recs.append({
            "term_name": meta.loc[key, "term_name"],
            "source": meta.loc[key, "source"],
            "n_sig": int(meta.loc[key, "n_sig"]),
            "pattern": pattern,
            "up_in": ", ".join(up),
            "down_in": ", ".join(dn),
            "opposite_in": "; ".join(rev),
        })
    out = pd.DataFrame(recs, index=matrix.index)
    return out.sort_values(["n_sig", "pattern"], ascending=[False, True])


def _jaccard_condensed(sets):
    n = len(sets)
    d = np.zeros(n * (n - 1) // 2)
    k = 0
    for i in range(n):
        si = sets[i]
        for j in range(i + 1, n):
            sj = sets[j]
            if not si and not sj:
                dist = 1.0
            else:
                inter = len(si & sj)
                union = len(si | sj)
                dist = 1.0 - (inter / union if union else 0.0)
            d[k] = dist
            k += 1
    return d


def redundancy_clusters(go_long: pd.DataFrame, matrix: pd.DataFrame, meta: pd.DataFrame,
                        jaccard_cut: float = 0.5):
    """Group redundant terms (same proteins) and elect a representative.

    jaccard_cut: cut distance (0.5 -> terms sharing >50% of their genes are merged).
    Only terms significant in >=1 contrast are considered (index of `matrix`).
    """
    df = go_long.copy()
    df["key"] = np.where(
        df["term_id"].astype(str).str.len() > 0,
        df["term_id"].astype(str),
        df["term_name"].astype(str) + "|" + df["source"].astype(str),
    )
    # accession set per term = union of intersections across all contrasts
    gene_sets = {}
    for key, grp in df.groupby("key"):
        s = set()
        for lst in grp["intersection"]:
            s.update(lst)
        gene_sets[key] = s

    keys = list(matrix.index)
    meta = meta.copy()
    meta["redundancy_cluster"] = -1
    meta["is_representative"] = True

    sets = [gene_sets.get(k, set()) for k in keys]
    if len(keys) >= 2 and any(sets):
        dcond = _jaccard_condensed(sets)
        Z = linkage(dcond, method="average")
        clusters = fcluster(Z, t=jaccard_cut, criterion="distance")
        meta.loc[keys, "redundancy_cluster"] = clusters
        # cluster representative = strongest mean |z|
        strength = matrix.abs().mean(axis=1).fillna(0)
        for cl in np.unique(clusters):
            members = [k for k, c in zip(keys, clusters) if c == cl]
            rep = max(members, key=lambda k: (strength.get(k, 0), meta.loc[k, "n_sig"]))
            for k in members:
                meta.loc[k, "is_representative"] = (k == rep)
    return meta


def _pair_relation(c1, c2):
    """Describe how two contrasts share conditions — shared references create
    mechanical correlation (same reference) or anti-correlation (flipped role),
    which must not be read as biology."""
    s1, s2 = str(c1), str(c2)
    if "_vs_" not in s1 or "_vs_" not in s2:
        return "disjoint", ""
    a1, b1 = s1.split("_vs_", 1)
    a2, b2 = s2.split("_vs_", 1)
    shared = {a1, b1} & {a2, b2}
    if not shared:
        return "disjoint", ""
    if len(shared) == 2:
        return "same pair", "same two conditions"
    cond = shared.pop()
    same_side = (cond == a1 and cond == a2) or (cond == b1 and cond == b2)
    if same_side:
        return "shared reference", f"both compare against {cond} on the same side"
    return "flipped reference", f"{cond} is on opposite sides"


def contrast_similarity(matrix: pd.DataFrame, min_shared: int = 3):
    """Cluster contrasts by their GO signature.

    Similarity = Pearson correlation between the signed z-score profiles of two
    contrasts, computed on the terms significant in BOTH (pairwise-complete). A
    high positive value means the two contrasts move the same functions the same
    way; a negative one means they move them in opposite directions.

    Returns (corr, order, n_shared) where `order` is the hierarchical leaf order.
    NaN where fewer than `min_shared` terms are shared (not comparable).
    """
    cols = list(matrix.columns)
    n = len(cols)
    corr = pd.DataFrame(np.nan, index=cols, columns=cols, dtype=float)
    shared = pd.DataFrame(0, index=cols, columns=cols, dtype=int)
    for i, a in enumerate(cols):
        corr.loc[a, a] = 1.0
        shared.loc[a, a] = int(matrix[a].notna().sum())
        for b in cols[i + 1:]:
            both = matrix[[a, b]].dropna()
            shared.loc[a, b] = shared.loc[b, a] = len(both)
            if len(both) >= min_shared and both[a].std() > 0 and both[b].std() > 0:
                r = float(both[a].corr(both[b]))
                corr.loc[a, b] = corr.loc[b, a] = r

    order = list(range(n))
    if n >= 3:
        try:
            d = 1.0 - corr.fillna(0.0).values          # NaN -> no similarity, for ordering only
            np.fill_diagonal(d, 0.0)
            d = (d + d.T) / 2.0
            Z = linkage(squareform(np.clip(d, 0, 2), checks=False), method="average")
            order = list(leaves_list(Z))
        except Exception:
            pass
    labels = [cols[i] for i in order]
    rel = pd.DataFrame("", index=labels, columns=labels, dtype=object)
    for a in labels:
        for b in labels:
            rel.loc[a, b] = _pair_relation(a, b)[0] if a != b else "self"
    return corr.loc[labels, labels], labels, shared.loc[labels, labels], rel


def run(go_long: pd.DataFrame, p_thr: float = GO_Q_THR, jaccard_cut: float = 0.5, pairs=None):
    matrix, meta, df = build_matrix(go_long, p_thr)
    meta = redundancy_clusters(go_long, matrix, meta, jaccard_cut)
    direction = directionality(matrix, meta, pairs)
    # accession sets per term (union of significant intersections) - for the GO UMAP
    sig = go_long[go_long["p_value"] < p_thr].copy()
    sig["key"] = np.where(sig["term_id"].astype(str).str.len() > 0,
                          sig["term_id"].astype(str),
                          sig["term_name"].astype(str) + "|" + sig["source"].astype(str))
    term_genes = {}
    for key, grp in sig.groupby("key"):
        s = set()
        for lst in grp["intersection"]:
            s.update(lst)
        term_genes[key] = s
    # contrast clustering on the GO signature
    ccorr, corder, cshared, crel = contrast_similarity(matrix)
    return {"matrix": matrix, "meta": meta, "directionality": direction,
            "term_genes": term_genes, "contrast_corr": ccorr,
            "contrast_order": corder, "contrast_shared": cshared,
            "contrast_relation": crel, "pairs": pairs}
