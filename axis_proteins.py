"""
Protein axis - cross-contrast comparison at the protein level.

The Meta-GO axis answers "which functions do these contrasts share"; this one answers
"which proteins carry the signature" - the table that ends up in a paper.

Core object: a `protein x contrast` matrix of log2 fold-changes, masked by
significance (NaN = not significant, never 0).

Ranking uses the Pi-score (|LFC| x -log10 p, computed upstream): it is a ranking
statistic, not a testable one, so it orders proteins but never selects them - the
selection stays on the adjusted p-value. Robustness, when the upstream module provides
it, is carried through as a confidence column.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _find_contrast_col(comp, contrast, kinds):
    """Locate '<contrast>_diff', 'Pi_Score_<contrast>'... tolerating naming variants."""
    cols = {str(c).lower(): c for c in comp.columns}
    c = str(contrast).lower()
    for k in kinds:
        for cand in (f"{c}_{k}", f"{k}_{c}"):
            if cand in cols:
                return cols[cand]
    # last resort: a column containing both the contrast and the token
    for lc, orig in cols.items():
        if c in lc and any(k in lc for k in kinds):
            return orig
    return None


def build_matrix(stats, p_thr: float = 0.05):
    """protein x contrast matrix of significant log2 fold-changes + per-protein meta."""
    comp = stats.comp
    idx = comp[stats.id_col].astype(str) if stats.id_col in comp.columns \
        else pd.Series([f"row{i}" for i in range(len(comp))])

    lfc, pi, rob = {}, {}, {}
    for c in stats.contrasts:
        p_col = _find_contrast_col(comp, c, ["p.adj", "padj", "p_adj"])
        d_col = _find_contrast_col(comp, c, ["diff", "logfc", "lfc"])
        if p_col is None or d_col is None:
            continue
        p = pd.to_numeric(comp[p_col], errors="coerce")
        d = pd.to_numeric(comp[d_col], errors="coerce")
        lfc[c] = np.where(p < p_thr, d, np.nan)
        pi_col = _find_contrast_col(comp, c, ["pi_score", "piscore", "pi"])
        if pi_col is not None:
            pi[c] = np.where(p < p_thr, pd.to_numeric(comp[pi_col], errors="coerce"), np.nan)
        rob_col = _find_contrast_col(comp, c, ["robustness_score", "robustness", "robust"])
        if rob_col is not None:
            rob[c] = np.where(p < p_thr, pd.to_numeric(comp[rob_col], errors="coerce"), np.nan)

    matrix = pd.DataFrame(lfc, index=idx.values)
    keep = matrix.notna().any(axis=1).values          # proteins significant somewhere
    matrix = matrix.loc[keep]
    comp_f = comp.loc[keep]
    pi_df = pd.DataFrame(pi, index=idx.values).loc[keep] if pi else None
    rob_df = pd.DataFrame(rob, index=idx.values).loc[keep] if rob else None

    meta = pd.DataFrame(index=matrix.index)
    if stats.gene_col in comp.columns:
        meta["gene"] = comp_f[stats.gene_col].astype(str).values
    if stats.desc_col in comp.columns:
        meta["description"] = comp_f[stats.desc_col].astype(str).values
    meta["n_sig"] = matrix.notna().sum(axis=1)
    meta["n_up"] = (matrix > 0).sum(axis=1)
    meta["n_down"] = (matrix < 0).sum(axis=1)
    meta["max_abs_lfc"] = matrix.abs().max(axis=1)
    if pi_df is not None:
        meta["best_pi"] = pi_df.abs().max(axis=1)
        meta["median_pi"] = pi_df.abs().median(axis=1)
    if rob_df is not None:
        meta["median_robustness"] = rob_df.median(axis=1)
    return matrix, meta, pi_df, rob_df


def directionality(matrix: pd.DataFrame, meta: pd.DataFrame):
    """Per protein: shared / specific / opposite regulation across contrasts."""
    up = meta["n_up"].values
    dn = meta["n_down"].values
    pattern = np.where((up > 0) & (dn > 0), "opposite",
                       np.where(up > 1, "shared up",
                                np.where(dn > 1, "shared down",
                                         np.where(up == 1, "specific up", "specific down"))))
    out = meta.copy()
    out["pattern"] = pattern
    out["up_in"] = [", ".join(c for c in matrix.columns if pd.notna(matrix.loc[i, c])
                              and matrix.loc[i, c] > 0) for i in matrix.index]
    out["down_in"] = [", ".join(c for c in matrix.columns if pd.notna(matrix.loc[i, c])
                                and matrix.loc[i, c] < 0) for i in matrix.index]
    return out


def rank_core(meta: pd.DataFrame, matrix: pd.DataFrame):
    """Rank proteins by recurrence first, then by Pi-score.

    Recurrence is what makes a protein a candidate for the core signature; the
    Pi-score breaks ties by effect size and confidence. The Pi-score is deliberately
    never used to include or exclude a protein - that stays with the adjusted p-value.
    """
    out = meta.copy()
    n_contrasts = matrix.shape[1] or 1
    out["recurrence"] = out["n_sig"] / n_contrasts
    if "best_pi" in out.columns:
        pi = out["best_pi"].fillna(0.0)
        rng = pi.max() - pi.min()
        out["pi_rank"] = (pi - pi.min()) / rng if rng > 0 else 0.0
    else:
        out["pi_rank"] = out["max_abs_lfc"].rank(pct=True).fillna(0.0)
    # recurrence dominates, Pi-score orders within it
    out["core_score"] = out["recurrence"] + 0.25 * out["pi_rank"]
    sort_cols = ["n_sig", "best_pi"] if "best_pi" in out.columns else ["n_sig", "max_abs_lfc"]
    return out.sort_values(sort_cols, ascending=False)


def run(stats, p_thr: float = 0.05):
    matrix, meta, pi_df, rob_df = build_matrix(stats, p_thr)
    if matrix.empty:
        return {"matrix": matrix, "meta": meta, "table": meta,
                "contrast_corr": pd.DataFrame(), "contrast_order": [],
                "contrast_shared": pd.DataFrame(), "contrast_relation": None,
                "has_pi": False, "n_significant": 0}
    direction = directionality(matrix, meta)
    table = rank_core(direction, matrix)
    # same clustering machinery as the GO axis, but on protein profiles
    from axis_meta_go import contrast_similarity
    ccorr, corder, cshared, crel = contrast_similarity(matrix, min_shared=5)
    return {"matrix": matrix, "meta": meta, "table": table,
            "contrast_corr": ccorr, "contrast_order": corder,
            "contrast_shared": cshared, "contrast_relation": crel,
            "has_pi": pi_df is not None, "has_robustness": rob_df is not None,
            "n_significant": int(len(matrix))}
