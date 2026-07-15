"""
Exploration axis - unsupervised assessment of the experimental design.

- PCA mirroring DEP::plot_pca (SVD on the centred samples x proteins matrix).
- Scree + loadings (proteins driving each PC).
- PC <-> factor association: Kruskal-Wallis of PC scores by condition + epsilon^2
  (effect size); Spearman for numeric factors (e.g. injection order). This is what
  reveals whether an axis captures the design or a batch effect.
- Inter-sample correlation (Pearson) + hierarchical ordering.
- UMAP: guard-rail n_samples < 4 x n_conditions -> sample UMAP skipped in favour of
  the PCA; protein UMAP always available. Optional (skipped if umap-learn absent).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, leaves_list


def pca(comp: pd.DataFrame, sample_int_cols, cond_map, id_labels=None):
    if len(sample_int_cols) < 3:
        raise ValueError(
            f"Only {len(sample_int_cols)} sample column(s) detected: the per-sample matrix "
            "(per-run intensities) could not be found in this workbook. Run --diagnose to locate "
            "the sheet, or provide the intensity file.")
    mat = comp[sample_int_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).values.T  # samples x proteins
    mat_c = mat - mat.mean(axis=0)
    U, S, Vt = np.linalg.svd(mat_c, full_matrices=False)
    scores = U * S                                    # samples x PC
    var_exp = (S ** 2 / np.sum(S ** 2) * 100)
    npc = min(10, scores.shape[1])

    scores_df = pd.DataFrame(
        scores[:, :npc], columns=[f"PC{i+1}" for i in range(npc)],
        index=sample_int_cols,
    )
    scores_df.insert(0, "condition", [cond_map[c] for c in sample_int_cols])

    loadings = Vt.T[:, :npc]                          # proteins x PC
    idx = id_labels if id_labels is not None else np.arange(loadings.shape[0])
    load_df = pd.DataFrame(loadings, columns=[f"PC{i+1}" for i in range(npc)], index=idx)

    return {"scores": scores_df, "loadings": load_df, "var_exp": var_exp[:npc]}


def pc_factor_association(scores_df: pd.DataFrame) -> pd.DataFrame:
    """Kruskal-Wallis of each PC score grouped by condition + epsilon^2 (effect size)."""
    conds = scores_df["condition"].values
    groups_idx = {g: np.where(conds == g)[0] for g in pd.unique(conds)}
    k = len(groups_idx)
    n = len(conds)
    recs = []
    pcs = [c for c in scores_df.columns if c.startswith("PC")]
    for pc in pcs:
        vals = scores_df[pc].values
        samples = [vals[ix] for ix in groups_idx.values() if len(ix) > 0]
        if k >= 2 and all(len(s) >= 1 for s in samples) and n > k:
            try:
                H, p = stats.kruskal(*samples)
            except ValueError:
                H, p = np.nan, np.nan
            eps2 = (H - k + 1) / (n - k) if np.isfinite(H) else np.nan   # epsilon^2 (taille d'effet)
        else:
            H, p, eps2 = np.nan, np.nan, np.nan
        recs.append({"PC": pc, "H": H, "p_value": p, "epsilon2": eps2})
    out = pd.DataFrame(recs)
    # epsilon^2 clipped to [0,1] (negative values mean no association, a formula artefact)
    out["epsilon2"] = out["epsilon2"].clip(lower=0.0, upper=1.0)
    # multiple-testing correction across all tested PCs (Benjamini-Hochberg)
    p = out["p_value"].values.astype(float)
    ok = np.isfinite(p)
    padj = np.full(len(p), np.nan)
    if ok.sum() > 0:
        pv = p[ok]; m = len(pv)
        order = np.argsort(pv)
        ranked = pv[order] * m / (np.arange(m) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        tmp = np.empty(m); tmp[order] = np.clip(ranked, 0, 1)
        padj[ok] = tmp
    out["p_adj_BH"] = padj
    return out[["PC", "H", "p_value", "p_adj_BH", "epsilon2"]]


def _bh(p):
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    out = np.full(len(p), np.nan)
    if ok.sum():
        pv = p[ok]; m = len(pv); order = np.argsort(pv)
        ranked = pv[order] * m / (np.arange(m) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        tmp = np.empty(m); tmp[order] = np.clip(ranked, 0, 1)
        out[ok] = tmp
    return out


def pc_factor_multi(scores_df, factors_df, max_pc=6):
    """Test every PC against every design factor.

    - categorical factor : Kruskal-Wallis + epsilon^2 (effect size)
    - numeric factor     : Spearman correlation (|rho| as effect size)
      -> typically the injection order: a high |rho| on a PC means technical drift.
    BH correction across all tests.
    """
    if factors_df is None or factors_df.empty:
        return pd.DataFrame()
    pcs = [c for c in scores_df.columns if c.startswith("PC")][:max_pc]
    recs = []
    for fac in factors_df.columns:
        vals = factors_df.loc[scores_df.index, fac]
        num = pd.to_numeric(vals, errors="coerce")
        is_num = num.notna().all() and num.nunique() > max(2, 0.5 * len(num))
        for pc in pcs:
            y = scores_df[pc].values.astype(float)
            if is_num:
                try:
                    rho, p = stats.spearmanr(num.values.astype(float), y)
                except Exception:
                    rho, p = np.nan, np.nan
                recs.append({"factor": str(fac), "type": "numeric", "PC": pc,
                             "test": "Spearman", "stat": rho, "p_value": p,
                             "effect": abs(rho) if np.isfinite(rho) else np.nan})
            else:
                groups = [y[np.where(vals.values == g)[0]] for g in pd.unique(vals.values)]
                groups = [g for g in groups if len(g) > 0]
                k, n = len(groups), len(y)
                if k >= 2 and n > k and min(len(g) for g in groups) >= 1:
                    try:
                        H, p = stats.kruskal(*groups)
                    except ValueError:
                        H, p = np.nan, np.nan
                    eps2 = (H - k + 1) / (n - k) if np.isfinite(H) else np.nan
                else:
                    H, p, eps2 = np.nan, np.nan, np.nan
                recs.append({"factor": str(fac), "type": "categorical", "PC": pc,
                             "test": "Kruskal-Wallis", "stat": H, "p_value": p,
                             "effect": min(max(eps2, 0.0), 1.0) if np.isfinite(eps2) else np.nan})
    out = pd.DataFrame(recs)
    if not out.empty:
        out["p_adj_BH"] = _bh(out["p_value"].values)
        out = out[["factor", "type", "PC", "test", "stat", "p_value", "p_adj_BH", "effect"]]
    return out
def sample_correlation(comp: pd.DataFrame, sample_int_cols):
    X = comp[sample_int_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    corr = X.corr(method="pearson")
    # hierarchical ordering
    try:
        d = 1 - corr.values
        np.fill_diagonal(d, 0.0)
        from scipy.spatial.distance import squareform
        Z = linkage(squareform(d, checks=False), method="average")
        order = leaves_list(Z)
        corr = corr.iloc[order, order]
    except Exception:
        pass
    return corr


def umap_samples(comp, sample_int_cols, cond_map):
    """Sample UMAP (design view), guarded: requires n_samples >= 4 x n_conditions."""
    try:
        import umap  # noqa
    except Exception:
        return None
    n = len(sample_int_cols)
    ncond = len(set(cond_map[c] for c in sample_int_cols))
    if n < max(4 * ncond, 10):
        return {"df": None, "warning": f"n={n} for {ncond} conditions: sample UMAP unreliable, PCA preferred."}
    X = comp[sample_int_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).values.T
    emb = umap.UMAP(n_neighbors=min(15, n - 1), min_dist=0.1, random_state=0).fit_transform(X)
    df = pd.DataFrame(emb, columns=["UMAP1", "UMAP2"], index=sample_int_cols)
    df["condition"] = [cond_map[c] for c in sample_int_cols]
    return {"df": df, "warning": None}


def umap_proteins(comp, sample_cols, contrasts, labels, ids, max_pts=6000):
    """Protein UMAP (functional structure) + colouring values.

    Colouring: mean abundance (log2), per-contrast LFC, and (via build_report_html)
    GO membership. `labels`=genes, `ids`=accessions (to match GO intersections).
    """
    try:
        import umap  # noqa
    except Exception:
        return {"df": None, "warning": "umap-learn not installed — protein UMAP skipped."}
    if len(sample_cols) < 3:
        return {"df": None, "warning": "no per-sample intensity matrix available for the protein UMAP."}
    X = comp[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).values  # proteins x samples
    idx = np.arange(X.shape[0])
    if X.shape[0] > max_pts:
        idx = np.sort(np.random.default_rng(0).choice(X.shape[0], max_pts, replace=False))
    emb = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=0).fit_transform(X[idx])
    df = pd.DataFrame(emb, columns=["UMAP1", "UMAP2"])
    df["abondance"] = X[idx].mean(axis=1)
    df["label"] = np.asarray(labels)[idx]
    df["pid"] = np.asarray(ids)[idx]
    lfc_used = []
    for c in contrasts:
        col = f"{c}_diff"
        if col in comp.columns:
            df[f"LFC:{c}"] = pd.to_numeric(comp[col], errors="coerce").values[idx]
            lfc_used.append(c)
    return {"df": df, "warning": None, "contrasts_lfc": lfc_used}


def run(stats_data, do_umap=True, factors_df=None):
    comp = stats_data.pca_df                       # per-sample matrix selected by the loader
    cols = stats_data.sample_int_cols
    cmap = stats_data.cond_map
    gcol = stats_data.gene_col if (stats_data.gene_col in comp.columns) else None
    icol = stats_data.id_col if (stats_data.id_col in comp.columns) else None
    if gcol:
        labels = comp[gcol].astype(str).values
    elif icol:
        labels = comp[icol].astype(str).values
    else:
        labels = np.arange(len(comp)).astype(str)

    pca_res = pca(comp, cols, cmap, id_labels=labels)
    assoc = pc_factor_association(pca_res["scores"])
    # QC correlation: RAW intensities when available (per-protein z-scores are
    # decorrelated by construction); otherwise fall back to the loaded PCA matrix.
    raw_df, raw_cols = stats_data.comp, [c for c in cols if c in stats_data.comp.columns]
    if len(raw_cols) < 3:
        raw_cols = [c for c in (stats_data.comp_sample_cols or []) if c in stats_data.comp.columns]
    if len(raw_cols) < 3:                      # comp holds no per-sample intensities
        raw_df, raw_cols = comp, cols
    corr = sample_correlation(raw_df, raw_cols)

    umap_s = umap_p = None
    if do_umap:
        umap_s = umap_samples(comp, cols, cmap)
        # protein UMAP base: the sheet that actually carries the sample columns
        base, raw = stats_data.comp, [c for c in (stats_data.comp_sample_cols or []) if c in stats_data.comp.columns]
        if len(raw) < 3:
            base, raw = comp, cols                  # fall back to the PCA matrix (e.g. Z-score)
        glab = base[stats_data.gene_col].astype(str).values if (stats_data.gene_col in base.columns) \
            else base[stats_data.id_col].astype(str).values if (stats_data.id_col in base.columns) \
            else np.arange(len(base)).astype(str)
        gids = base[stats_data.id_col].astype(str).values if (stats_data.id_col in base.columns) \
            else glab
        umap_p = umap_proteins(base, raw, stats_data.contrasts, glab, gids)

    return {"pca": pca_res, "pc_factor": assoc, "correlation": corr,
            "pc_factor_multi": pc_factor_multi(pca_res["scores"], factors_df),
            "umap_samples": umap_s, "umap_proteins": umap_p}
