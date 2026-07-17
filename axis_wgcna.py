"""
WGCNA axis - co-expression structure, read from the upstream module.

Deimos already computes the modules; Pathfinder does not recompute them. It answers
the questions the upstream report leaves open:

  - which module responds to which condition (eigengene per sample vs design, the same
    Kruskal-Wallis + epsilon squared used for the PCs);
  - where the core protein signature sits (Fisher exact test per module: are the
    recurrent proteins concentrated somewhere, or spread out?);
  - which proteins are the hubs of each module.

Sheets are detected by CONTENT: 'grey' is WGCNA's bin for unassigned proteins and is
reported but never treated as a module.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd
from scipy import stats as sstats

GREY = ("grey", "grey0", "gray", "ME_grey", "kME_grey")


def _norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def _is_text(series):
    """Text column, whatever pandas calls it: 2.x says 'object', 3.x says 'str'."""
    return not pd.api.types.is_numeric_dtype(series)


def _is_grey(module):
    return _norm(module) in {_norm(g) for g in GREY}


def load_wgcna(xls, sample_cols=None):
    """Pull whatever WGCNA produced out of the workbook. Everything is optional."""
    out = {"assignment": None, "eigengenes": None, "kme": None, "modules": None,
           "module_trait": None, "id_col": None, "module_col": None}
    for name, df in xls.items():
        if df is None or df.empty:
            continue
        low = str(name).lower()
        cols = [str(c) for c in df.columns]

        # eigengenes: one row per sample, columns ME_<module>
        me = [c for c in cols if re.match(r"^me[_\-. ]?\w+", str(c).lower())]
        if me and (sample_cols is None or len(df) == len(sample_cols)):
            eg = df[me].apply(pd.to_numeric, errors="coerce")
            idx = None
            for c in cols:
                if c not in me and _is_text(df[c]):
                    idx = df[c].astype(str)
                    break
            eg.index = idx.values if idx is not None else range(len(eg))
            out["eigengenes"] = eg
            continue

        # kME: protein x module membership
        kme = [c for c in cols if str(c).lower().startswith("kme")]
        if kme:
            k = df[kme].apply(pd.to_numeric, errors="coerce")
            idc = next((c for c in cols if c not in kme and _is_text(df[c])), None)
            if idc is not None:
                k.index = df[idc].astype(str).values
                out["id_col"] = idc
            out["kme"] = k
            continue

        # protein -> module assignment
        mod_col = next((c for c in cols if _norm(c) in ("module", "modules", "color", "colour",
                                                        "modulecolor", "modulecolour")), None)
        if mod_col is not None and len(df) > 20:
            idc = next((c for c in cols
                        if c != mod_col and _is_text(df[c])), None)
            if idc is not None:
                a = pd.DataFrame({"protein": df[idc].astype(str).values,
                                  "module": df[mod_col].astype(str).values})
                hub = next((c for c in cols if "hub" in _norm(c) or _norm(c) == "kme"), None)
                if hub is not None:
                    a["hub_score"] = pd.to_numeric(df[hub], errors="coerce").values
                gene = next((c for c in cols if _norm(c) in ("genes", "gene", "genename")), None)
                if gene is not None:
                    a["gene"] = df[gene].astype(str).values
                out["assignment"] = a
                out["module_col"] = mod_col
            continue

        # module summary / module-trait: kept as-is for the report
        if "module" in low and len(df) <= 50:
            if "trait" in low:
                out["module_trait"] = df
            else:
                out["modules"] = df
    return out


def module_condition_association(eigengenes, cond_map, sample_cols=None):
    """Test each module eigengene against the condition (Kruskal-Wallis + epsilon^2).

    The eigengene is the module's first principal component: one value per sample, so
    it is testable exactly like a PC score.
    """
    if eigengenes is None or eigengenes.empty or not cond_map:
        return pd.DataFrame()
    eg = eigengenes.copy()
    # align eigengene rows with the samples we know the condition of
    if not set(eg.index) & set(cond_map):
        if sample_cols is not None and len(eg) == len(sample_cols):
            eg.index = list(sample_cols)            # positional fallback
        else:
            return pd.DataFrame()
    conds = [cond_map.get(s) for s in eg.index]
    recs = []
    for m in eg.columns:
        y = pd.to_numeric(eg[m], errors="coerce").values
        groups = [y[[i for i, c in enumerate(conds) if c == g]] for g in pd.unique(pd.Series(conds))]
        groups = [g[~np.isnan(g)] for g in groups if len(g)]
        k, n = len(groups), int(sum(len(g) for g in groups))
        if k >= 2 and n > k and min(len(g) for g in groups) >= 1:
            try:
                H, p = sstats.kruskal(*groups)
            except ValueError:
                H, p = np.nan, np.nan
            eps2 = (H - k + 1) / (n - k) if np.isfinite(H) else np.nan
        else:
            H, p, eps2 = np.nan, np.nan, np.nan
        recs.append({"module": str(m), "H": H, "p_value": p,
                     "effect": min(max(eps2, 0.0), 1.0) if np.isfinite(eps2) else np.nan,
                     "is_grey": _is_grey(m)})
    out = pd.DataFrame(recs)
    if not out.empty:
        from axis_explore import _bh
        out["p_adj_BH"] = _bh(out["p_value"].values)
        out = out[["module", "H", "p_value", "p_adj_BH", "effect", "is_grey"]]
        out = out.sort_values("effect", ascending=False)
    return out


def module_contrast_association(eigengenes, cond_map, contrasts, sample_cols=None):
    """Which module separates which contrast.

    The module-vs-condition test says "this module responds to the design"; it does not
    say *which comparison* it follows. Here each module eigengene is compared between
    the two groups of every contrast (Mann-Whitney U), with the rank-biserial
    correlation as effect size: 1 means the two groups are perfectly separated.
    """
    if eigengenes is None or eigengenes.empty or not cond_map or not contrasts:
        return pd.DataFrame()
    eg = eigengenes.copy()
    if not set(eg.index) & set(cond_map):
        if sample_cols is not None and len(eg) == len(sample_cols):
            eg.index = list(sample_cols)
        else:
            return pd.DataFrame()
    conds = np.array([cond_map.get(s) for s in eg.index])

    def _pick(name):
        """Samples of a condition, tolerating spelling differences (CAM203_CTR vs CAM203CTR)."""
        hit = conds == name
        if hit.any():
            return hit
        return np.array([_norm(c) == _norm(name) for c in conds])

    recs = []
    for m in eg.columns:
        y = pd.to_numeric(eg[m], errors="coerce").values
        for c in contrasts:
            if "_vs_" not in str(c):
                continue
            a, b = str(c).split("_vs_", 1)
            ya, yb = y[_pick(a)], y[_pick(b)]
            ya, yb = ya[~np.isnan(ya)], yb[~np.isnan(yb)]
            if len(ya) < 2 or len(yb) < 2:
                continue
            try:
                U, p = sstats.mannwhitneyu(ya, yb, alternative="two-sided")
                rb = 2 * U / (len(ya) * len(yb)) - 1        # rank-biserial correlation
            except ValueError:
                U, p, rb = np.nan, np.nan, np.nan
            recs.append({"module": str(m), "contrast": str(c), "U": U, "p_value": p,
                         "effect": abs(rb) if np.isfinite(rb) else np.nan,
                         "direction": "up" if np.isfinite(rb) and rb > 0 else "down",
                         "is_grey": _is_grey(m)})
    out = pd.DataFrame(recs)
    if not out.empty:
        from axis_explore import _bh
        out["p_adj_BH"] = _bh(out["p_value"].values)
        out = out[["module", "contrast", "U", "p_value", "p_adj_BH", "effect",
                   "direction", "is_grey"]]
        out = out.sort_values(["effect", "p_value"], ascending=[False, True])
        # With few replicates the exact test has a hard floor: 3 vs 3 allows only 20
        # permutations, so p can never drop below 0.1 whatever the separation. Saying
        # so is the difference between "not significant" and "not testable".
        sizes = pd.Series(conds).value_counts()
        n1 = int(sizes.min()) if len(sizes) else 0
        if n1 >= 1:
            from math import comb
            floor = 2 / comb(2 * n1, n1) if comb(2 * n1, n1) else 1.0
            out.attrs["p_floor"] = float(min(1.0, floor))
            out.attrs["n_per_group"] = n1
    return out


def core_enrichment(assignment, core_ids, all_ids=None):
    """Fisher exact test per module: is the core signature concentrated there?

    Answers a question neither the Meta-GO nor the protein axis can: does the recurrent
    signature map onto one co-expression module, or is it spread across the network?
    """
    if assignment is None or assignment.empty or not len(core_ids):
        return pd.DataFrame()
    a = assignment.copy()
    a["_n"] = a["protein"].map(_norm)
    core = {_norm(x) for x in core_ids}
    universe = a["_n"] if all_ids is None else pd.Series(sorted({_norm(x) for x in all_ids}))
    total = len(universe)
    in_core = a["_n"].isin(core)
    recs = []
    for m, grp in a.groupby("module"):
        k = int(in_core[grp.index].sum())          # core proteins in this module
        n_mod = len(grp)
        n_core = int(in_core.sum())
        table = [[k, n_mod - k], [n_core - k, total - n_mod - n_core + k]]
        try:
            odds, p = sstats.fisher_exact(table, alternative="greater")
        except Exception:
            odds, p = np.nan, np.nan
        recs.append({"module": str(m), "n_proteins": n_mod, "n_core": k,
                     "core_fraction": k / n_mod if n_mod else np.nan,
                     "odds_ratio": odds, "p_value": p, "is_grey": _is_grey(m)})
    out = pd.DataFrame(recs)
    if not out.empty:
        from axis_explore import _bh
        out["p_adj_BH"] = _bh(out["p_value"].values)
        out = out.sort_values("p_value")
    return out


def hubs(assignment, kme=None, top=5):
    """Top proteins per module, by hub score when available, else by kME."""
    if assignment is None or assignment.empty:
        return pd.DataFrame()
    a = assignment.copy()
    if "hub_score" not in a.columns and kme is not None and not kme.empty:
        best = {}
        for p in a["protein"]:
            if p in kme.index:
                row = kme.loc[p]
                row = row if isinstance(row, pd.Series) else row.iloc[0]
                best[p] = float(np.nanmax(np.abs(row.values)))
        a["hub_score"] = a["protein"].map(best)
    if "hub_score" not in a.columns:
        return pd.DataFrame()
    out = (a.sort_values("hub_score", ascending=False)
             .groupby("module", group_keys=False).head(top))
    return out.sort_values(["module", "hub_score"], ascending=[True, False])


def run(xls, stats, core_ids=None):
    w = load_wgcna(xls, sample_cols=stats.sample_int_cols)
    if w["assignment"] is None and w["eigengenes"] is None:
        return None
    assoc = module_condition_association(w["eigengenes"], stats.cond_map,
                                         stats.sample_int_cols)
    per_contrast = module_contrast_association(w["eigengenes"], stats.cond_map,
                                               stats.contrasts, stats.sample_int_cols)
    all_ids = stats.comp[stats.id_col].astype(str).tolist() if stats.id_col in stats.comp.columns else None
    enrich = core_enrichment(w["assignment"], core_ids or [], all_ids)
    hub = hubs(w["assignment"], w["kme"])
    n_mod = 0
    if w["assignment"] is not None:
        n_mod = int(sum(1 for m in w["assignment"]["module"].unique() if not _is_grey(m)))
    return {"assignment": w["assignment"], "eigengenes": w["eigengenes"],
            "kme": w["kme"], "modules": w["modules"], "module_trait": w["module_trait"],
            "association": assoc, "contrast_association": per_contrast,
            "core_enrichment": enrich, "hubs": hub, "n_modules": n_mod}
