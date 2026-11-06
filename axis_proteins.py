"""
Protein axis - cross-contrast comparison at the protein level.

The Meta-GO axis answers "which functions do these contrasts share"; this one answers
"which proteins carry the signature" - the table that ends up in a paper.

Core object: a `protein x contrast` matrix of log2 fold-changes, masked by
significance (NaN = not significant, never 0).

Selection follows the upstream rule (io_layer.SelectionRule: the Deimos volcano rule
read from the manifest), so the proteins are exactly those Deimos reports. Ranking uses
the Pi-score (|LFC| x -log10 p, computed upstream): it is a ranking statistic, not a
testable one, so it orders proteins but never selects them. Robustness, when the
upstream module provides it, is carried through as a confidence column.
"""
from __future__ import annotations

import re

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


_P_KINDS = {"p.adj": ["p.adj", "padj", "p_adj"], "p.val": ["p.val", "pvalue", "p_val", "p.value"]}


def build_matrix(stats, p_thr: float = 0.05, rule=None):
    """protein x contrast matrix of significant log2 fold-changes + per-protein meta.
    rule: io_layer.SelectionRule; None = p.adj < p_thr without fold-change criterion."""
    from io_layer import SelectionRule
    rule = rule or SelectionRule(p_thr=p_thr)
    comp = stats.comp
    idx = comp[stats.id_col].astype(str) if stats.id_col in comp.columns \
        else pd.Series([f"row{i}" for i in range(len(comp))])

    lfc, pi, rob = {}, {}, {}
    for c in stats.contrasts:
        p_col = _find_contrast_col(comp, c, _P_KINDS.get(rule.p_column, _P_KINDS["p.adj"]))
        d_col = _find_contrast_col(comp, c, ["diff", "logfc", "lfc"])
        if p_col is None or d_col is None:
            continue
        p = pd.to_numeric(comp[p_col], errors="coerce")
        d = pd.to_numeric(comp[d_col], errors="coerce")
        sel = rule.mask(p, d)
        lfc[c] = np.where(sel, d, np.nan)
        pi_col = _find_contrast_col(comp, c, ["pi_score", "piscore", "pi"])
        if pi_col is not None:
            pi[c] = np.where(sel, pd.to_numeric(comp[pi_col], errors="coerce"), np.nan)
        rob_col = _find_contrast_col(comp, c, ["robustness_score", "robustness", "robust"])
        if rob_col is not None:
            rob[c] = np.where(sel, pd.to_numeric(comp[rob_col], errors="coerce"), np.nan)

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


# |r| between the log2 fold-change profiles of two contrasts (all proteins) above which
# they are taken to measure the same comparison in different strata (e.g. 24 h vs 0 h in
# two cultivars). Below it, the two contrasts compare different things and opposite
# signs are not a contradiction.
PARALLEL_R_MIN = 0.3


def lfc_table(stats) -> pd.DataFrame:
    """Unmasked log2 fold-changes, one column per contrast (all proteins)."""
    comp = stats.comp
    out = {}
    for c in stats.contrasts:
        d_col = _find_contrast_col(comp, c, ["diff", "logfc", "lfc"])
        if d_col is not None:
            out[c] = pd.to_numeric(comp[d_col], errors="coerce").values
    return pd.DataFrame(out)


def comparable_pairs(stats, r_min: float = PARALLEL_R_MIN, min_n: int = 50):
    """Pairs of contrasts whose directions can be compared protein by protein.

    Two contrasts are comparable when they share no condition AND their fold-change
    profiles over all proteins are correlated (|r| >= r_min): they then measure the
    same comparison in different strata, and a protein going up in one and down in the
    other is a genuine reversal. Contrasts sharing a condition are never comparable:
    A_vs_B up and A_vs_C down only says B < A < C, an ordering, not a contradiction;
    nor are uncorrelated contrasts (they compare different things).

    Returns {(a, b): +1 | -1} (sign of r: -1 when one contrast is written the other way
    round), {} when nothing is comparable, None when fold-changes are unavailable.
    """
    from axis_meta_go import _pair_relation
    L = lfc_table(stats)
    if L.shape[1] < 2:
        return None
    pairs = {}
    cols = list(L.columns)
    for i, a in enumerate(cols):
        for b in cols[i + 1:]:
            if _pair_relation(a, b)[0] != "disjoint":
                continue
            both = L[[a, b]].dropna()
            if len(both) < min_n or both[a].std() == 0 or both[b].std() == 0:
                continue
            r = float(both[a].corr(both[b]))
            if abs(r) >= r_min:
                pairs[(a, b)] = 1 if r > 0 else -1
    return pairs


def _reversals(values: dict, pairs: dict) -> list:
    """Comparable pairs in which this protein/term moves against the pair's trend."""
    out = []
    for (a, b), sgn in pairs.items():
        va, vb = values.get(a), values.get(b)
        if va is None or vb is None or not (np.isfinite(va) and np.isfinite(vb)):
            continue
        if np.sign(va) * np.sign(vb) == -sgn:
            out.append(f"{a} ({va:+.2f}) vs {b} ({vb:+.2f})")
    return out


def directionality(matrix: pd.DataFrame, meta: pd.DataFrame, pairs=None):
    """Per protein: shared / specific / opposite / mixed regulation across contrasts.

    pairs (see comparable_pairs): 'opposite' = reversed between two comparable
    contrasts; 'mixed' = up in some contrasts and down in others that compare
    different things (not a contradiction). Without pairs (fold-changes unavailable),
    any protein both up and down is called 'opposite', as before."""
    up = meta["n_up"].values
    dn = meta["n_down"].values
    out = meta.copy()
    if pairs is None:
        both = np.where((up > 0) & (dn > 0), "opposite", "")
        out["opposite_in"] = ""
    else:
        rev = []
        for i in matrix.index:
            row = matrix.loc[i]
            vals = {c: float(v) for c, v in row.items() if pd.notna(v)}
            rev.append("; ".join(_reversals(vals, pairs)) if (len(vals) > 1) else "")
        out["opposite_in"] = rev
        has_rev = np.array([bool(r) for r in rev])
        both = np.where(has_rev, "opposite", np.where((up > 0) & (dn > 0), "mixed", ""))
    pattern = np.where(both != "", both,
                       np.where(up > 1, "shared up",
                                np.where(dn > 1, "shared down",
                                         np.where(up == 1, "specific up", "specific down"))))
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


def _keys(pid, gene) -> set:
    """Identifiers a protein can be matched on: its accession and its gene symbol(s)
    (gProfiler / STRING may return symbols rather than accessions)."""
    out = {str(pid).strip().upper()}
    g = "" if gene is None else str(gene)
    if g and g.lower() not in ("nan", "none"):
        out |= {x.strip().upper() for x in re.split(r"[;,]", g) if x.strip()}
    return out


def annotate(table: pd.DataFrame, go_res=None, wgcna=None, max_terms: int = 8):
    """The join no single tab provides: which proteins, doing what, in which module.

    go_res (axis_meta_go.run): adds `go_terms` (names of the significant terms whose
    genes include the protein, most recurrent first, at most `max_terms` listed) and
    `n_go_terms` (all of them). wgcna (axis_wgcna.run): adds `module`. Proteins are
    matched on accession or gene symbol. Without GO nor WGCNA the table is returned
    unchanged."""
    out = table.copy()
    meta = (go_res or {}).get("meta")
    term_genes = (go_res or {}).get("term_genes") or {}
    has_go = bool(term_genes) and meta is not None and len(meta) > 0
    assign = (wgcna or {}).get("assignment")
    has_wg = assign is not None and not assign.empty
    if not (has_go or has_wg):
        return out
    genes = out["gene"].tolist() if "gene" in out.columns else [None] * len(out)
    keys = [_keys(i, g) for i, g in zip(out.index, genes)]
    if has_go:
        inv = {}
        for tkey, accs in term_genes.items():
            for acc in accs:
                inv.setdefault(str(acc).strip().upper(), set()).add(tkey)
        n_sig = meta["n_sig"] if "n_sig" in meta.columns else pd.Series(0, index=meta.index)
        name = meta["term_name"] if "term_name" in meta.columns else pd.Series(meta.index, index=meta.index)

        def _order(t):
            return (-int(n_sig.get(t, 0) or 0), str(name.get(t, t)))
        labels, counts = [], []
        for ks in keys:
            terms = set().union(*(inv.get(k, set()) for k in ks))
            names = list(dict.fromkeys(str(name.get(t, t)) for t in sorted(terms, key=_order)))
            counts.append(len(terms))
            labels.append("; ".join(names[:max_terms]) + ("; …" if len(names) > max_terms else ""))
        out["n_go_terms"] = counts
        out["go_terms"] = labels
    if has_wg:
        mod = {str(p).strip().upper(): str(m)
               for p, m in zip(assign["protein"], assign["module"])}
        out["module"] = [next((mod[k] for k in ks if k in mod), "") for ks in keys]
    return out


def run(stats, p_thr: float = 0.05, pairs="auto", rule=None):
    """pairs: 'auto' (computed from the fold-changes), a dict from comparable_pairs,
    or None (no comparability information: legacy 'opposite' rule).
    rule: io_layer.SelectionRule (upstream volcano rule); None = p.adj < p_thr."""
    if isinstance(pairs, str) and pairs == "auto":
        pairs = comparable_pairs(stats)
    matrix, meta, pi_df, rob_df = build_matrix(stats, p_thr, rule=rule)
    if matrix.empty:
        return {"matrix": matrix, "meta": meta, "table": meta,
                "contrast_corr": pd.DataFrame(), "contrast_order": [],
                "contrast_shared": pd.DataFrame(), "contrast_relation": None,
                "has_pi": False, "n_significant": 0, "pairs": pairs}
    direction = directionality(matrix, meta, pairs)
    table = rank_core(direction, matrix)
    # same clustering machinery as the GO axis, but on protein profiles
    from axis_meta_go import contrast_similarity
    ccorr, corder, cshared, crel = contrast_similarity(matrix, min_shared=5)
    return {"matrix": matrix, "meta": meta, "table": table,
            "contrast_corr": ccorr, "contrast_order": corder,
            "contrast_shared": cshared, "contrast_relation": crel,
            "has_pi": pi_df is not None, "has_robustness": rob_df is not None,
            "n_significant": int(len(matrix)), "pairs": pairs}
