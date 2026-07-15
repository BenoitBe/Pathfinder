"""
Automatic interpretive summary (deterministic, offline).

Reads the objects already computed (no statistics are recomputed) and produces a
plain-English synthesis: design capture by the PCA, batch suspicion, outlier
samples, functional synthesis per source. Emitted as Markdown, as rows (Excel)
and as HTML.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

# effect-size thresholds (epsilon^2 / eta^2, conventional cut-offs)
EFF_LARGE, EFF_MED = 0.14, 0.06
# BH-adjusted significance threshold, shared by the design and multi-factor sections
# (FDR 10%: PCs are numerous and correlated, and this is exploratory QC, not inference)
P_ADJ_THR = 0.10

SRC_ORDER = ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"]
SRC_LABEL = {"GO:BP": "GO:BP (biological process)", "GO:MF": "GO:MF (molecular function)",
             "GO:CC": "GO:CC (cellular component)", "KEGG": "KEGG", "REAC": "Reactome"}

TECH = ("batch", "lot", "plate", "plaque", "order", "ordre", "injection", "run",
        "date", "operator", "operateur", "column", "colonne")


def _fmt_p(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "NA"
    return f"{p:.1e}" if p < 1e-3 else f"{p:.3f}"


def _sources_present(meta):
    present = [s for s in SRC_ORDER if (meta["source"] == s).any()]
    others = [s for s in meta["source"].dropna().unique() if s not in SRC_ORDER]
    return present + sorted(others)


def _design_section(explore):
    pca = explore["pca"]; assoc = explore["pc_factor"].set_index("PC")
    ve = {f"PC{i+1}": float(v) for i, v in enumerate(pca["var_exp"])}
    scores = pca["scores"]
    sizes = scores["condition"].value_counts()
    k, n, min_n = len(sizes), len(scores), int(sizes.min())
    lines = []

    # no usable replicates -> test not computable, no batch alarm
    if n <= k or min_n < 2 or assoc["p_value"].isna().all():
        lines.append(f"{k} groups for {n} samples (min {min_n}/group): PC ↔ condition test "
                     f"not computable without replicates — inspect the PCA visually "
                     f"(Exploration tab).")
        top = sorted(ve.items(), key=lambda kv: -kv[1])[:2]
        lines.append("Variance concentrated on " +
                     ", ".join(f"{p} ({v:.0f}%)" for p, v in top) +
                     (" — a strongly dominant PC may reflect an atypical sample or a technical gradient."
                      if top and top[0][1] >= 80 else "."))
        return lines

    captured, batch = [], []
    pcol = "p_adj_BH" if "p_adj_BH" in assoc.columns else "p_value"
    for pc in assoc.index:
        eps2, p, var = assoc.loc[pc, "epsilon2"], assoc.loc[pc, pcol], ve.get(pc, 0)
        if pd.notna(p) and p < P_ADJ_THR and eps2 >= EFF_MED:
            captured.append((pc, var, eps2, p))
        elif var >= 10 and (pd.isna(p) or p > 0.2):
            batch.append((pc, var, eps2, p))

    if captured:
        parts = [f"{pc} ({var:.0f}% var., ε²={eps2:.2f}, p={_fmt_p(p)})" for pc, var, eps2, p in captured]
        strong = any(e >= EFF_LARGE for _, _, e, _ in captured)
        lines.append(f"Condition structures the variance on {', '.join(parts)} — "
                     f"group separation is {'clear' if strong else 'moderate'}.")
    else:
        lines.append("No PC is significantly associated with condition: weak group separation, "
                     "or variance dominated by other factors.")
    for pc, var, eps2, p in batch:
        lines.append(f"⚠ {pc} carries {var:.0f}% of the variance with no link to condition "
                     f"(p={_fmt_p(p)}) — suspect a batch effect or a technical covariate.")
    return lines


def _multifactor_section(explore):
    multi = explore.get("pc_factor_multi")
    if multi is None or multi.empty:
        return []
    ve = {f"PC{i+1}": float(v) for i, v in enumerate(explore["pca"]["var_exp"])}
    lines = []
    hits = multi[(multi["p_adj_BH"] < P_ADJ_THR) & (multi["effect"] >= EFF_LARGE)]
    if hits.empty:
        return ["No design factor is significantly associated with a PC (after BH correction)."]
    for fac, grp in hits.groupby("factor"):
        g = grp.sort_values("effect", ascending=False)
        parts = [f"{r.PC} ({ve.get(r.PC, 0):.0f}% var., effect={r.effect:.2f})"
                 for _, r in g.head(3).iterrows()]
        tech = any(t in str(fac).lower() for t in TECH)
        prefix = "⚠ Technical factor" if tech else "Factor"
        tail = (" — technical variance: consider correcting it (limma::removeBatchEffect, ComBat) "
                "or including it as a covariate." if tech else "")
        lines.append(f"{prefix} “{fac}” associated with {', '.join(parts)}{tail}")
    return lines


def _qc_section(explore, cond_map):
    corr = explore["correlation"]
    zscored = bool(corr.attrs.get("zscored"))
    samples = list(corr.columns)
    vals = {}
    for s in samples:
        peers = [o for o in samples if o != s and cond_map.get(o) == cond_map.get(s)]
        if len(peers) >= 2:
            vals[s] = float(np.mean([corr.loc[s, o] for o in peers]))
    if not vals:
        return ["Not enough replicates per condition to assess sample homogeneity."]

    v = np.array(list(vals.values()))
    med = float(np.median(v))
    mad = float(np.median(np.abs(v - med))) * 1.4826
    lines, outliers = [], []
    for s, r in vals.items():
        # robust rule (works on z-scored data, where absolute r is low by construction)
        rel = mad > 1e-9 and r < med - 3 * mad
        absolute = (not zscored) and r < 0.90
        if rel or absolute:
            outliers.append((s, cond_map.get(s), r))
    if outliers:
        for s, c, r in sorted(outliers, key=lambda t: t[2]):
            lines.append(f"⚠ Sample {s} ({c}) correlates poorly with its replicates "
                         f"(mean r={r:.2f} vs median {med:.2f} across samples) — check the "
                         "run/preparation before interpreting.")
    else:
        lines.append(f"Within-condition correlation is homogeneous (median mean r={med:.2f}); "
                     "no outlier sample detected.")
    if zscored:
        lines.append("Correlations are computed on per-protein z-scores (no raw intensity matrix "
                     "available): absolute values are low by construction, only relative "
                     "differences between samples are informative.")
    return lines


def _go_section(meta_go):
    matrix, meta, direction = meta_go["matrix"], meta_go["meta"], meta_go["directionality"]
    if matrix.empty:
        return ["GO enrichment not provided — functional synthesis unavailable."]
    lines = [f"{matrix.shape[0]} significant terms across {matrix.shape[1]} contrast(s), by source:"]
    for src in _sources_present(meta):
        keys = meta.index[meta["source"] == src]
        d = direction[direction["source"] == src] if "source" in direction.columns else direction.loc[keys]
        vc = d["pattern"].value_counts().to_dict()
        shared = vc.get("shared up", 0) + vc.get("shared down", 0)
        mixed = vc.get("mixed", 0)
        rec = meta.loc[keys].sort_values("n_sig", ascending=False)
        dom = (f" — top: {rec.iloc[0]['term_name']} ({int(rec.iloc[0]['n_sig'])} contrasts)"
               if len(rec) else "")
        lines.append(f"{SRC_LABEL.get(src, src)}: {len(keys)} terms "
                     f"({shared} shared, {mixed} with opposite direction){dom}.")
    mx = direction[direction["pattern"] == "mixed"].head(3)
    if len(mx):
        lines.append("Opposite regulation depending on contrast (worth investigating): "
                     + "; ".join(f"{r.term_name} [{r.source}]" for _, r in mx.iterrows()) + ".")
    if (meta["redundancy_cluster"] >= 0).any():
        nc = meta["redundancy_cluster"].nunique()
        lines.append(f"After redundancy grouping: {nc} functional groups.")
    return lines


def _contrast_section(meta_go):
    matrix = meta_go["matrix"]
    if matrix.empty:
        return []
    counts = matrix.notna().sum(axis=0).sort_values(ascending=False)
    return [f"{c}: {int(n)} significantly enriched terms." for c, n in counts.items()]


def generate(stats, meta_go, explore, provenance):
    sections = [
        ("Experimental design (PCA)", _design_section(explore)),
        ("Design factors (PC × factors)", _multifactor_section(explore)),
        ("Quality & homogeneity", _qc_section(explore, stats.cond_map)),
        ("Functional synthesis (Meta-GO)", _go_section(meta_go)),
        ("Per-contrast detail", _contrast_section(meta_go)),
    ]
    sections = [(t, l) for t, l in sections if l]

    md = ["# Automatic summary — Pathfinder",
          f"_{provenance.get('date','')} · {provenance.get('contrasts','')}_", ""]
    lines_flat = []
    for title, lines in sections:
        md.append(f"## {title}")
        for l in lines:
            md.append(f"- {l}")
            lines_flat.append({"section": title, "point": l})
        md.append("")
    return {"markdown": "\n".join(md),
            "sections": sections,
            "lines": pd.DataFrame(lines_flat)}
