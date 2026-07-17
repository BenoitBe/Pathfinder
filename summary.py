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
# A correlation drop only matters if it is material, not merely statistical: with a
# tight replicate distribution, 3*MAD can flag a 0.02 difference.
QC_MIN_GAP = 0.05

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
        best_p = min(p for _, _, _, p in captured)
        strong = any(e >= EFF_LARGE for _, _, e, _ in captured)
        if best_p < 0.05 and strong:
            verdict = "group separation is clear"
        elif best_p < 0.05:
            verdict = "group separation is moderate"
        else:
            # effect size can be large while the test misses significance: few
            # replicates per group leave little power, especially with many groups.
            verdict = (f"the effect is large but does not reach significance after "
                       f"correction (best p={_fmt_p(best_p)}, {min_n} replicates per group "
                       f"for {k} groups): treat it as a consistent trend, not a proof")
        lines.append(f"Condition structures the variance on {', '.join(parts)} — {verdict}.")
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
        # A tight distribution makes 3*MAD tiny, so the robust rule alone flags noise:
        # require the gap to be material too, and keep the classic absolute rule on
        # raw intensities (replicates below r=0.90 are genuinely worth a look).
        material = (med - r) >= QC_MIN_GAP
        rel = mad > 1e-9 and r < med - 3 * mad and material
        absolute = (not zscored) and r < 0.90 and material
        if rel or absolute:
            outliers.append((s, cond_map.get(s), r))
    if outliers:
        frac = len(outliers) / len(vals)
        if frac > 0.5:
            lines.append(f"⚠ {len(outliers)}/{len(vals)} samples correlate poorly with their "
                         f"replicates (median mean r={med:.2f}) — this is a global pattern, not "
                         "isolated outliers: check normalisation, imputation or the design mapping "
                         "rather than individual runs.")
        else:
            for s, c, r in sorted(outliers, key=lambda t: t[2])[:5]:
                lines.append(f"⚠ Sample {s} ({c}) correlates poorly with its replicates "
                             f"(mean r={r:.2f} vs median {med:.2f} across samples) — check the "
                             "run/preparation before interpreting.")
            if len(outliers) > 5:
                lines.append(f"⚠ {len(outliers) - 5} further sample(s) flagged — see the "
                             "Sample_correlation sheet.")
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
    has_rep = "is_representative" in meta.columns and bool(meta["is_representative"].any())
    n_rep = int(meta["is_representative"].sum()) if has_rep else matrix.shape[0]
    lines = [f"{matrix.shape[0]} significant terms across {matrix.shape[1]} contrast(s)"
             + (f", collapsing to {n_rep} non-redundant terms" if has_rep and n_rep < matrix.shape[0] else "")
             + ", by source:"]
    for src in _sources_present(meta):
        keys = meta.index[meta["source"] == src]
        sub = meta.loc[keys]
        rep_keys = keys[sub["is_representative"]] if has_rep else keys
        d = direction.loc[rep_keys]
        vc = d["pattern"].value_counts().to_dict()
        shared = vc.get("shared up", 0) + vc.get("shared down", 0)
        mixed = vc.get("mixed", 0)
        rec = meta.loc[rep_keys].sort_values("n_sig", ascending=False)
        dom = (f" — top: {rec.iloc[0]['term_name']} ({int(rec.iloc[0]['n_sig'])} contrasts)"
               if len(rec) else "")
        count = f"{len(rep_keys)}" + (f"/{len(keys)}" if has_rep and len(rep_keys) < len(keys) else "")
        lines.append(f"{SRC_LABEL.get(src, src)}: {count} non-redundant terms "
                     f"({shared} shared, {mixed} with opposite direction){dom}.")
    reps = meta.index[meta["is_representative"]] if has_rep else meta.index
    mx = direction.loc[reps]
    mx = mx[mx["pattern"] == "mixed"].head(3)
    if len(mx):
        lines.append("Opposite regulation depending on contrast (worth investigating): "
                     + "; ".join(f"{r.term_name} [{r.source}]" for _, r in mx.iterrows()) + ".")
    if (meta["redundancy_cluster"] >= 0).any():
        nc = meta["redundancy_cluster"].nunique()
        lines.append(f"After redundancy grouping: {nc} functional groups "
                     "(terms hitting the same proteins are collapsed onto one representative).")
    return lines


def _contrast_section(meta_go):
    matrix = meta_go["matrix"]
    if matrix.empty:
        return []
    R_MIN = 0.3          # below this, correlation is not worth a claim
    N_MIN = 10           # a correlation over a handful of terms is noise, not signal
    lines = []
    cc = meta_go.get("contrast_corr")
    rel = meta_go.get("contrast_relation")
    sh = meta_go.get("contrast_shared")
    if cc is not None and not cc.empty and cc.shape[0] >= 2:
        off = cc.where(~np.eye(cc.shape[0], dtype=bool)).stack().dropna()
        if len(off):
            # only pairs sharing no condition are biologically interpretable: a common
            # reference inflates r, a reversed one inflates -r, regardless of biology.
            if rel is not None:
                keep = [rel.loc[a, b] == "disjoint" for a, b in off.index]
                disj = off[keep]
                n_shared_cond = len(off) - len(disj)
            else:
                disj, n_shared_cond = off, 0

            # and only pairs with enough shared terms are worth a headline
            n_thin = 0
            if sh is not None and len(disj):
                enough = [int(sh.loc[a, b]) >= N_MIN for a, b in disj.index]
                n_thin = len(disj) - sum(enough)
                disj = disj[enough]

            def _n(a, b):
                return int(sh.loc[a, b]) if sh is not None else 0

            if len(disj):
                a, b = disj.idxmax()
                if disj.max() >= R_MIN:
                    lines.append(f"Closest independent contrasts: {a} and {b} "
                                 f"(r={disj.max():+.2f} over {_n(a, b)} shared terms) — "
                                 "they move the same functions in the same direction.")
                na, nb = disj.idxmin()
                if disj.min() <= -R_MIN:
                    lines.append(f"Most opposed independent contrasts: {na} and {nb} "
                                 f"(r={disj.min():+.2f} over {_n(na, nb)} shared terms) — "
                                 "shared functions regulated in opposite directions.")
                if disj.max() < R_MIN and disj.min() > -R_MIN:
                    lines.append(f"No clear relationship between independent contrasts "
                                 f"(all |r| < {R_MIN}): their GO signatures are largely unrelated.")
            else:
                lines.append(f"No contrast pair is both independent and based on at least {N_MIN} "
                             "shared terms — no reliable relationship can be reported.")
            notes = []
            if n_shared_cond:
                notes.append(f"{n_shared_cond}/{len(off)} pairs share a condition (a common "
                             "reference mechanically inflates correlation, a reversed one inflates "
                             "anti-correlation)")
            if n_thin:
                notes.append(f"{n_thin} pair(s) rest on fewer than {N_MIN} shared terms")
            if notes:
                lines.append("Excluded from the statements above: " + "; ".join(notes) + ".")
    counts = matrix.notna().sum(axis=0).sort_values(ascending=False)
    lines += [f"{c}: {int(n)} significantly enriched terms." for c, n in counts.items()]
    return lines


def _matrix_caveat(stats):
    if not getattr(stats, "pca_is_zscored", False):
        return []
    return [f"⚠ The PCA runs on '{stats.pca_source}' ({len(stats.pca_df)} proteins), a per-protein "
            "z-scored matrix. If it holds a selected subset (significant or most-variable "
            "proteins), group separation and effect sizes are inflated by construction — the "
            "direction of the result stands, the magnitude does not. Provide the full log2 "
            "intensity matrix for an unbiased assessment."]


def _protein_section(prots):
    """Which proteins carry the signature - the table that ends up in a paper."""
    if prots is None or prots["matrix"].empty:
        return ["No protein reaches significance in any contrast."]
    t, matrix = prots["table"], prots["matrix"]
    n_c = matrix.shape[1]
    lines = [f"{len(matrix)} proteins significant in at least one of the {n_c} contrasts."]
    vc = t["pattern"].value_counts().to_dict()
    shared = vc.get("shared up", 0) + vc.get("shared down", 0)
    specific = vc.get("specific up", 0) + vc.get("specific down", 0)
    opposite = vc.get("opposite", 0)
    lines.append(f"{shared} recurrent (same direction in several contrasts), {specific} "
                 f"specific to one contrast, {opposite} regulated in opposite directions.")

    core = t[t["n_sig"] >= max(2, int(0.5 * n_c))]
    if len(core):
        names = []
        for i, r in core.head(5).iterrows():
            label = r.get("gene") if isinstance(r.get("gene"), str) and r.get("gene") not in ("", "nan") else i
            pi = f", Pi={r['best_pi']:.1f}" if "best_pi" in core.columns and pd.notna(r["best_pi"]) else ""
            names.append(f"{label} ({int(r['n_sig'])}/{n_c} contrasts{pi})")
        lines.append(f"Core signature - {len(core)} proteins significant in at least half the "
                     f"contrasts, top: {'; '.join(names)}.")
    else:
        lines.append("No protein is significant in half the contrasts or more: the signatures "
                     "are largely contrast-specific.")

    opp = t[t["pattern"] == "opposite"]
    if len(opp):
        names = []
        for i, r in opp.head(3).iterrows():
            label = r.get("gene") if isinstance(r.get("gene"), str) and r.get("gene") not in ("", "nan") else i
            names.append(f"{label} (up in {r['up_in'].split(',')[0]}, down in {r['down_in'].split(',')[0]})")
        lines.append(f"Opposite regulation, often the most informative: {'; '.join(names)}.")

    if prots.get("has_pi"):
        lines.append("Proteins are ranked by recurrence then Pi-score (|LFC| x -log10 p). The "
                     "Pi-score orders them; inclusion still rests on the adjusted p-value.")
    return lines


def generate(stats, meta_go, explore, provenance, prots=None):
    sections = [
        ("Experimental design (PCA)", _matrix_caveat(stats) + _design_section(explore)),
        ("Design factors (PC × factors)", _multifactor_section(explore)),
        ("Quality & homogeneity", _qc_section(explore, stats.cond_map)),
        ("Protein-level signature", _protein_section(prots)),
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
