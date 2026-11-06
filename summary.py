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


DISPERSION_RATIO = 2.5   # within-condition spread vs median condition (null: ~1-2 % false alarms)


def _dispersion_flags(explore):
    """Heterogeneous conditions: replicates spread along the design axes of the PCA.

    The correlation check above works on all proteins and misses a condition whose
    replicates sit at different points of the biological response (asynchronous
    development, partial treatment...): r stays ~0.95 because most proteins do not
    move. Here the spread of each condition around its centroid is measured on the
    PCs that carry the condition effect and compared with the median condition."""
    try:
        sc = explore["pca"]["scores"]
        assoc = explore["pc_factor"].set_index("PC")
    except (KeyError, TypeError):
        return []
    pcol = "p_adj_BH" if "p_adj_BH" in assoc.columns else "p_value"
    pcs = [pc for pc in assoc.index[:6]
           if pd.notna(assoc.loc[pc, pcol]) and assoc.loc[pc, pcol] < P_ADJ_THR
           and assoc.loc[pc, "epsilon2"] >= EFF_MED][:4]
    if len(pcs) < 2:
        pcs = [c for c in ("PC1", "PC2") if c in sc.columns]
    sizes = sc["condition"].value_counts()
    if len(pcs) < 2 or len(sizes) < 3 or sizes.min() < 3:
        return []
    X = sc[pcs].astype(float)
    cent = X.groupby(sc["condition"]).mean()
    spread, far, members = {}, {}, {}
    for cond, g in X.groupby(sc["condition"]):
        d = np.sqrt(((g - cent.loc[cond]) ** 2).sum(axis=1))
        spread[cond] = float(d.mean())
        far[cond] = d.idxmax()
        members[cond] = list(g.index)
    med = float(np.median(list(spread.values())))
    if med <= 0:
        return []
    out = []
    for cond, v in sorted(spread.items(), key=lambda t: -t[1]):
        ratio = v / med
        if ratio < DISPERSION_RATIO:
            break
        # replicates sitting nearer another condition's centroid than their own
        misplaced = []
        for smp in members[cond]:
            dist = np.sqrt(((cent - X.loc[smp]) ** 2).sum(axis=1))
            other = dist.drop(cond)
            if len(other) and other.min() < dist[cond]:
                misplaced.append(f"{smp} lies closer to {other.idxmin()}")
        where = "; ".join(misplaced) if misplaced else f"{far[cond]} is the farthest"
        out.append(f"⚠ Condition {cond}: replicates {ratio:.1f}× more dispersed than the median "
                   f"condition along the design axes ({', '.join(pcs)}) — heterogeneous "
                   f"response (e.g. asynchronous development) or a sample issue; {where}. "
                   "Global correlations cannot show it; expect weaker, less robust calls "
                   "in the contrasts involving this condition.")
    return out


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
    lines.extend(_dispersion_flags(explore))
    if zscored:
        lines.append("Correlations are computed on per-protein z-scores (no raw intensity matrix "
                     "available): absolute values are low by construction, only relative "
                     "differences between samples are informative.")
    return lines


def _go_section(meta_go):
    matrix, meta, direction = meta_go["matrix"], meta_go["meta"], meta_go["directionality"]
    legacy = meta_go.get("pairs") is None
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
        rec = meta.loc[rep_keys].sort_values("n_sig", ascending=False)
        dom = (f" — top: {rec.iloc[0]['term_name']} ({int(rec.iloc[0]['n_sig'])} contrasts)"
               if len(rec) else "")
        count = f"{len(rep_keys)}" + (f"/{len(keys)}" if has_rep and len(rep_keys) < len(keys) else "")
        if legacy:
            dirn = f"{vc.get('mixed', 0)} with opposite direction"
        else:
            dirn = (f"{vc.get('opposite', 0)} reversed between comparable contrasts, "
                    f"{vc.get('mixed', 0)} mixed")
        lines.append(f"{SRC_LABEL.get(src, src)}: {count} non-redundant terms "
                     f"({shared} shared, {dirn}){dom}.")
    reps = meta.index[meta["is_representative"]] if has_rep else meta.index
    mx = direction.loc[reps]
    mx = mx[mx["pattern"] == ("mixed" if legacy else "opposite")].head(3)
    if len(mx):
        lines.append(("Opposite regulation depending on contrast (worth investigating): "
                      if legacy else "Reversed between comparable contrasts (worth investigating): ")
                     + "; ".join(f"{r.term_name} [{r.source}]" for _, r in mx.iterrows()) + ".")
    if not legacy:
        lines.append("'Mixed' terms go up in some contrasts and down in others that compare "
                     "different conditions: an ordering of the conditions, not a reversal.")
    if (meta["redundancy_cluster"] >= 0).any():
        nc = meta["redundancy_cluster"].nunique()
        lines.append(f"After redundancy grouping: {nc} functional groups "
                     "(terms hitting the same proteins are collapsed onto one representative).")
    return lines


MIN_CONTRASTS_FOR_FAMILIES = 8   # above this, per-contrast lines are grouped into families
FAMILY_R = 0.5                   # GO-signature correlation within a family (average linkage)


def _contrast_families(cc, contrasts):
    """Groups of contrasts with near-identical GO signatures: average-linkage clustering
    on 1 - r (pairs without enough shared terms count as unrelated), cut at r = FAMILY_R."""
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    cc = cc.reindex(index=contrasts, columns=contrasts)
    d = 1.0 - cc.values.astype(float)
    d = np.where(np.isfinite(d), d, 1.0)
    np.fill_diagonal(d, 0.0)
    d = np.clip((d + d.T) / 2.0, 0.0, 2.0)
    labels = fcluster(linkage(squareform(d, checks=False), method="average"),
                      t=1.0 - FAMILY_R, criterion="distance")
    fams = {}
    for c, l in zip(contrasts, labels):
        fams.setdefault(l, []).append(c)
    return list(fams.values())


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
    literal = [f"{c}: {int(n)} significantly enriched terms." for c, n in counts.items()]
    if len(counts) <= MIN_CONTRASTS_FOR_FAMILIES or cc is None or cc.empty:
        return lines + literal
    fams = _contrast_families(cc, list(counts.index))
    if not any(len(f) > 1 for f in fams):
        return lines + literal
    fams.sort(key=lambda f: (-len(f), -float(counts[f].mean())))
    lines.append(f"{len(fams)} families summarising {len(counts)} contrasts (GO signatures "
                 f"correlated at r ≥ {FAMILY_R} within a family; a family may gather contrasts "
                 "sharing a condition, whose similarity is partly mechanical):")
    k = 0
    for fam in fams:
        if len(fam) == 1:
            c = fam[0]
            lines.append(f"{c}: {int(counts[c])} significantly enriched terms.")
            continue
        k += 1
        sub = cc.reindex(index=fam, columns=fam).values.astype(float)
        r = np.nanmean(sub[~np.eye(len(fam), dtype=bool)])
        n = counts[fam]
        lines.append(f"Family {k} — {len(fam)} contrasts (mean r={r:+.2f}, "
                     f"{int(n.min())}–{int(n.max())} terms each): {', '.join(fam)}.")
    return lines


def _matrix_caveat(stats):
    if not getattr(stats, "pca_is_zscored", False):
        return []
    return [f"⚠ The PCA runs on '{stats.pca_source}' ({len(stats.pca_df)} proteins), a per-protein "
            "z-scored matrix. If it holds a selected subset (significant or most-variable "
            "proteins), group separation and effect sizes are inflated by construction — the "
            "direction of the result stands, the magnitude does not. Provide the full log2 "
            "intensity matrix for an unbiased assessment."]


def _significance_rule(provenance):
    """The rule behind 'significant', stated once. By default it is the upstream volcano
    rule, so Pathfinder's proteins are the ones of the Deimos volcanos and dashboard;
    when it was overridden the two counts differ and the summary says so."""
    if not provenance or not provenance.get("significance_rule"):
        return []
    line = f"Significance rule: {provenance['significance_rule']}."
    same = provenance.get("significance_rule_matches_upstream")
    if same == "yes":
        line += " Same proteins as the upstream volcanos and dashboard."
    elif provenance.get("volcano_rule_upstream"):
        line += (f" Upstream volcanos and dashboard use {provenance['volcano_rule_upstream']}: "
                 "the two counts differ.")
    return [line]


def _protein_section(prots, provenance=None):
    """Which proteins carry the signature - the table that ends up in a paper."""
    rule = _significance_rule(provenance)
    if prots is None or prots["matrix"].empty:
        return ["No protein reaches significance in any contrast."] + rule
    t, matrix = prots["table"], prots["matrix"]
    n_c = matrix.shape[1]
    lines = [f"{len(matrix)} proteins significant in at least one of the {n_c} contrasts."] + rule
    vc = t["pattern"].value_counts().to_dict()
    shared = vc.get("shared up", 0) + vc.get("shared down", 0)
    specific = vc.get("specific up", 0) + vc.get("specific down", 0)
    opposite = vc.get("opposite", 0)
    mixed = vc.get("mixed", 0)
    if prots.get("pairs") is None:
        lines.append(f"{shared} recurrent (same direction in several contrasts), {specific} "
                     f"specific to one contrast, {opposite} regulated in opposite directions.")
    else:
        lines.append(f"{shared} recurrent (same direction in several contrasts), {specific} "
                     f"specific to one contrast, {opposite} reversed between comparable contrasts "
                     f"(same comparison in different strata), {mixed} up in some contrasts and "
                     "down in others that compare different conditions (an ordering, not a "
                     "contradiction).")

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
            where = (str(r.get("opposite_in", "")).split(";")[0] if r.get("opposite_in")
                     else f"up in {r['up_in'].split(',')[0]}, down in {r['down_in'].split(',')[0]}")
            names.append(f"{label} ({where})")
        lines.append(f"Reversed regulation, often the most informative: {'; '.join(names)}.")
    elif prots.get("pairs") is not None:
        lines.append("No protein is reversed between comparable contrasts"
                     + (" (no pair of contrasts measures the same comparison in different "
                        "strata)." if not prots["pairs"] else "."))

    if prots.get("has_pi"):
        lines.append("Proteins are ranked by recurrence then Pi-score (|LFC| x -log10 p). The "
                     "Pi-score orders them; inclusion still rests on the adjusted p-value.")
    return lines


def _wgcna_section(wgcna):
    """Co-expression structure: which module tracks the design, and where the core sits."""
    if not wgcna or not wgcna.get("n_modules"):
        return []
    lines = [f"{wgcna['n_modules']} co-expression modules read from the upstream WGCNA "
             "(not recomputed)."]
    a = wgcna.get("association")
    if a is not None and not a.empty:
        hits = a[(~a["is_grey"]) & (a["p_adj_BH"] < P_ADJ_THR) & (a["effect"] >= EFF_LARGE)]
        if len(hits):
            parts = [f"{r['module']} (effect={r['effect']:.2f}, p={_fmt_p(r['p_adj_BH'])})"
                     for _, r in hits.head(3).iterrows()]
            lines.append(f"Module eigengenes tracking the condition: {', '.join(parts)} — "
                         "these modules respond to the design.")
        else:
            lines.append("No module eigengene is significantly associated with the condition "
                         "after correction: co-expression structure and design do not align.")
    ca = wgcna.get("contrast_association")
    if ca is not None and not ca.empty:
        best = ca[~ca["is_grey"]].head(1)
        if len(best):
            r = best.iloc[0]
            floor = ca.attrs.get("p_floor")
            note = ""
            if floor and floor > 0.05:
                note = (f" (with {ca.attrs.get('n_per_group')} replicates per group the exact "
                        f"test cannot go below p={floor:.2f}: read the separation, not the p)")
            lines.append(f"Module {r['module']} separates {r['contrast']} best "
                         f"(separation={r['effect']:.2f}, p={_fmt_p(r['p_value'])}){note}.")
    e = wgcna.get("core_enrichment")
    if e is not None and not e.empty:
        hits = e[(~e["is_grey"]) & (e["p_adj_BH"] < 0.05)]
        if len(hits):
            r = hits.iloc[0]
            rule = wgcna.get("core_set_rule", "the core protein signature")
            lines.append(f"The top-ranked proteins ({rule}) concentrate in module "
                         f"{r['module']} ({int(r['n_core'])}/{int(r['n_proteins'])} of its "
                         f"proteins, p={_fmt_p(r['p_adj_BH'])}) — they are co-expressed, not "
                         "scattered across the network.")
        else:
            lines.append("The core signature is not concentrated in any module: recurrent "
                         "proteins are spread across the co-expression network.")
    return lines


def _headline(stats, meta_go, explore, prots, wgcna):
    """The five facts a reader should get at a glance, as structured data.

    level: 'ok' (nothing to flag), 'warn' (needs attention), 'info' (neutral).
    """
    cards = []

    # --- design
    assoc = explore["pc_factor"].set_index("PC")
    ve = {f"PC{i+1}": float(v) for i, v in enumerate(explore["pca"]["var_exp"])}
    pcol = "p_adj_BH" if "p_adj_BH" in assoc.columns else "p_value"
    good = [(pc, ve.get(pc, 0), assoc.loc[pc, "epsilon2"], assoc.loc[pc, pcol])
            for pc in assoc.index
            if pd.notna(assoc.loc[pc, pcol]) and assoc.loc[pc, pcol] < P_ADJ_THR
            and assoc.loc[pc, "epsilon2"] >= EFF_MED]
    if assoc[pcol].isna().all():
        cards.append({"label": "Design", "value": "not testable", "level": "warn",
                      "detail": "no replicates"})
    elif good:
        pc, var, eps2, p = max(good, key=lambda t: t[1])
        sig = p < 0.05
        cards.append({"label": "Design", "value": f"{pc} · {var:.0f}% var.",
                      "detail": f"ε²={eps2:.2f}, p={_fmt_p(p)}" + ("" if sig else " (trend)"),
                      "level": "ok" if sig else "warn"})
    else:
        cards.append({"label": "Design", "value": "not captured", "level": "warn",
                      "detail": "no PC follows the condition"})

    # --- quality
    qc = _qc_section(explore, stats.cond_map)
    flagged = [l for l in qc if l.startswith("⚠")]
    cards.append({"label": "Quality",
                  "value": "clean" if not flagged else f"{len(flagged)} flag(s)",
                  "detail": "replicates homogeneous" if not flagged else flagged[0][:60],
                  "level": "ok" if not flagged else "warn"})

    # --- core proteins
    if prots and not prots["matrix"].empty:
        t, n_c = prots["table"], prots["matrix"].shape[1]
        core = t[t["n_sig"] >= max(2, int(0.5 * n_c))]
        cards.append({"label": "Core proteins", "value": f"{len(core)}",
                      "detail": f"significant in ≥ {max(2, int(0.5*n_c))}/{n_c} contrasts",
                      "level": "info"})
    else:
        cards.append({"label": "Core proteins", "value": "0", "level": "warn",
                      "detail": "nothing significant"})

    # --- module holding the core signature
    mod_card = {"label": "Core module", "value": "—", "level": "info",
                "detail": "no WGCNA module in the workbook"}
    if wgcna and wgcna.get("n_modules"):
        e = wgcna.get("core_enrichment")
        hits = e[(~e["is_grey"]) & (e["p_adj_BH"] < 0.05)] if e is not None and not e.empty \
            else pd.DataFrame()
        if len(hits):
            r = hits.iloc[0]
            n_set = wgcna.get("core_set_n")
            mod_card = {"label": "Core module", "value": str(r["module"]),
                        "detail": (f"{int(r['n_core'])} of the top {n_set} · {int(r['n_proteins'])} "
                                   f"proteins, p={_fmt_p(r['p_adj_BH'])}" if n_set else
                                   f"{int(r['n_core'])}/{int(r['n_proteins'])} proteins, "
                                   f"p={_fmt_p(r['p_adj_BH'])}"), "level": "ok"}
        else:
            mod_card = {"label": "Core module", "value": "scattered", "level": "info",
                        "detail": "top-ranked proteins span several modules"}
    cards.append(mod_card)

    # --- dominant function
    if meta_go and not meta_go["matrix"].empty:
        meta = meta_go["meta"]
        rep = meta[meta["is_representative"]] if "is_representative" in meta.columns else meta
        top = rep.sort_values("n_sig", ascending=False).head(1)
        if len(top):
            r = top.iloc[0]
            cards.append({"label": "Top function", "value": str(r["term_name"])[:40],
                          "detail": f"{int(r['n_sig'])} contrasts · {r['source']}",
                          "level": "info"})
    return cards


def _convergence(prots, wgcna):
    """The chain a reader is really after: core proteins -> module -> comparison.

    Returns None when a link is missing, rather than implying one.
    """
    if not (prots and not prots["matrix"].empty and wgcna and wgcna.get("n_modules")):
        return None
    e = wgcna.get("core_enrichment")
    if e is None or e.empty:
        return None
    hits = e[(~e["is_grey"]) & (e["p_adj_BH"] < 0.05)]
    n_c = prots["matrix"].shape[1]
    # the set actually tested for module enrichment (top-ranked proteins), not the
    # stricter "significant in half the contrasts" count of the headline card
    core_n = int(wgcna.get("core_set_n") or
                 (prots["table"]["n_sig"] >= max(2, int(0.5 * n_c))).sum())
    rule = wgcna.get("core_set_rule", "recurrent proteins")
    if not len(hits):
        return {"core_n": core_n, "core_title": "top-ranked signature",
                "core_detail": "recurrence, then Pi-score", "module": None, "contrast": None,
                "verdict": f"The top-ranked proteins ({rule}) are spread across the "
                           "co-expression network: no module concentrates them."}
    r = hits.iloc[0]
    out = {"core_n": core_n, "core_title": "top-ranked signature",
           "core_detail": "recurrence, then Pi-score", "module": str(r["module"]),
           "module_detail": f"{int(r['n_core'])}/{int(r['n_proteins'])} proteins "
                            f"(p={_fmt_p(r['p_adj_BH'])})",
           "contrast": None, "contrast_detail": None, "verdict": None}

    ca = wgcna.get("contrast_association")
    if ca is not None and not ca.empty:
        # eigengene column is 'ME_<module>' while the assignment says '<module>'
        m = ca[ca["module"].str.replace("^ME[_.-]?", "", regex=True).str.lower()
               == out["module"].lower()]
        if len(m):
            best = m.sort_values(["effect", "p_value"], ascending=[False, True]).iloc[0]
            floor = ca.attrs.get("p_floor")
            note = ""
            if floor and floor > 0.05:
                note = (f" — with {ca.attrs.get('n_per_group')} replicates per group the exact "
                        f"test cannot go below p={floor:.2f}, so read the effect, not the p-value")
            out["contrast"] = str(best["contrast"])
            out["contrast_detail"] = (f"separation={best['effect']:.2f}, "
                                      f"p={_fmt_p(best['p_value'])}{note}")
            out["verdict"] = (f"The top-ranked proteins ({rule}) concentrate in module "
                              f"{out['module']} ({int(r['n_core'])} of them), whose eigengene "
                              f"separates {best['contrast']} best.")
    if not out["verdict"]:
        out["verdict"] = (f"The top-ranked proteins ({rule}) concentrate in module "
                          f"{out['module']}, but no eigengene was available to link it "
                          "to a comparison.")
    return out


def generate(stats, meta_go, explore, provenance, prots=None, wgcna=None):
    sections = [
        ("Experimental design (PCA)", _matrix_caveat(stats) + _design_section(explore)),
        ("Design factors (PC × factors)", _multifactor_section(explore)),
        ("Quality & homogeneity", _qc_section(explore, stats.cond_map)),
        ("Protein-level signature", _protein_section(prots, provenance)),
        ("Co-expression modules (WGCNA)", _wgcna_section(wgcna)),
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
            "headline": _headline(stats, meta_go, explore, prots, wgcna),
            "convergence": _convergence(prots, wgcna),
            "lines": pd.DataFrame(lines_flat)}
