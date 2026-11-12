"""Orchestration of the biological reading: response programmes, functional groups
(Meta-GO), SVG figures and numeric facts. Every step is guarded: a failing figure is
reported and stops neither the others nor the rest of the report."""
from __future__ import annotations

import html
import os
import traceback

from .factors import detect_design
from .groups import build_groups
from .io import id_lookup, load_run
from .programmes import THEME_Q_THR, build_programmes, programme_stability, stability_class
from .svg import (fig_factor_effects, fig_factor_profiles, fig_factors, fig_groups, fig_landscape,
                  fig_bump_chart, fig_dynamic_range, fig_fc_vs_fc,
                  fig_go_treemap, fig_map, fig_programme_groups, fig_programmes,
                  fig_radar, fig_slope_graph, fig_spotlight, fig_term_tree,
                  fig_variance_violins, fig_volcano, fc_pairs)
from .tables import reading_tables


def _safe(label, fn, *a, **kw):
    try:
        return fn(*a, **kw)
    except Exception as e:                                            # noqa: BLE001
        print(f"[Pathfinder] ⚠ reading: {label} skipped ({type(e).__name__}: {e})")
        if os.environ.get("PATHFINDER_DEBUG"):
            traceback.print_exc()
        return None


def _names(lst, n=5):
    return ", ".join(f"<em>{html.escape(str(x))}</em>" for x in lst[:n]) or "—"


def _gsea_groups(stats_path, run, manifest, lookup, jaccard, forced=False):
    """Functional groups from the GSEA leading edges (ORA empty, or --groups-from gsea):
    pathways grouped like the Meta-GO terms (average-linkage Jaccard, --jaccard cut)."""
    from .gsea import load_gsea
    import axis_meta_go
    gl = _safe("GSEA sheets", load_gsea, stats_path, run.contrasts, manifest, run.de, lookup)
    if gl is None or gl.empty:
        if forced:
            print("[Pathfinder] \u26a0 reading: --groups-from gsea but no usable GSEA_* sheet.")
        return None, ""
    res = _safe("GSEA grouping", axis_meta_go.run, gl, p_thr=1.01, jaccard_cut=jaccard)
    if res is None or res["matrix"].empty:
        return None, ""
    g = _safe("GSEA groups", build_groups, res, gl, lookup, 1.01, "GSEA")
    if g is None:
        return None, ""
    why = "requested (--groups-from gsea)" if forced else "the ORA gave no usable term"
    print(f"[Pathfinder] reading: functional groups from the GSEA leading edges ({why}): "
          f"{gl.term_name.nunique()} pathways -> {len(g.names)} groups")
    return g, why


def _design_factors(run, design, factors_df, factor_regex):
    """Factors of a factorial design (or None), with roles and comparison classes."""
    from .factorial import (condition_order, contrast_factors, design_factors, factor_blocks,
                            roles)
    conds = [c for c in (design.conditions or run.conditions)]
    try:
        F, src = design_factors(run, conds, factors_df, factor_regex)
    except ValueError as e:                    # bad --factor-regex: say so, read automatically
        print(f"[Pathfinder] \u26a0 reading: {e}; factors read automatically instead.")
        F, src = design_factors(run, conds, factors_df)
    if F is None:
        return None
    r = roles(F)
    cf = contrast_factors(run.pairs, F, run.contrasts)
    order = condition_order(F, r)
    blk = {c: F.at[c, r["lines"]] for c in F.index} if r["lines"] else None
    print(f"[Pathfinder] reading: factorial design, factors "
          + ", ".join(f"{f} ({F[f].nunique()})" for f in F.columns)
          + f" [{src}]; x = {r['x']}, lines = {r['lines']}, panels = {', '.join(r['facets']) or '—'}")
    return {"F": F, "source": src, "roles": r, "cf": cf, "blocks": factor_blocks(cf, F),
            "order": order, "blockmap": blk}


def _factorial(run, prog, fd):
    """Variance partition when the design has >= 2 factors (Hordeum-like)."""
    from .factorial import summarise, variance_partition
    if not fd:
        return None
    F, src = fd["F"], fd["source"]
    sig = sum(run.sig(c).fillna(False) for c in run.contrasts).gt(0)
    ids = list(prog.z.index) if prog is not None else list(run.de.index[sig])
    S, terms = variance_partition(run, ids, F)
    if S is None or S.empty:
        return None
    summ = summarise(S, prog)
    print(f"[Pathfinder] reading: variance partition of {len(S)} proteins")
    allr = summ.iloc[0]
    parts = [t for t in terms + ["higher-order", "replicates"] if t in summ.columns]
    top = [t for t in sorted(terms, key=lambda t: -float(allr[t])) if float(allr[t]) >= 0.01]
    from .factorial import MISSING, levels as _levels

    def _desc(c):
        lv_ = [v for v in _levels(F, c) if v != MISSING]
        miss = int((F[c] == MISSING).sum())
        return (f"<b>{html.escape(c)}</b> ({len(lv_)} levels: {html.escape(', '.join(lv_[:9]))}"
                + ("…" if len(lv_) > 9 else "")
                + (f"; absent from {miss} condition(s)" if miss else "") + ")")
    b = [f"Design factors (from {src}): " + "; ".join(_desc(c) for c in F.columns) + "."
         + (" Name them with columns of the design file (e.g. genotype, light, time) or "
            "--factor-regex." if src == "condition names" else ""),
         f"All {int(allr['n_proteins'])} differential proteins, median share of their variance: "
         + (", ".join(f"{html.escape(t)} {float(allr[t]) * 100:.0f} %" for t in top) or "no term ≥ 1 %")
         + f"; higher-order condition effects {float(allr['higher-order']) * 100:.0f} %, "
           f"replicates (noise) {float(allr['replicates']) * 100:.0f} %."]
    for _, d in summ.iloc[1:].iterrows():       # column names hold spaces / '×'
        best = [t for t in sorted(terms, key=lambda t: -float(d[t]))[:2] if float(d[t]) >= 0.05]
        b.append(f"<b>{d['set']}</b>: " + (("mostly " + " and ".join(
            f"{html.escape(t)} ({float(d[t]) * 100:.0f} %)" for t in best)) if best else
            "no design factor carries 5 % of the variance")
            + f" · replicates {float(d['replicates']) * 100:.0f} %.")
    if float(S["sum"].median()) > 1.1 or float(S["sum"].median()) < 0.9:
        b.append(f"Unbalanced design: the shares of a protein sum to {float(S['sum'].median()):.2f} "
                 "(median), correlated terms share part of their variance.")
    # background: the same partition on every quantified protein (the selected ones were
    # picked on between-condition differences, which inflates their factor shares)
    S_ref, _ = variance_partition(run, list(run.abundance.index), F)
    return {"F": F, "source": src, "S": S, "S_ref": S_ref, "terms": terms, "summary": summ,
            "bullets": b}


def design_text(design) -> str:
    if design.ordered:
        return (f"{len(design.groups) or 1} group(s) × ordered {design.level_name} "
                f"({', '.join(design.label(l) for l in design.levels)})")
    return f"{len(design.conditions)} conditions, no order"


def run_reading(stats, stats_path, rule, manifest=None, go_res=None, go_long=None,
                k=None, volcanos=None, max_volcanos=4, design_overrides=None,
                reference_fasta=None, go_q_thr=0.05, jaccard=0.5, groups_from="auto",
                n_boot=50, factors_df=None, factor_regex=None) -> dict | None:
    """stats / rule / manifest / go_res: objects already built by Pathfinder;
    go_long: GOData.long (term intersections). groups_from: 'auto' (ORA groups, else the
    GSEA leading edges when the ORA is empty), 'ora' or 'gsea'. n_boot: replicate
    resamplings for the programme stability (0 = off). factors_df: design factors per
    sample (Pathfinder --design); factor_regex: named groups read in the condition names
    (--factor-regex). Factorial designs (>= 2 factors) get the factor views."""
    run = load_run(stats, stats_path, rule, manifest, reference_fasta=reference_fasta)
    design = detect_design(run.conditions or list(stats.conditions), run.samples,
                           design_overrides)
    kind = design_text(design)
    print(f"[Pathfinder] reading: design {kind} [{design.source}]")
    lookup = id_lookup(run.de)
    groups = None
    if groups_from in ("auto", "ora"):
        groups = _safe("functional groups", build_groups, go_res, go_long, lookup, go_q_thr)
    gsea_note = ""
    if groups is None and groups_from in ("auto", "gsea"):
        groups, gsea_note = _gsea_groups(stats_path, run, manifest, lookup, jaccard,
                                         forced=groups_from == "gsea")
    prog = _safe("programmes", build_programmes, run, design, groups, k)
    if prog is not None and n_boot:
        prog = _safe("programme stability", programme_stability, run, prog, n_boot) or prog
    if prog is None and run.means is not None and len(run.conditions) >= 3:
        print("[Pathfinder] reading: no response programmes (fewer than 20 selected "
              "proteins with abundances).")
    fd = _safe("design factors", _design_factors, run, design, factors_df, factor_regex)

    sections = []

    r = _safe("landscape", fig_landscape, run)
    if r:
        svgs, f = r
        b = [f"{f['n_contrasts']} comparisons; {f['empty']} without any significant protein."]
        b += [f"<b>{html.escape(n)}</b>: {u} up, {d} down." for n, u, d in f["top"]]
        sections.append(dict(id="landscape", title="How much changes", svgs=svgs, bullets=b))

    if fd and fd["blocks"]:
        r = _safe("factor effects", fig_factor_effects, run, fd["F"], fd["cf"], fd["blocks"])
        if r:
            svgs, facts = r
            b = []
            for f in facts:
                lo_, hi_ = f["min"], f["max"]
                name = ("several factors at once" if f["factor"] == "several"
                        else html.escape(str(f["factor"])))
                b.append(f"<b>{name}</b> ({f['n']} comparison(s)): median {f['median']:.0f} "
                         f"differential proteins; most in {html.escape(run.short(hi_[1]))} "
                         f"({hi_[2]} up, {hi_[3]} down)"
                         + (f", fewest in {html.escape(run.short(lo_[1]))} ({lo_[0]})"
                            if f["n"] > 1 else "") + ".")
            n_out = len(run.contrasts) - len(fd["cf"])
            if n_out:
                b.append(f"{n_out} comparison(s) involve a condition outside the factor table "
                         "and are not arranged here.")
            sections.append(dict(id="factor-effects", title="Effects of each design factor",
                                 svgs=svgs, bullets=b))

    r = _safe("condition map", fig_map, run, design,
              fd["F"] if fd else None, fd["roles"] if fd else None)
    if r:
        svgs, f = r
        b = []
        v1, v2 = f.get("var", (None, None))
        if v1 is not None and v2 is not None:
            b.append(f"The first two components carry {v1:.0f} % and {v2:.0f} % of the variance.")
        if f.get("farthest"):
            a, c = f["farthest"]
            b.append(f"Most distant conditions: <b>{html.escape(str(a))}</b> and "
                     f"<b>{html.escape(str(c))}</b>.")
        sp = f.get("factor_spread") or {}
        if sp:
            b.append("Spread of the factor levels on the map (mean distance of their centroids "
                     "to the centre): " + ", ".join(
                         f"{html.escape(k)} {v:.1f}" for k, v in sorted(sp.items(), key=lambda kv: -kv[1]))
                     + ".")
        for g, a, c, _d in f.get("steps", [])[:2]:
            b.append(f"Largest step{'' if g == 'all' else ' in ' + html.escape(str(g))}: "
                     f"{html.escape(str(a))} → {html.escape(str(c))}.")
        sections.append(dict(id="map", title="Where the conditions sit", svgs=svgs, bullets=b))

    r = _safe("programmes", fig_programmes, run, design, prog,
              fd["order"] if fd else None, fd["blockmap"] if fd else None)
    if r:
        svgs, _ = r
        b = []
        for p in range(1, prog.k + 1):
            n = int((prog.labels == p).sum())
            th = prog.themes.get(p, [])
            ths = "; ".join(html.escape(t[0]) for t in th) if th else "no enriched functional group"
            stab = ""
            if prog.stability:
                j = prog.stability[p]
                stab = f" · stability {j:.2f} ({stability_class(j)})"
            b.append(f"<b>P{p}</b> ({n} proteins{stab}): highest in "
                     f"{html.escape(str(prog.high[p]))}, lowest in "
                     f"{html.escape(str(prog.low[p]))} · groups: {ths} · "
                     f"e.g. {_names(prog.names.get(p, []), 4)}.")
        if prog.stability:
            js = list(prog.stability.values())
            n_st = sum(1 for j in js if j >= 0.75)
            b.insert(0, f"Stability under replicate resampling ({prog.n_boot} resamplings, "
                        f"k-means re-run, Jaccard with the original programme): {n_st}/{prog.k} "
                        "stable (≥ 0.75); below 0.6 a programme dissolves and its profile "
                        "should not be read on its own.")
        sections.append(dict(id="programmes", title="Response programmes", svgs=svgs, bullets=b))

    if fd and prog is not None:
        r = _safe("programme profiles by factor", fig_factor_profiles, run, prog, fd["F"],
                  fd["roles"], fd["cf"])
        if r:
            svgs, f = r
            ro = fd["roles"]
            b = [f"x = {html.escape(ro['x'])}, one line per {html.escape(ro['lines'])}"
                 + (f", one panel per {html.escape(ro['facets'][0])}" if ro["facets"] else "")
                 + (f"; {html.escape(f['ref'])} (the reference of the {html.escape(ro['lines'])} "
                    "comparisons) as the neutral line (white at night, black by day)" if f["ref"] else "")
                 + (f"; {', '.join(html.escape(v) for v in f['shared'])} = shared start of every "
                    "panel" if f["shared"] else "") + "."]
            for p, (cmax, vmax_, cmin, vmin_) in f["extremes"].items():
                b.append(f"<b>P{p}</b>: highest in {html.escape(str(cmax))} (z {vmax_:+.1f}), "
                         f"lowest in {html.escape(str(cmin))} (z {vmin_:+.1f}).")
            sections.append(dict(id="factor-profiles", title="Programmes by design factor",
                                 svgs=svgs, bullets=b))

    if fd:
        sl = _safe("slope graph", fig_slope_graph, run, fd["F"], fd["roles"], prog)
        bp = _safe("bump chart", fig_bump_chart, run, prog, fd["F"], fd["roles"])
        if sl or bp:
            svgs = [x[0] for x in (sl, bp) if x]
            b = []
            if sl:
                f = sl[1]
                b.append(f"End to end on {html.escape(f['factor'])}: the {f['n']} proteins "
                         f"moving most from {html.escape(f['from'])} to "
                         f"{html.escape(f['to'])}, {f['up']} up and {f['n'] - f['up']} down. "
                         "Abundances averaged over the other factors, not a test.")
            if bp:
                f = bp[1]
                b.append(f"Rank of the programmes along {html.escape(f['factor'])}"
                         + (f", one panel per level ({', '.join(f['panels'])})"
                            if f["panels"] else "")
                         + ": crossing lines mean the programmes swap order, which an "
                           "amplitude plot hides.")
            sections.append(dict(id="factor-order", title="Order and end-to-end change",
                                 svgs=svgs, bullets=b))

    fact = _safe("variance partition", _factorial, run, prog, fd)
    if fact and fact.get("summary") is not None:
        r = _safe("factor figure", fig_factors, fact["summary"], fact["terms"], fact["F"],
                  fact["source"])
        if r:
            sections.append(dict(id="factors", title="What drives the response",
                                 svgs=r[0], bullets=fact["bullets"]))
        v = _safe("variance violins", fig_variance_violins, fact["S"], fact.get("S_ref"),
                  fact["terms"], fact["source"])
        if v:
            svgs, f = v
            b = ["One pair of half-violins per design term: the distribution over the proteins "
                 "that the median bars above summarise."]
            if f["n_ref"]:
                md = f["medians"]
                gain = sorted(((t, ms - mr) for t, (ms, mr) in md.items()
                               if t not in ("replicates", "higher-order")), key=lambda x: -x[1])
                b.append(f"Background = all {f['n_ref']} quantified proteins without missing "
                         "value. The selected proteins were picked on differences between "
                         "conditions, so a factor matters where its share exceeds the background.")
                up = [f"{html.escape(t)} ({md[t][0] * 100:.0f} % vs {md[t][1] * 100:.0f} %)"
                      for t, g in gain if g >= 0.05][:3]
                b.append(("Above the background (median): " + ", ".join(up) + ".") if up else
                         "No design term exceeds the background by 5 points (median).")
                if "replicates" in md:
                    ms, mr = md["replicates"]
                    b.append(f"Replicate noise: {ms * 100:.0f} % of the variance of a selected "
                             f"protein, {mr * 100:.0f} % of a background protein (median).")
            sections.append(dict(id="variance-violins", title="Variance shares per protein",
                                 svgs=svgs, bullets=b))

    r = _safe("dynamic range", fig_dynamic_range, run, prog)
    if r:
        svg, f = r
        b = [f"{f['n']} proteins over {f['orders']:.1f} orders of magnitude "
             f"({f['range_log2']:.1f} log2).",
             f"{f['frac'] * 100:.0f} % of the proteins are differential somewhere; "
             f"{f['top_decile'] * 100:.0f} % in the most abundant tenth and "
             f"{f['bottom_decile'] * 100:.0f} % in the least abundant tenth."]
        if f["frac"] > 0 and f["top_decile"] > 2.5 * max(f["bottom_decile"], 1e-9):
            b.append("The differential proteins sit mostly at the abundant end: part of the "
                     "signal is detectability, not biology.")
        sections.append(dict(id="dynamic-range", title="Dynamic range", svgs=[svg], bullets=b))

    r = _safe("term treemap", fig_go_treemap, run, groups)
    tree = _safe("term tree", fig_term_tree, run, groups, jaccard)
    if r or tree:
        svgs = [x[0] for x in (r, tree) if x]
        b = []
        if r:
            f = r[1]
            b.append(f"{f['n_terms']} enriched terms, each tile sized by its −log10(q) and "
                     f"coloured by its Meta-GO group ({f['n_groups']} groups shown).")
            th, nm, share = f["largest"]
            b.append(f"Largest block: G{th + 1} · {html.escape(nm)}, {share * 100:.0f} % of the "
                     "total area — its terms share most of their proteins, so that area is one "
                     "piece of evidence, not many.")
        if tree:
            b.append(f"The tree below is the one the groups come from: terms joined by "
                     f"1 − Jaccard on their proteins, average linkage, cut at {jaccard:g}.")
        sections.append(dict(id="go-terms", title="Enriched terms and their redundancy",
                             svgs=svgs, bullets=b))

    r = _safe("functional groups", fig_groups, run, groups)
    if r:
        svgs, f = r
        if groups.origin == "GSEA":
            lead_in = ("Groups requested from the GSEA (--groups-from gsea)"
                       if gsea_note.startswith("requested") else
                       "The ORA gave no usable term: the groups come from the <b>GSEA</b>")
            b = [f"{lead_in} of Deimos. {f['n_terms']} pathways (FDR as in the GSEA_* sheets) in "
                 f"{f['n_groups']} groups of pathways sharing their leading-edge proteins "
                 f"(Jaccard, cut {jaccard:g}); the {f['shown']} most significant are drawn. "
                 "Direction = NES; members = leading-edge proteins."]
        else:
            b = [f"{f['n_terms']} enriched terms in {f['n_groups']} functional groups (the "
                 f"Meta-GO redundancy groups: terms sharing proteins); the {f['shown']} most "
                 "significant are drawn."]
        b += [f"<b>{html.escape(n)}</b>: strongest in {html.escape(', '.join(cs))}."
              for n, cs in f["top"]]
        if groups is not None and groups.n_unmapped:
            b.append(f"{groups.n_unmapped} accession(s) of the GO intersections match no protein "
                     "of the workbook and are left out of the groups.")
        sections.append(dict(id="groups", title="Functional groups", svgs=svgs, bullets=b))

    if volcanos:
        vc = [c for c in volcanos if c in run.contrasts]
        miss = [c for c in volcanos if c not in run.contrasts]
        if miss:
            print(f"[Pathfinder] ⚠ reading: --volcano not found: {', '.join(miss)}")
    else:
        counts = {c: int(run.sig(c).fillna(False).sum()) for c in run.contrasts}
        vc = [c for c, n in sorted(counts.items(), key=lambda x: -x[1]) if n > 0][:max_volcanos]
    for i, c in enumerate(vc, 1):
        r = _safe(f"volcano {c}", fig_volcano, run, c, groups, prog, i)
        if r:
            svgs, f = r
            b = [f"{f['up']} up, {f['down']} down."]
            if f["up_names"]:
                b.append(f"Most marked up: {_names(f['up_names'])}.")
            if f["down_names"]:
                b.append(f"Most marked down: {_names(f['down_names'])}.")
            if f["groups"]:
                b.append("Coloured functional groups: "
                         + "; ".join(html.escape(t) for t in f["groups"]) + ".")
            sections.append(dict(id=f"volcano{i}", title=f"Volcano · {f['contrast']}", svgs=svgs,
                                 bullets=b))

    P = _safe("contrast pairs", fc_pairs, run) if len(run.contrasts) >= 2 else None
    if P is not None and not P.empty:
        r = _safe("FC vs FC", fig_fc_vs_fc, run, P, prog)
        if r:
            svgs, f = r
            b = ["log2FC against log2FC: a protein away from the diagonal responds to one "
                 "contrast more than to the other."]
            best = f["best"]
            if best is not None and best.shared >= 5 and best.slope == best.slope:
                sl = float(best.slope)
                b.append(f"{run.short(best.contrast_y)} against {run.short(best.contrast_x)}: "
                         f"{best.shared} shared DEP, major-axis slope {sl:+.2f}"
                         + (" — same amplitude." if 0.8 <= abs(sl) <= 1.25 else
                            " — the y contrast moves them less." if abs(sl) < 0.8 else
                            " — the y contrast moves them further."))
                b.append(f"{best.up_up + best.down_down} of the shared DEP move the same way, "
                         f"{best.up_down + best.down_up} the opposite way.")
            elif best is not None:
                b.append(f"The widest pair shares {best.shared} DEP, too few for a slope: "
                         "the contrasts hit different proteins.")
            b.append("The slope is fitted on the shared DEP only. On the union, proteins "
                     "significant in one contrast alone form a cross along the axes and pull "
                     "the slope towards that axis.")
            if f["n_pairs"] > len(f["picked"]):
                b.append(f"{f['n_pairs']} pairs in all (sheet FC_vs_FC); the matrix shows every "
                         f"pair, the scatters the {len(f['picked'])} sharing the most DEP.")
            sections.append(dict(id="fc-vs-fc", title="Shared and divergent responses",
                                 svgs=svgs, bullets=b))

    r = _safe("radar", fig_radar, run, prog, design, fd["F"] if fd else None)
    if r:
        svg, f = r
        sections.append(dict(
            id="radar", title="Shape of each programme",
            svgs=[svg],
            bullets=[f"The {f['k']} programmes drawn over the {f['n_cond']} conditions on one "
                     "common z scale: the shape says where each one is high and where it is "
                     "flat.",
                     "Only drawn for a design with at most eight conditions; beyond that the "
                     "programme profiles above read better."]))

    pg = _safe("programmes and groups", fig_programme_groups, run, prog, groups)
    if pg:
        svg, f = pg
        p_, g_, k_, q_ = f["best"]
        b = [f"{f['n_links']} programme–group pairs sharing at least two proteins, "
             f"{f['n_sig']} more than expected by chance.",
             f"Strongest link: P{p_} with {html.escape(g_)} ({k_} proteins, q {q_:.1g})."]
        multi = (f["table"][f["table"].q < 0.05].groupby("group").programme.nunique())
        spread = [g for g, n in multi.items() if n > 1]
        if spread:
            b.append("Tied to more than one programme: "
                     + "; ".join(html.escape(str(g)) for g in spread[:3])
                     + " — the clustering does not split them.")
        sections.append(dict(id="programme-groups", title="Programmes and functional groups",
                             svgs=[svg], bullets=b))
        prog_groups = f["table"]
    else:
        prog_groups = None

    r = _safe("spotlight", fig_spotlight, run, design, prog)
    if r:
        svgs, f = r
        sections.append(dict(id="spotlight", title="Spotlight proteins", svgs=svgs,
                             bullets=[f"Shown: {_names(f['proteins'], 8)}."]))

    n_sig = int(sum(run.sig(c).fillna(False) for c in run.contrasts).gt(0).sum()) \
        if run.contrasts else 0
    kpis = [("proteins quantified", len(run.de)), ("comparisons", len(run.contrasts)),
            ("differential (any comparison)", n_sig)]
    if prog is not None:
        kpis.append(("response programmes", prog.k))
    if groups is not None:
        kpis.append(("functional groups", len(groups.names)))
    meta = {
        "reading_design": f"{kind} [{design.source}]",
        "reading_abundance": run.abundance_kind or "n/a",
        "reading_protein_names": run.label_source,
        "reading_programmes": (f"k = {prog.k} ({'set by --k' if k else 'silhouette, 3-8'}); "
                               f"k-means on z-scored condition means of the {len(prog.z)} "
                               "selected proteins; functional groups enriched by hypergeometric "
                               f"test, BH q < {THEME_Q_THR:g}") if prog is not None else "none",
        "reading_programme_stability": (
            f"{prog.n_boot} replicate resamplings (bootstrap within condition), k-means re-run, "
            "Hungarian matching, mean Jaccard per programme (>= 0.75 stable, < 0.6 dissolved)"
            if prog is not None and prog.stability else "not computed"),
        "reading_groups": ("none (no GO, no GSEA)" if groups is None else
                           ("GSEA leading edges (" + gsea_note + "), pathways grouped by "
                            f"Jaccard on their leading-edge proteins (cut {jaccard:g})"
                            if groups.origin == "GSEA" else
                            "Meta-GO redundancy groups (Jaccard, --jaccard), named by their "
                            "representative term")),
    }
    if fact and fact.get("summary") is not None:
        ro = fd["roles"]
        meta["reading_factor_roles"] = (f"x = {ro['x']}, lines = {ro['lines']}, panels = "
                                        f"{', '.join(ro['facets']) or '-'}")
        meta["reading_factors"] = (f"{', '.join(fact['F'].columns)} (from {fact['source']}); "
                                   "type-II sums of squares of a main-effect + two-factor "
                                   "interaction model on per-sample log2 abundances; "
                                   "'higher-order' = what condition means add; 'replicates' = "
                                   "within-condition residual")
    tables = _safe("tables", reading_tables, run, design, prog, groups, sections, kpis,
                   fact, P, prog_groups) or {}
    return {"sections": sections, "kpis": kpis, "meta": meta, "tables": tables,
            "rule_text": run.rule_text(),
            "prog": prog, "groups": groups, "design": design, "design_text": kind, "run": run}
