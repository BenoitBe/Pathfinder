#!/usr/bin/env python3
"""
Pathfinder - exploration & synthesis layer on top of Deimos/Phobos outputs.

Reads the upstream Excel report (no statistics are recomputed) and produces:
  - Meta-GO      : term x contrast matrix, cross-directionality, redundancy grouping
  - Exploration  : PCA, PC <-> factor association (multi-factor), correlation, UMAP
  - Publications : STRING PMID enrichment per contrast (online, opt-in)
  - Summary      : deterministic interpretive summary

Outputs: pathfinder_report.html + pathfinder_synthesis.xlsx + pathfinder_summary.md

Usage:
  python pathfinder.py [--stats FILE] [--go FILE] [--design FILE] [--outdir DIR]
                       [--pthr P] [--jaccard D] [--no-umap] [--diagnose]
                       [--org hsapiens | --no-online] [--refresh]
"""
from __future__ import annotations
import argparse
import os
import sys
import time

import pandas as pd

import io_layer
import axis_meta_go
import axis_explore
import axis_proteins
import axis_wgcna
import axis_publications
import summary as summary_mod
import species as species_mod
import build_report_html as build_html
import build_report_xlsx as build_xlsx

__version__ = "1.1.1"


def main(argv=None):
    ap = argparse.ArgumentParser(prog="pathfinder")
    ap.add_argument("--version", action="version", version=f"Pathfinder {__version__}")
    ap.add_argument("--stats", default=None,
                    help="statistics workbook (auto-detected by content if omitted)")
    ap.add_argument("--go", default=None,
                    help="GO enrichment workbook (optional; may be the same file). "
                         "'none' disables the Meta-GO axis and any GO auto-detection.")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--pthr", type=float, default=0.05,
                    help="adjusted-p threshold: protein signature (p.adj, no fold-change "
                         "criterion) and GO matrix mask")
    ap.add_argument("--jaccard", type=float, default=0.5, help="cut distance for GO redundancy grouping")
    ap.add_argument("--no-umap", action="store_true")
    ap.add_argument("--design", default=None,
                    help="ExperimentalDesign (CSV/xlsx) mapping sample -> condition (+ extra factors)")
    ap.add_argument("--org", default=None,
                    help="species for the Publications axis: gProfiler organism (hsapiens...), "
                         "'Homo sapiens', 'human', or an NCBI taxid. Skips the interactive prompt.")
    ap.add_argument("--no-online", action="store_true",
                    help="fully offline: no species prompt, no STRING query")
    ap.add_argument("--refresh", action="store_true", help="bypass the cache and re-issue online calls")
    ap.add_argument("--write-design-template", action="store_true",
                    help="write a pre-filled ExperimentalDesign_template.csv, then exit")
    ap.add_argument("--diagnose", action="store_true",
                    help="report detected sheets/columns and group sizes, then exit")
    args = ap.parse_args(argv)

    # --go none: no Meta-GO, and above all no hunting for a GO workbook around
    # (it could belong to another project).
    go_disabled = str(args.go).strip().lower() == "none" if args.go else False
    if go_disabled:
        args.go = None

    os.makedirs(args.outdir, exist_ok=True)

    # CONTENT-based auto-detection (independent of filenames)
    import glob as _glob

    def _list_xlsx():
        seen, out = set(), []
        for d in [args.outdir, ".", ".."]:
            for p in sorted(_glob.glob(os.path.join(d, "*.xlsx"))):
                rp = os.path.realpath(p)
                if rp not in seen:
                    seen.add(rp); out.append(p)
        for p in sorted(_glob.glob(os.path.join(".", "**", "*.xlsx"), recursive=True)):
            rp = os.path.realpath(p)
            if rp not in seen:
                seen.add(rp); out.append(p)
        # exclude our own outputs
        return [p for p in out if not os.path.basename(p).lower().startswith("pathfinder_")]

    def _peek(path):
        try:
            import openpyxl
            wb = openpyxl.load_workbook(path, read_only=True)
        except Exception:
            return False, False
        stats = go = False
        for ws in wb.worksheets:
            hdr = next(ws.iter_rows(min_row=1, max_row=1, values_only=True), ())
            hdr = [str(h) for h in hdr if h is not None]
            if any(h.endswith("_p.val") for h in hdr):
                stats = True
            if str(ws.title).lower().startswith("go") and "term_name" in hdr:
                go = True
        wb.close()
        return stats, go

    candidates = None
    if not args.stats or (not args.go and not go_disabled):
        candidates = [(p, _peek(p)) for p in _list_xlsx()]

    if not args.stats:
        for p, (is_stats, _) in candidates:
            if is_stats:
                args.stats = p
                print(f"[Pathfinder] statistics auto-detected (by content): {p}")
                break
        if not args.stats:
            ap.error("no workbook with '<C>_p.val' columns found - specify --stats <file.xlsx>")
    if go_disabled:
        print("[Pathfinder] GO disabled (--go none): Meta-GO axis skipped.")
    elif not args.go:
        # the GO sheets of the analysed workbook always win over any other file around
        if args.stats and _peek(args.stats)[1]:
            args.go = args.stats
            print("[Pathfinder] GO found inside the statistics workbook.")
        else:
            for p, (_, is_go) in (candidates or []):
                if is_go:
                    args.go = p
                    print(f"[Pathfinder] GO auto-detected (by content): {p}")
                    break

    provenance = {"tool": f"Pathfinder {__version__}",
                  "date": time.strftime("%Y-%m-%d %H:%M:%S")}

    # ---- Experimental design (sample -> condition) ----
    design_map, design_src, design_frame = None, None, None
    if args.design and os.path.exists(args.design):
        design_map, design_src = io_layer.load_design(args.design), args.design
        design_frame = io_layer.load_design_frame(args.design)
    else:
        dcsv = None
        for d in [args.outdir, ".", ".."]:
            for pat in ("ExperimentalDesign*.csv", "*experimentaldesign*.csv",
                        "*design*.csv", "*design*.tsv", "*metadata*.csv"):
                hits = sorted(_glob.glob(os.path.join(d, pat)))
                if hits:
                    dcsv = hits[0]; break
            if dcsv:
                break
        if dcsv:
            design_map, design_src = io_layer.load_design(dcsv), dcsv
            design_frame = io_layer.load_design_frame(dcsv)
        else:                                   # design inside a sheet of the workbook?
            try:
                for nm, df in pd.read_excel(args.stats, sheet_name=None).items():
                    cols = [str(c).lower() for c in df.columns]
                    by_name = any(k in str(nm).lower() for k in ("design", "metadata", "experiment", "samples"))
                    # by content: a 'condition/group...' column + few rows (~ n samples)
                    by_content = (any(any(k in c for k in ("condition", "group", "groupe", "class", "treatment"))
                                      for c in cols) and 1 <= len(df) <= 500 and df.shape[1] <= 12)
                    if by_name or by_content:
                        design_map, design_src = io_layer.load_design(df), f"{os.path.basename(args.stats)}#{nm}"
                        design_frame = io_layer.load_design_frame(df)
                        break
            except Exception:
                pass
    if design_src:
        print(f"[Pathfinder] design detected: {design_src}")

    # ---- Loading ----
    stats = io_layer.load_stats(args.stats, design_map=design_map)
    if not stats.conditions_resolved:
        print("[Pathfinder] ⚠ conditions NOT resolved: the sample columns "
              f"({', '.join(map(str, stats.sample_int_cols[:4]))}...) do not encode the condition "
              "and no design was found. Provide --design ExperimentalDesign.csv "
              "(sample + condition columns) to colour the PCA and test PC <-> condition.")

    if args.write_design_template:
        tpl = os.path.join(args.outdir, "ExperimentalDesign_template.csv")
        pd.DataFrame({"sample": stats.sample_int_cols,
                      "condition": ["TO_FILL"] * len(stats.sample_int_cols)}).to_csv(tpl, index=False)
        print(f"[Pathfinder] template written: {tpl}\n"
              f"  -> fill in the 'condition' column ({len(stats.sample_int_cols)} rows), then run:\n"
              f"     python pathfinder.py --design {os.path.basename(tpl)}")
        return None

    if args.diagnose:
        from collections import Counter
        print("\n=== Pathfinder DIAGNOSTIC ===")
        print(f"Statistics file: {args.stats}\n")
        # structure of the WHOLE workbook
        allsheets = pd.read_excel(args.stats, sheet_name=None)
        print("Workbook sheets:")
        for name, df in allsheets.items():
            num = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
            has_pval = any(str(c).endswith("_p.val") for c in df.columns)
            has_term = "term_name" in [str(c) for c in df.columns]
            # numeric columns outside known stats patterns (sample candidates)
            statish = [c for c in df.columns if any(str(c).endswith(s) for s in ("_diff","_p.val","_p.adj"))
                       or str(c).startswith(("Robustness_Score_","Pi_Score_"))]
            sample_like = [c for c in num if c not in statish and str(c) not in io_layer._META]
            flags = []
            if has_pval: flags.append("STATS(_p.val)")
            if has_term: flags.append("GO(term_name)")
            print(f"  · {name:<28} lignes={len(df):>5}  cols={len(df.columns):>4}  "
                  f"numeric={len(num):>4}  sample_candidates={len(sample_like):>4}  {' '.join(flags)}")
            if sample_like:
                ex = ', '.join(str(c) for c in sample_like[:8])
                print(f"      e.g. sample columns: {ex}{' ...' if len(sample_like)>8 else ''}")
        print()
        # what the loader actually kept
        sizes = Counter(stats.cond_map[c] for c in stats.sample_int_cols)
        print(f"[Loader] comparison sheet selected: {len(stats.comp)} rows")
        print(f"[Loader] PCA MATRIX: sheet '{stats.pca_source}' -> {len(stats.sample_int_cols)} "
              f"sample columns, {len(stats.pca_df)} proteins"
              + (" [Z-SCORED]" if stats.pca_is_zscored else " [raw intensities]"))
        if getattr(stats, "candidates", None):
            print("         candidate matrices considered:")
            for c in sorted(stats.candidates, key=lambda c: -c["n_rows"]):
                tags = []
                if c["zscored"]: tags.append("z-scored")
                if not c["log_scale"]: tags.append("not log-scaled")
                if c["match"]: tags.append(f"design match {c['match']:.0%}")
                mark = " <== used" if c["sheet"] == stats.pca_source else ""
                print(f"           {c['sheet']:<26} {c['n_cols']:>3} cols x {c['n_rows']:>5} rows"
                      f"  {', '.join(tags) or 'raw log2'}{mark}")
        if stats.sample_int_cols:
            print("         e.g.: " + ", ".join(str(c) for c in stats.sample_int_cols[:10]))
        print(f"[Loader] Z-score sheet: {'present' if stats.zscore is not None else 'absent'}")
        print(f"[Loader] contrasts ({len(stats.contrasts)}): {', '.join(stats.contrasts[:6])}"
              f"{' ...' if len(stats.contrasts)>6 else ''}")
        print(f"[Loader] PCA intensity columns kept: {len(stats.sample_int_cols)}")
        print(f"[Loader] conditions ({len(sizes)}) group sizes: "
              + ", ".join(f"{c}={n}" for c, n in sizes.items()))
        from collections import Counter
        sizes = Counter(stats.cond_map[c] for c in stats.sample_int_cols)
        print(f"[Loader] design: {design_src or 'not provided'} -> conditions "
              f"{'RESOLVED' if stats.conditions_resolved else 'NOT resolved'}")
        print(f"[Loader] conditions ({len(sizes)}) group sizes: "
              + ", ".join(f"{c}={n}" for c, n in list(sizes.items())[:12]))
        # GO sheet -> contrast mapping audit
        if args.go and os.path.exists(args.go):
            try:
                g = io_layer.load_go(args.go, stats.contrasts, conditions=stats.conditions, verbose=False)
                print(f"\n[GO mapping] {len(g.mapping)} GO sheet(s) -> contrasts")
                for sh, c in g.mapping.items():
                    flag = "OK " if c else "!! "
                    print(f"    {flag}{sh:<28} -> {c or 'UNMATCHED'}   [{g.mapping_how[sh]}]")
                if g.unmatched:
                    # explain WHY, instead of leaving the user guessing
                    aliases = io_layer._condition_aliases(stats.conditions)
                    print("\n    Why unmatched? Pathfinder splits a sheet name into two known")
                    print("    conditions. Known spellings (aliases -> condition):")
                    for a, c in sorted(aliases.items()):
                        print(f"      {a:<24} -> {c}")
                    print("    Contrasts to match against:")
                    for c in stats.contrasts[:4]:
                        print(f"      {c}")
                    print("    If a sheet name uses tokens absent from that list, the design "
                          "conditions and the stats columns disagree: fix the design, or "
                          "rename the sheets.")
                if g.duplicates:
                    for c, shs in g.duplicates.items():
                        print(f"    !! COLLISION on '{c}': {', '.join(shs)}")
                if g.missing:
                    print(f"    note: contrasts without a GO sheet: {', '.join(map(str, g.missing))}")
                print(f"    summary: {len(g.mapping) - len(g.unmatched)}/{len(g.mapping)} mapped, "
                      f"{len(g.duplicates)} collision(s), {len(g.missing)} contrast(s) without GO")
            except Exception as e:
                print(f"\n[GO mapping] not available ({type(e).__name__}: {e})")
        print("=============================\n")
        return None

    provenance["stats_file"] = os.path.basename(args.stats)
    # ---- Upstream manifest (Deimos contract): version check + what it declares ----
    manifest = io_layer.read_manifest(args.stats)
    if manifest:
        sch, run = manifest.get("schema", {}), manifest.get("run", {})
        provenance["upstream_manifest"] = (f"{sch.get('name', '?')} {sch.get('version', '?')}"
                                           + (f" · commit {run['deimos_commit']}"
                                              if run.get("deimos_commit") else ""))
        warn = io_layer.manifest_version_warning(manifest)
        if warn:
            print(f"[Pathfinder] ⚠ {warn}")
            provenance["upstream_manifest_warning"] = warn
        declared = manifest.get("matrix", {}).get("sheet")
        if declared and declared != stats.pca_source:
            print(f"[Pathfinder] ⚠ the upstream manifest names '{declared}' as the sample "
                  f"matrix, but '{stats.pca_source}' was selected - check --diagnose.")
        if manifest.get("stats", {}).get("p_adjust"):
            provenance["p_adjust_upstream"] = manifest["stats"]["p_adjust"]
        if manifest.get("stats", {}).get("volcano_rule"):
            provenance["volcano_rule_upstream"] = manifest["stats"]["volcano_rule"]
    provenance["contrasts"] = ", ".join(stats.contrasts)
    provenance["conditions"] = ", ".join(stats.conditions)
    provenance["design"] = design_src or "not provided"
    provenance["conditions_resolved"] = "yes" if stats.conditions_resolved else "NO - provide --design"
    provenance["n_samples"] = len(stats.sample_int_cols)
    provenance["n_proteins"] = len(stats.comp)
    provenance["pca_matrix"] = (f"{stats.pca_source} ({len(stats.pca_df)} proteins x "
                                f"{len(stats.sample_int_cols)} samples"
                                + (", z-scored subset" if stats.pca_is_zscored else ", raw log2") + ")")

    # Contrasts whose directions can be compared (same comparison, different strata):
    # the only pairs where "up in one, down in the other" is a genuine reversal.
    pairs = axis_proteins.comparable_pairs(stats)
    if pairs is not None:
        provenance["opposite_rule"] = (
            f"reversal only between contrasts sharing no condition with |r| >= "
            f"{axis_proteins.PARALLEL_R_MIN} between their log2FC profiles "
            f"({len(pairs)} comparable pair(s)); up in one contrast and down in another "
            "that compares different conditions is reported as 'mixed'")

    go_res = None
    if args.go and os.path.exists(args.go):
        try:
            go = io_layer.load_go(args.go, stats.contrasts, conditions=stats.conditions)
        except Exception as e:
            print(f"[Pathfinder] ⚠ GO file unreadable ({type(e).__name__}: {e}) — "
                  "Meta-GO axis skipped.")
            provenance["go_file"] = f"{os.path.basename(args.go)} - unreadable, axis skipped"
            go = None
        if go is not None:
            provenance["go_sheets_mapped"] = (f"{len(go.mapping) - len(go.unmatched)}/{len(go.mapping)}"
                                              + (f" ({go.from_manifest} via the upstream manifest)"
                                                 if getattr(go, "from_manifest", 0) else ""))
            if go.unmatched:
                provenance["go_unmatched"] = ", ".join(go.unmatched)
            if go.duplicates:
                provenance["go_collisions"] = "; ".join(f"{c}: {', '.join(sh)}"
                                                        for c, sh in go.duplicates.items())
            if go.missing:
                provenance["go_contrasts_without_sheet"] = ", ".join(map(str, go.missing))
            if len(go.long) == 0:
                print("[Pathfinder] \u26a0 no GO sheet matches these contrasts - Meta-GO axis "
                      "skipped. Check that the GO file belongs to this analysis (--go).")
                provenance["go_file"] = f"{os.path.basename(args.go)} - no matching sheet, axis skipped"
            else:
                provenance["go_file"] = os.path.basename(args.go)
                provenance["go_rows"] = len(go.long)
                go_res = axis_meta_go.run(go.long, p_thr=args.pthr, jaccard_cut=args.jaccard,
                                          pairs=pairs)
    if go_res is None:
        # empty matrix -> informative Meta-GO tab
        empty = pd.DataFrame(columns=["term_name", "term_id", "source", "n_sig",
                                      "redundancy_cluster", "is_representative"])
        go_res = {"matrix": pd.DataFrame(), "meta": empty,
                  "directionality": pd.DataFrame(), "term_genes": {},
                  "contrast_corr": pd.DataFrame(), "contrast_order": [],
                  "contrast_shared": pd.DataFrame(), "contrast_relation": None}
        provenance.setdefault("go_file", "disabled (--go none) - Meta-GO axis skipped"
                              if go_disabled else "absent - Meta-GO axis skipped")

    # ---- Exploration axis ----
    # design factors (multi-factor: condition, batch, injection order...)
    factors_df = None
    if design_frame is not None:
        try:
            factors_df = io_layer.map_design_factors(stats.sample_int_cols, design_frame[0], design_frame[1])
            if factors_df is not None:
                print(f"[Pathfinder] design factors tested: {', '.join(map(str, factors_df.columns))}")
            dropped = [str(c) for c in design_frame[0].columns if io_layer.is_replicate_index(c)]
            if dropped:
                provenance["design_not_tested"] = (f"{', '.join(dropped)} (replicate index "
                                                   "within a condition, not a factor)")
        except Exception as e:
            print(f"[Pathfinder] ⚠ design factors unusable ({type(e).__name__})")

    explore = axis_explore.run(stats, do_umap=not args.no_umap, factors_df=factors_df)

    # ---- Protein axis (which proteins carry the signature) ----
    prots = axis_proteins.run(stats, p_thr=args.pthr, pairs=pairs)
    # ---- WGCNA axis (co-expression structure computed upstream) ----
    wgcna = None
    try:
        _xls = pd.read_excel(args.stats, sheet_name=None)
        core_ids = list(prots["table"].head(200).index) if not prots["matrix"].empty else []
        wgcna = axis_wgcna.run(_xls, stats, core_ids=core_ids)
        if wgcna:
            # the set actually tested for module enrichment (reported as such)
            wgcna["core_set_n"] = len(core_ids)
            wgcna["core_set_rule"] = (f"top {len(core_ids)} proteins ranked by recurrence "
                                      "then Pi-score")
    except Exception as e:
        print(f"[Pathfinder] \u26a0 WGCNA sheets unusable ({type(e).__name__}) - axis skipped.")
    if wgcna:
        print(f"[Pathfinder] WGCNA: {wgcna['n_modules']} module(s) read from the workbook")
        provenance["wgcna_modules"] = f"{wgcna['n_modules']} modules (read, not recomputed)"

    provenance["proteins_significant"] = (
        f"{prots['n_significant']} proteins significant in >=1 contrast (p.adj < {args.pthr})"
        + (" - ranked by recurrence then Pi-score" if prots.get("has_pi") else ""))
    provenance["significance_rule"] = (
        f"p.adj < {args.pthr}, no fold-change criterion"
        + (f" ({provenance['p_adjust_upstream']}, as computed upstream)"
           if provenance.get("p_adjust_upstream") else ""))

    # ---- Species (taxid) + Publications axis (online, opt-in) ----
    org_arg = args.org
    excel_org = None
    if not org_arg and not args.no_online:
        excel_org = species_mod.detect_org_in_workbook(args.stats)
        if excel_org:
            print(f"[Pathfinder] species detected in workbook: {excel_org}")

    taxid, org, src = species_mod.resolve_species(org_arg=org_arg)

    def _via_string(query):
        """Resolve a free-text organism against STRING's own species list (~12k),
        which also proves the Publications axis can run for it."""
        try:
            cands = species_mod.lookup_species(query)
        except Exception as e:
            print(f"[Pathfinder] \u26a0 could not reach the STRING species list "
                  f"({type(e).__name__}) - pass an NCBI taxid instead.")
            return None
        return species_mod.choose_candidate(cands)

    if not taxid and org_arg:
        # explicit --org: never prompt for the species, just try to resolve it
        taxid = _via_string(org_arg)
        src = "STRING species list" if taxid else f"unresolved '{org_arg}'"
        if not taxid:
            print(f"[Pathfinder] \u26a0 organism '{org_arg}' not resolved - Publications "
                  "axis skipped. Give the NCBI taxid (e.g. --org 90675) to be sure.")

    elif not taxid and not args.no_online:
        answer = species_mod.prompt_species(default=excel_org)
        if answer is None:                       # not a TTY -> use detection silently
            taxid, org, src = species_mod.resolve_species(excel_org=excel_org)
            answer = excel_org or ""
        elif answer:
            taxid, org, src = species_mod.resolve_species(org_arg=answer)
        else:
            src = "skipped by user"
        if not taxid and answer:
            taxid = _via_string(answer)
            src = "STRING species list" if taxid else f"unresolved '{answer}'"
            if not taxid:
                print(f"[Pathfinder] \u26a0 organism '{answer}' not resolved - Publications "
                      "axis skipped. Give the NCBI taxid (e.g. --org 90675) to be sure.")

    provenance["species_taxid"] = f"{taxid} ({org or 'taxid'} · {src})" if taxid else src

    pubs = None
    run_online = bool(taxid) and not args.no_online
    if run_online:
        print(f"[Pathfinder] Publications axis: querying STRING PMID (taxid {taxid})...")
        pubs = axis_publications.run(stats, taxid, p_thr=args.pthr, refresh=args.refresh)
        n = sum(len(v) for v in pubs["results"].values())
        if n == 0 and not pubs["results"]:
            print("[Pathfinder] ⚠ no publication returned (network issue or no significant proteins).")
        provenance["publications"] = f"STRING PMID · {len(pubs['results'])} contrasts · {n} PMIDs"
        for c, info in pubs["provenance"].items():
            provenance["publications_criterion" if c == "_criterion" else
                       "publications_skipped" if c == "_skipped" else f"pub[{c}]"] = info

    # ---- KPIs ----
    n_go_sig = int(go_res["matrix"].notna().sum().sum()) if not go_res["matrix"].empty else 0
    kpis = {
        "Contrasts": len(stats.contrasts),
        "Conditions": len(stats.conditions),
        "Samples": len(stats.sample_int_cols),
        "Proteins": len(stats.comp),
        "Significant GO cells": n_go_sig,
        "Significant proteins": prots["n_significant"],
        "PC1 variance": f"{explore['pca']['var_exp'][0]:.0f}%",
    }

    # ---- Interpretive summary ----
    resume = summary_mod.generate(stats, go_res, explore, provenance, prots, wgcna)
    md_path = os.path.join(args.outdir, "pathfinder_summary.md")
    try:
        with open(md_path, "w", encoding="utf-8") as f:
            f.write(resume["markdown"])
    except PermissionError:
        import time as _t
        md_path = md_path.replace(".md", f"_{_t.strftime('%H%M%S')}.md")
        with open(md_path, "w", encoding="utf-8") as f:
            f.write(resume["markdown"])

    # ---- Outputs (robust to Windows file locks) ----
    def _safe_write(writer_fn, path, label):
        try:
            writer_fn(path)
            return path
        except PermissionError:
            import time as _t
            alt = path.replace(".", f"_{_t.strftime('%H%M%S')}.", 1)
            print(f"[Pathfinder] ⚠ {os.path.basename(path)} is locked (open in another program?) "
                  f"-> writing to {os.path.basename(alt)}")
            writer_fn(alt)
            return alt

    xlsx_path = _safe_write(lambda p: build_xlsx.build(p, go_res, explore, provenance, resume, pubs, prots, wgcna),
                            os.path.join(args.outdir, "pathfinder_synthesis.xlsx"), "xlsx")
    # branding: inline SVG mark (PNG base64 fallback)
    import base64
    _here = os.path.dirname(os.path.abspath(__file__))
    icon_svg, icon_b64 = None, None
    _svg = os.path.join(_here, "pathfinder_mark.svg")
    if os.path.exists(_svg):
        icon_svg = open(_svg, encoding="utf-8").read()
        icon_svg = icon_svg[icon_svg.index("<svg"):]
    else:
        for _nm in ("pathfinder_icon_preview.png", "pathfinder_icon_128.png"):
            _p = os.path.join(_here, _nm)
            if os.path.exists(_p):
                icon_b64 = base64.b64encode(open(_p, "rb").read()).decode()
                break
    branding = {"icon_b64": icon_b64, "icon_svg": icon_svg}

    html_path = _safe_write(lambda p: build_html.build(p, go_res, explore, provenance, kpis, resume, branding, pubs, prots, wgcna),
                            os.path.join(args.outdir, "pathfinder_report.html"), "html")

    print("\n" + resume["markdown"] + "\n")
    print(f"[Pathfinder] OK\n  - {xlsx_path}\n  - {html_path}\n  - {md_path}")
    return xlsx_path, html_path


if __name__ == "__main__":
    sys.exit(0 if main() else 0)
