#!/usr/bin/env python3
"""
Pathfinder — couche d'exploration & de synthèse au-dessus des sorties Deimos.

Étape 1 (offline, déterministe) :
  - Axe B : méta-comparaison GO (matrice terme × contraste, directionnalité, redondance)
  - Axe C : exploration (PCA enrichie, association PC↔facteur, corrélation, UMAP gardé)
  - Sorties : pathfinder_synthese.xlsx + pathfinder_report.html

Axe A (publications, online + cache) : branché à l'étape 2 via species.py / cache.py.

Usage :
  python pathfinder.py --stats rapportstatistique.xlsx [--go rapportenrichissementGO.xlsx]
                       [--outdir .] [--pthr 0.05] [--jaccard 0.5] [--no-umap]
                       [--org rnorvegicus]   # transmis par Deimos, réservé axe A
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
import axis_publications
import summary as summary_mod
import species as species_mod
from report import build_xlsx, build_html


def main(argv=None):
    ap = argparse.ArgumentParser(prog="pathfinder")
    ap.add_argument("--stats", default=None,
                    help="rapportstatistique.xlsx (auto-détecté dans --outdir si absent)")
    ap.add_argument("--go", default=None, help="rapportenrichissementGO.xlsx (optionnel)")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--pthr", type=float, default=0.05, help="seuil p pour masquer la matrice GO")
    ap.add_argument("--jaccard", type=float, default=0.5, help="distance de coupe redondance GO")
    ap.add_argument("--no-umap", action="store_true")
    ap.add_argument("--design", default=None,
                    help="ExperimentalDesign (CSV/xlsx) mappant échantillon→condition")
    ap.add_argument("--org", default=None, help="org gprofiler (hsapiens…) ou taxid — requis pour --online")
    ap.add_argument("--online", action="store_true",
                    help="active l'axe A (publications STRING PMID) — nécessite un réseau et --org")
    ap.add_argument("--refresh", action="store_true", help="ignore le cache et refait les appels online")
    ap.add_argument("--write-design-template", action="store_true",
                    help="écrit un ExperimentalDesign_template.csv pré-rempli avec les colonnes détectées, puis quitte")
    ap.add_argument("--diagnose", action="store_true",
                    help="affiche les colonnes/conditions détectées puis quitte")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)

    # auto-détection par CONTENU (indépendante du nom de fichier)
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
        # exclure nos propres sorties
        return [p for p in out if "pathfinder_synthese" not in os.path.basename(p).lower()]

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
    if not args.stats or not args.go:
        candidates = [(p, _peek(p)) for p in _list_xlsx()]

    if not args.stats:
        for p, (is_stats, _) in candidates:
            if is_stats:
                args.stats = p
                print(f"[Pathfinder] stats auto-détecté (par contenu) : {p}")
                break
        if not args.stats:
            ap.error("aucun classeur avec colonnes '<C>_p.val' trouvé — précise --stats <fichier.xlsx>")
    if not args.go:
        for p, (_, is_go) in (candidates or []):
            if is_go:
                args.go = p
                print(f"[Pathfinder] GO auto-détecté (par contenu) : {p}")
                break
        # GO possiblement dans le même classeur que les stats
        if not args.go and args.stats and _peek(args.stats)[1]:
            args.go = args.stats
            print("[Pathfinder] GO détecté dans le classeur stats.")

    provenance = {"outil": "Pathfinder v1 (offline core)",
                  "date": time.strftime("%Y-%m-%d %H:%M:%S")}

    # ---- Design expérimental (échantillon→condition) ----
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
        else:                                   # design dans une sheet du classeur ?
            try:
                for nm, df in pd.read_excel(args.stats, sheet_name=None).items():
                    cols = [str(c).lower() for c in df.columns]
                    by_name = any(k in str(nm).lower() for k in ("design", "metadata", "experiment", "samples"))
                    # par contenu : une colonne 'condition/group…' + peu de lignes (≈ n échantillons)
                    by_content = (any(any(k in c for k in ("condition", "group", "groupe", "class", "treatment"))
                                      for c in cols) and 1 <= len(df) <= 500 and df.shape[1] <= 12)
                    if by_name or by_content:
                        design_map, design_src = io_layer.load_design(df), f"{os.path.basename(args.stats)}#{nm}"
                        design_frame = io_layer.load_design_frame(df)
                        break
            except Exception:
                pass
    if design_src:
        print(f"[Pathfinder] design détecté : {design_src}")

    # ---- Chargement ----
    stats = io_layer.load_stats(args.stats, design_map=design_map)
    if not stats.conditions_resolved:
        print("[Pathfinder] ⚠ conditions NON résolues : les colonnes échantillon "
              f"({', '.join(map(str, stats.sample_int_cols[:4]))}…) ne codent pas la condition et "
              "aucun design n'a été trouvé. Fournis --design ExperimentalDesign.csv "
              "(colonnes échantillon + condition) pour colorer la PCA et tester PC↔condition.")

    if args.write_design_template:
        tpl = os.path.join(args.outdir, "ExperimentalDesign_template.csv")
        pd.DataFrame({"sample": stats.sample_int_cols,
                      "condition": ["A_REMPLIR"] * len(stats.sample_int_cols)}).to_csv(tpl, index=False)
        print(f"[Pathfinder] template écrit : {tpl}\n"
              f"  → remplis la colonne 'condition' ({len(stats.sample_int_cols)} lignes), puis :\n"
              f"     python pathfinder.py --design {os.path.basename(tpl)}")
        return None

    if args.diagnose:
        from collections import Counter
        print("\n=== DIAGNOSTIC Pathfinder ===")
        print(f"Fichier stats : {args.stats}\n")
        # structure de TOUT le classeur
        allsheets = pd.read_excel(args.stats, sheet_name=None)
        print("Sheets du classeur :")
        for name, df in allsheets.items():
            num = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
            has_pval = any(str(c).endswith("_p.val") for c in df.columns)
            has_term = "term_name" in [str(c) for c in df.columns]
            # colonnes numériques hors motifs stats connus (candidates échantillons)
            statish = [c for c in df.columns if any(str(c).endswith(s) for s in ("_diff","_p.val","_p.adj"))
                       or str(c).startswith(("Robustness_Score_","Pi_Score_"))]
            sample_like = [c for c in num if c not in statish and str(c) not in io_layer._META]
            flags = []
            if has_pval: flags.append("STATS(_p.val)")
            if has_term: flags.append("GO(term_name)")
            print(f"  · {name:<28} lignes={len(df):>5}  cols={len(df.columns):>4}  "
                  f"numeriques={len(num):>4}  candidates_echantillon={len(sample_like):>4}  {' '.join(flags)}")
            if sample_like:
                ex = ', '.join(str(c) for c in sample_like[:8])
                print(f"      ex. colonnes échantillon : {ex}{' …' if len(sample_like)>8 else ''}")
        print()
        # ce que le loader a effectivement retenu
        sizes = Counter(stats.cond_map[c] for c in stats.sample_int_cols)
        print(f"[Loader] sheet comparaison retenue : {len(stats.comp)} lignes")
        print(f"[Loader] MATRICE PCA : sheet '{stats.pca_source}' → {len(stats.sample_int_cols)} colonnes échantillon")
        if stats.sample_int_cols:
            print("         ex. : " + ", ".join(str(c) for c in stats.sample_int_cols[:10]))
        print(f"[Loader] Z-score : {'présente' if stats.zscore is not None else 'absente'}")
        print(f"[Loader] contrastes ({len(stats.contrasts)}) : {', '.join(stats.contrasts[:6])}"
              f"{' …' if len(stats.contrasts)>6 else ''}")
        print(f"[Loader] colonnes intensité PCA retenues : {len(stats.sample_int_cols)}")
        print(f"[Loader] conditions ({len(sizes)}) tailles : "
              + ", ".join(f"{c}={n}" for c, n in sizes.items()))
        from collections import Counter
        sizes = Counter(stats.cond_map[c] for c in stats.sample_int_cols)
        print(f"[Loader] design : {design_src or 'non fourni'} → conditions "
              f"{'RÉSOLUES' if stats.conditions_resolved else 'NON résolues'}")
        print(f"[Loader] conditions ({len(sizes)}) tailles : "
              + ", ".join(f"{c}={n}" for c, n in list(sizes.items())[:12]))
        print("=============================\n")
        return None

    provenance["fichier_stats"] = os.path.basename(args.stats)
    provenance["contrastes"] = ", ".join(stats.contrasts)
    provenance["conditions"] = ", ".join(stats.conditions)
    provenance["design"] = design_src or "non fourni"
    provenance["conditions_resolues"] = "oui" if stats.conditions_resolved else "NON — fournir --design"
    provenance["n_echantillons"] = len(stats.sample_int_cols)
    provenance["n_proteines"] = len(stats.comp)

    go_res = None
    if args.go and os.path.exists(args.go):
        go = io_layer.load_go(args.go, stats.contrasts)
        provenance["fichier_go"] = os.path.basename(args.go)
        provenance["n_lignes_go"] = len(go.long)
        go_res = axis_meta_go.run(go.long, p_thr=args.pthr, jaccard_cut=args.jaccard)
    else:
        # matrice vide -> onglet Méta-GO informatif
        empty = pd.DataFrame(columns=["term_name", "term_id", "source", "n_sig",
                                      "redundancy_cluster", "is_representative"])
        go_res = {"matrix": pd.DataFrame(), "meta": empty,
                  "directionality": pd.DataFrame(), "term_genes": {}}
        provenance["fichier_go"] = "absent — axe méta-GO ignoré"

    # ---- Axe C ----
    # facteurs du design (multi-facteurs : condition, batch, ordre d'injection…)
    factors_df = None
    if design_frame is not None:
        try:
            factors_df = io_layer.map_design_factors(stats.sample_int_cols, design_frame[0], design_frame[1])
            if factors_df is not None:
                print(f"[Pathfinder] facteurs du design testés : {', '.join(map(str, factors_df.columns))}")
        except Exception as e:
            print(f"[Pathfinder] ⚠ facteurs du design non exploitables ({type(e).__name__})")

    explore = axis_explore.run(stats, do_umap=not args.no_umap, factors_df=factors_df)

    # ---- Espèce (taxid) + Axe A publications (online, opt-in) ----
    taxid, org, src = species_mod.resolve_species(org_arg=args.org)
    provenance["espece_taxid"] = f"{taxid} ({src})" if taxid else src

    pubs = None
    if args.online:
        if not taxid:
            print("[Pathfinder] ⚠ --online demandé mais taxid inconnu : précise --org "
                  "(hsapiens, mmusculus, rnorvegicus… ou un taxid numérique).")
        else:
            print(f"[Pathfinder] axe A (publications) : requêtes STRING PMID (taxid {taxid})…")
            pubs = axis_publications.run(stats, taxid, p_thr=args.pthr, refresh=args.refresh)
            n = sum(len(v) for v in pubs["results"].values())
            provenance["publications"] = f"STRING PMID · {len(pubs['results'])} contrastes · {n} PMID"
            for c, info in pubs["provenance"].items():
                provenance[f"pub[{c}]"] = info

    # ---- KPIs ----
    n_go_sig = int(go_res["matrix"].notna().sum().sum()) if not go_res["matrix"].empty else 0
    kpis = {
        "Contrastes": len(stats.contrasts),
        "Conditions": len(stats.conditions),
        "Échantillons": len(stats.sample_int_cols),
        "Protéines": len(stats.comp),
        "Cellules GO sig.": n_go_sig,
        "Var. PC1": f"{explore['pca']['var_exp'][0]:.0f}%",
    }

    # ---- Résumé interprétatif ----
    resume = summary_mod.generate(stats, go_res, explore, provenance)
    md_path = os.path.join(args.outdir, "pathfinder_resume.md")
    try:
        with open(md_path, "w", encoding="utf-8") as f:
            f.write(resume["markdown"])
    except PermissionError:
        import time as _t
        md_path = md_path.replace(".md", f"_{_t.strftime('%H%M%S')}.md")
        with open(md_path, "w", encoding="utf-8") as f:
            f.write(resume["markdown"])

    # ---- Sorties (robustes au verrou Windows si un fichier est ouvert) ----
    def _safe_write(writer_fn, path, label):
        try:
            writer_fn(path)
            return path
        except PermissionError:
            import time as _t
            alt = path.replace(".", f"_{_t.strftime('%H%M%S')}.", 1)
            print(f"[Pathfinder] ⚠ {os.path.basename(path)} verrouillé (ouvert dans un logiciel ?) "
                  f"→ écriture dans {os.path.basename(alt)}")
            writer_fn(alt)
            return alt

    xlsx_path = _safe_write(lambda p: build_xlsx.build(p, go_res, explore, provenance, resume, pubs),
                            os.path.join(args.outdir, "pathfinder_synthese.xlsx"), "xlsx")
    # branding (icône embarquée dans le header, style Deimos)
    import base64
    _here = os.path.dirname(os.path.abspath(__file__))
    icon_b64 = None
    _dirs = [_here, os.path.join(_here, "assets"), ".", os.path.join(".", "assets")]
    for _d in _dirs:
        for _nm in ("pathfinder_icon_128.png", "pathfinder_icon_256.png", "pathfinder_icon_200.png"):
            _p = os.path.join(_d, _nm)
            if os.path.exists(_p):
                icon_b64 = base64.b64encode(open(_p, "rb").read()).decode()
                break
        if icon_b64:
            break
    branding = {"icon_b64": icon_b64}

    html_path = _safe_write(lambda p: build_html.build(p, go_res, explore, provenance, kpis, resume, branding, pubs),
                            os.path.join(args.outdir, "pathfinder_report.html"), "html")

    print("\n" + resume["markdown"] + "\n")
    print(f"[Pathfinder] OK\n  - {xlsx_path}\n  - {html_path}\n  - {md_path}")
    return xlsx_path, html_path


if __name__ == "__main__":
    sys.exit(0 if main() else 0)
