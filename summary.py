"""
Résumé interprétatif automatique (déterministe, offline).

Lit les objets déjà calculés (aucun recalcul de stat) et produit une synthèse
en français : capture du design par la PCA, suspicion de batch, échantillons
aberrants, synthèse fonctionnelle méta-GO. Sortie exploitable en Markdown,
en lignes (Excel) et en HTML.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

# seuils de taille d'effet (epsilon^2 / eta^2, conventions usuelles)
EFF_LARGE, EFF_MED = 0.14, 0.06


def _fmt_p(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "NA"
    return f"{p:.1e}" if p < 1e-3 else f"{p:.3f}"


def _design_section(explore):
    pca = explore["pca"]; assoc = explore["pc_factor"].set_index("PC")
    ve = {f"PC{i+1}": float(v) for i, v in enumerate(pca["var_exp"])}
    scores = pca["scores"]
    sizes = scores["condition"].value_counts()
    k, n, min_n = len(sizes), len(scores), int(sizes.min())
    lines = []

    # design sans réplicats suffisants -> test non calculable, pas d'alerte batch
    if n <= k or min_n < 2 or assoc["p_value"].isna().all():
        lines.append(f"{k} groupes pour {n} échantillons (min {min_n}/groupe) : "
                     f"test PC↔condition non calculable faute de réplicats — "
                     f"interpréter la PCA visuellement (onglet Exploration).")
        top = sorted(ve.items(), key=lambda kv: -kv[1])[:2]
        lines.append("Variance concentrée sur " +
                     ", ".join(f"{p} ({v:.0f}%)" for p, v in top) +
                     (" — un PC très dominant peut trahir un échantillon atypique ou un gradient technique."
                      if top and top[0][1] >= 80 else "."))
        return lines

    captured, batch = [], []
    pcol = "p_adj_BH" if "p_adj_BH" in assoc.columns else "p_value"
    for pc in assoc.index:
        eps2, p, var = assoc.loc[pc, "epsilon2"], assoc.loc[pc, pcol], ve.get(pc, 0)
        if pd.notna(p) and p < 0.05 and eps2 >= EFF_MED:
            captured.append((pc, var, eps2, p))
        elif var >= 10 and (pd.isna(p) or p > 0.2):
            batch.append((pc, var, eps2, p))

    if captured:
        parts = [f"{pc} ({var:.0f}% var., ε²={eps2:.2f}, p={_fmt_p(p)})" for pc, var, eps2, p in captured]
        strong = any(e >= EFF_LARGE for _, _, e, _ in captured)
        lines.append(f"La condition structure la variance sur {', '.join(parts)} — "
                     f"séparation des groupes {'nette' if strong else 'modérée'}.")
    else:
        lines.append("Aucune PC n'est significativement associée à la condition : "
                     "séparation des groupes faible ou variance dominée par d'autres facteurs.")
    for pc, var, eps2, p in batch:
        lines.append(f"⚠ {pc} porte {var:.0f}% de variance sans lien avec la condition "
                     f"(p={_fmt_p(p)}) — suspecter un effet batch / covariable technique.")
    return lines


def _multifactor_section(explore):
    multi = explore.get("pc_factor_multi")
    if multi is None or multi.empty:
        return []
    ve = {f"PC{i+1}": float(v) for i, v in enumerate(explore["pca"]["var_exp"])}
    TECH = ("batch", "lot", "plate", "plaque", "ordre", "order", "injection", "run",
            "date", "operateur", "operator", "colonne", "column")
    lines = []
    hits = multi[(multi["p_adj_BH"] < 0.1) & (multi["effet"] >= 0.14)]
    if hits.empty:
        return ["Aucun facteur du design n'est significativement associé à une PC (après correction BH)."]
    for fac, grp in hits.groupby("facteur"):
        g = grp.sort_values("effet", ascending=False)
        parts = [f"{r.PC} ({ve.get(r.PC, 0):.0f}% var., effet={r.effet:.2f})" for _, r in g.head(3).iterrows()]
        tech = any(t in str(fac).lower() for t in TECH)
        prefix = "⚠ Facteur technique" if tech else "Facteur"
        tail = (" — variance d'origine technique : envisager une correction (limma::removeBatchEffect, "
                "ComBat) ou l'inclure comme covariable." if tech else "")
        lines.append(f"{prefix} « {fac} » associé à {', '.join(parts)}{tail}")
    return lines
def _qc_section(explore, cond_map):
    corr = explore["correlation"]
    samples = list(corr.columns)
    lines, outliers = [], []
    for s in samples:
        peers = [o for o in samples if o != s and cond_map.get(o) == cond_map.get(s)]
        if len(peers) >= 2:
            r = float(np.mean([corr.loc[s, o] for o in peers]))
            if r < 0.90:
                outliers.append((s, cond_map.get(s), r))
    if outliers:
        for s, c, r in sorted(outliers, key=lambda t: t[2]):
            lines.append(f"⚠ Échantillon {s} ({c}) faiblement corrélé à ses réplicats (r̄={r:.2f}) "
                         "— vérifier run/préparation avant interprétation.")
    else:
        lines.append("Corrélation intra-condition homogène, pas d'échantillon aberrant détecté.")
    return lines


SRC_ORDER = ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"]
SRC_LABEL = {"GO:BP": "GO:BP (processus biologique)", "GO:MF": "GO:MF (fonction moléculaire)",
             "GO:CC": "GO:CC (composant cellulaire)", "KEGG": "KEGG", "REAC": "Reactome"}


def _sources_present(meta):
    present = [s for s in SRC_ORDER if (meta["source"] == s).any()]
    others = [s for s in meta["source"].dropna().unique() if s not in SRC_ORDER]
    return present + sorted(others)


def _go_section(meta_go):
    matrix, meta, direction = meta_go["matrix"], meta_go["meta"], meta_go["directionality"]
    if matrix.empty:
        return ["Enrichissement GO non fourni — synthèse fonctionnelle non disponible."]
    lines = [f"{matrix.shape[0]} termes significatifs sur {matrix.shape[1]} contraste(s), ventilés par source :"]
    for src in _sources_present(meta):
        keys = meta.index[meta["source"] == src]
        d = direction[direction["source"] == src] if "source" in direction.columns else direction.loc[keys]
        vc = d["pattern"].value_counts().to_dict()
        shared = vc.get("up partagé", 0) + vc.get("down partagé", 0)
        mixte = vc.get("mixte", 0)
        rec = meta.loc[keys].sort_values("n_sig", ascending=False)
        if len(rec):
            top = rec.iloc[0]
            dom = f" — dominant : {top['term_name']} ({int(top['n_sig'])} contrastes)"
        else:
            dom = ""
        lines.append(f"{SRC_LABEL.get(src, src)} : {len(keys)} termes "
                     f"({shared} partagés, {mixte} à direction inverse){dom}.")
    # mixtes globaux notables (régulation inverse)
    mx = direction[direction["pattern"] == "mixte"].head(3)
    if len(mx):
        lines.append("Régulation inverse selon le contraste (à investiguer) : "
                     + "; ".join(f"{r.term_name} [{r.source}]" for _, r in mx.iterrows()) + ".")
    if (meta["redundancy_cluster"] >= 0).any():
        nc = meta["redundancy_cluster"].nunique()
        lines.append(f"Après regroupement de redondance : {nc} groupes fonctionnels.")
    return lines


def _contrast_section(meta_go):
    matrix = meta_go["matrix"]
    if matrix.empty:
        return []
    counts = matrix.notna().sum(axis=0).sort_values(ascending=False)
    return [f"{c} : {int(n)} termes enrichis significatifs." for c, n in counts.items()]


def generate(stats, meta_go, explore, provenance):
    sections = [
        ("Design expérimental (PCA)", _design_section(explore)),
        ("Facteurs du design (PC × facteurs)", _multifactor_section(explore)),
        ("Qualité & homogénéité", _qc_section(explore, stats.cond_map)),
        ("Synthèse fonctionnelle (méta-GO)", _go_section(meta_go)),
        ("Détail par contraste", _contrast_section(meta_go)),
    ]
    sections = [(t, l) for t, l in sections if l]

    md = [f"# Résumé automatique — Pathfinder",
          f"_{provenance.get('date','')} · {provenance.get('contrastes','')}_", ""]
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
