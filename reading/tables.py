"""Tables of the biological reading: what the figures show, as reusable data, added to
the pathfinder_synthesis.xlsx workbook.

  Reading_facts       numeric facts of the Biological reading tab, section by section
  Programmes          one row per programme protein: programme, fit to the programme
                      profile (r, silhouette), functional groups, significant contrasts,
                      observed fraction, log2 means and z-scores per condition
  Programme_summary   one row per programme: size, high / low conditions, enriched
                      functional groups (q), example proteins
  Groups              one row per functional group (= Meta-GO redundancy group)
  Group_terms         one row per term: group, source, FDR, proteins
  Group_x_contrast    long format group x contrast (FDR of the best term, direction)
  Reading_columns     column dictionary

pathfinder_programmes.tsv: the Programmes sheet as tab-separated text (read directly in R).
"""
from __future__ import annotations

import html
import re

import numpy as np
import pandas as pd

TSV_NAME = "pathfinder_programmes.tsv"
_TAG = re.compile(r"<[^>]+>")


def _plain(s: str) -> str:
    """Text of a fact, without HTML tags."""
    return html.unescape(_TAG.sub("", str(s))).strip()


def _theme_id(t) -> str:
    return f"G{int(t) + 1}"


def _stab_class(j):
    from .programmes import stability_class
    return stability_class(j)


def _join(xs, sep="; ") -> str:
    return sep.join(str(x) for x in xs)


# --------------------------------------------------------------------------- programmes

def programmes_table(run, design, prog, themes=None) -> pd.DataFrame:
    """One row per protein assigned to a programme."""
    from sklearn.metrics import silhouette_samples
    Z = prog.z
    lab = prog.labels.reindex(Z.index)
    centroids = {p: Z[lab == p].mean() for p in range(1, prog.k + 1)}
    r = pd.Series({pid: float(np.corrcoef(Z.loc[pid].values, centroids[lab[pid]].values)[0, 1])
                   for pid in Z.index})
    sil = pd.Series(silhouette_samples(Z.values, lab.values), index=Z.index) \
        if prog.k > 1 and len(Z) > prog.k else pd.Series(np.nan, index=Z.index)

    # groups: the programme's enriched groups containing the protein, and every group
    th_genes = {t: themes.genes(t) for t in themes.names} if themes is not None else {}
    enriched = {p: [name for name, _q, _n in prog.themes.get(p, [])] for p in range(1, prog.k + 1)}
    name_to_t = {v: t for t, v in themes.names.items()} if themes is not None else {}

    sig_cols = {c: run.sig(c).fillna(False) for c in run.contrasts}
    rows = []
    for pid in Z.index:
        p = int(lab[pid])
        mine = [t for t, g in th_genes.items() if pid in g]
        in_prog = [n for n in enriched[p] if pid in th_genes.get(name_to_t.get(n), set())]
        sig_in = [run.short(c) for c in run.contrasts if bool(sig_cols[c].get(pid, False))]
        lfc = {c: run.de.at[pid, c + "_diff"] for c in run.contrasts if pid in run.de.index}
        peak_c = max(lfc, key=lambda c: abs(lfc[c]) if np.isfinite(lfc[c]) else -1) if lfc else None
        rows.append(dict(
            protein=pid, name=run.labels.get(pid, pid) if run.labels is not None else pid,
            programme=f"P{p}", r_to_programme=round(r[pid], 3),
            silhouette=round(float(sil[pid]), 3) if np.isfinite(sil[pid]) else np.nan,
            high_in=prog.high[p], low_in=prog.low[p],
            programme_groups=_join(in_prog),
            all_groups=_join(themes.names[t] for t in mine) if mine else "",
            n_significant=len(sig_in), significant_in=_join(sig_in),
            peak_contrast=run.short(peak_c) if peak_c else "",
            peak_log2FC=round(float(lfc[peak_c]), 3) if peak_c else np.nan))
    out = pd.DataFrame(rows).set_index("protein")
    if prog.protein_stability is not None:
        out.insert(out.columns.get_loc("silhouette") + 1, "stability",
                   prog.protein_stability.reindex(out.index).round(3))

    # observed fraction (DPC): a profile carried by non-detections is flagged
    if run.observed is not None and run.samples is not None:
        obs = run.observed.reindex(Z.index)
        cond = run.samples.set_index("label")["condition"]
        cols = [c for c in obs.columns if c in cond.index]
        frac = obs[cols].T.groupby(cond.loc[cols].values).mean().T
        frac = frac.reindex(columns=[c for c in Z.columns if c in frac.columns])
        out["observed_fraction"] = obs[cols].mean(axis=1).round(3)
        out["min_observed_fraction"] = frac.min(axis=1).round(3)
        out["least_observed_in"] = frac.idxmin(axis=1)
    conds = list(Z.columns)
    means = run.means.reindex(index=Z.index, columns=conds)
    out = out.join(means.round(3).add_prefix("log2_"))
    out = out.join(Z.round(3).add_prefix("z_"))
    order = {f"P{p}": p for p in range(1, prog.k + 1)}
    out = (out.assign(_o=out.programme.map(order))
              .sort_values(["_o", "r_to_programme"], ascending=[True, False])
              .drop(columns="_o"))
    return out.reset_index()


def programme_summary(prog, ptab: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for p in range(1, prog.k + 1):
        sub = ptab[ptab.programme == f"P{p}"]
        th = prog.themes.get(p, [])
        rows.append(dict(
            programme=f"P{p}", n_proteins=len(sub), high_in=prog.high[p], low_in=prog.low[p],
            mean_r_to_programme=round(float(sub.r_to_programme.mean()), 3),
            mean_silhouette=round(float(sub.silhouette.mean()), 3)
            if sub.silhouette.notna().any() else np.nan,
            stability_jaccard=round(prog.stability[p], 3) if prog.stability else np.nan,
            stability=_stab_class(prog.stability[p]) if prog.stability else "not computed",
            enriched_groups=_join(f"{n} (q={q:.2g}, {k} prot.)" for n, q, k in th) or "none",
            examples=_join(prog.names.get(p, []), ", ")))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- groups

def themes_tables(run, themes, prog=None, labels=None):
    """(Groups, Group_terms, Group_x_contrast)."""
    terms = themes.terms.copy()
    size = terms.groupby("theme").term.nunique()
    best = themes.table.groupby("theme").p.min()
    in_prog = {}
    if prog is not None:
        for p, lst in prog.themes.items():
            for n, q, _k in lst:
                in_prog.setdefault(n, []).append(f"P{p} (q={q:.2g})")
    lab = labels if labels is not None else {}
    rows = []
    for t in sorted(themes.names, key=lambda t: best.get(t, 1.0)):
        sub = terms[terms.theme == t].sort_values("p")
        genes = sorted(themes.genes(t))
        tc = themes.table[themes.table.theme == t].sort_values("p")
        rows.append(dict(
            group_id=_theme_id(t), group=themes.names[t],
            metago_cluster=themes.cluster_of.get(t), n_terms=int(size.get(t, 0)),
            n_proteins=len(genes), best_fdr=float(best.get(t, np.nan)),
            sources=_join(sorted(set(sub.source)), ", "),
            enriched_in=_join(f"{run.short(c)} (FDR={p:.2g})" for c, p in zip(tc.contrast, tc.p)),
            programmes=_join(in_prog.get(themes.names[t], [])),
            terms=_join(f"{i} {n}" for i, n in zip(sub.term, sub.name)),
            proteins=_join((lab.get(g, g) for g in genes), ", "),
            protein_ids=_join(genes, ", ")))
    tt = terms.assign(group_id=terms.theme.map(_theme_id),
                      group_name=terms.theme.map(themes.names),
                      n_proteins=terms.genes.map(len),
                      proteins=terms.genes.map(lambda g: _join(sorted(g), ", ")))
    tt = tt.sort_values(["theme", "p"])[["group_id", "group_name", "term", "name", "source",
                                          "rep", "p", "n_proteins", "proteins"]]
    tt = tt.rename(columns={"term": "term_id", "name": "term_name", "p": "best_fdr",
                            "rep": "metago_representative"})
    tx = themes.table.assign(group_id=themes.table.theme.map(_theme_id),
                             group_name=themes.table.theme.map(themes.names),
                             comparison=themes.table.contrast.map(run.short))
    tx = tx.sort_values(["theme", "p"])[["group_id", "group_name", "contrast", "comparison",
                                          "p", "z", "n"]]
    tx = tx.rename(columns={"p": "best_fdr", "z": "direction_z", "n": "n_terms"})
    return pd.DataFrame(rows), tt.reset_index(drop=True), tx.reset_index(drop=True)


# --------------------------------------------------------------------------- assembly

COLUMNS = [
    ("Programmes", "protein", "Deimos identifier (name, else head of the Protein.Group)"),
    ("Programmes", "name", "readable name (gene symbol, orthologue, description or accession)"),
    ("Programmes", "programme", "response programme (k-means on the z-scored condition means)"),
    ("Programmes", "r_to_programme", "Pearson correlation between the protein z profile and the "
                                     "mean profile of its programme"),
    ("Programmes", "silhouette", "silhouette of the protein (-1 to 1): < 0 = closer to another "
                                 "programme"),
    ("Programmes", "stability", "fraction of the replicate resamplings where the protein stays "
                                "in (the match of) its programme"),
    ("Programmes", "high_in / low_in", "conditions where the programme is highest / lowest"),
    ("Programme_summary", "stability_jaccard", "mean Jaccard between the programme and its match "
                                               "over the replicate resamplings (k-means re-run)"),
    ("Programme_summary", "stability", ">= 0.75 stable, 0.6-0.75 partial, < 0.6 unstable "
                                       "(Hennig's clusterboot thresholds)"),
    ("Programmes", "programme_groups", "functional groups enriched in the programme that contain "
                                       "the protein"),
    ("Programmes", "all_groups", "every group with a term containing the protein"),
    ("Programmes", "n_significant / significant_in", "comparisons where the protein passes the "
                                                     "selection rule (Proteins tab)"),
    ("Programmes", "peak_contrast / peak_log2FC", "comparison with the largest |log2FC|"),
    ("Programmes", "observed_fraction", "fraction of samples where the protein is detected (DPC)"),
    ("Programmes", "min_observed_fraction / least_observed_in",
     "lowest observed fraction over the conditions, and which one: a profile carried by "
     "non-detections (estimated values) shows here"),
    ("Programmes", "log2_<condition>", "log2 mean per condition (abundance matrix)"),
    ("Programmes", "z_<condition>", "z-score of the mean (profile used by the k-means)"),
    ("Programmes", "var_<term>", "share of the protein's variance carried by a design factor or "
                                 "two-factor interaction (type-II sums of squares), by the "
                                 "higher-order condition effects, or by the replicates"),
    ("Variance_partition", "<term>", "median share over the proteins of the set (all selected, "
                                     "then each programme)"),
    ("Groups", "group_id", "G1 = group with the best term; same groups as the Functional groups "
                           "figure"),
    ("Groups", "metago_cluster", "Meta-GO redundancy_cluster (MetaGO_matrix sheet)"),
    ("Groups", "group", "Meta-GO representative term (Non-redundant view)"),
    ("Groups", "best_fdr", "FDR of the best term of the group (q < 0.05, computed upstream)"),
    ("Groups", "enriched_in", "comparisons where at least one term of the group is enriched"),
    ("Groups", "programmes", "programmes where the group is over-represented (q < 0.05)"),
    ("Group_x_contrast", "direction_z", "median z-score of the terms (> 0: proteins rather up)"),
]


def reading_tables(run, design, prog=None, groups=None, sections=(), kpis=(), fact=None) -> dict:
    """{sheet name: DataFrame} for the synthesis workbook."""
    sheets = {}
    labels = run.labels.to_dict() if run.labels is not None else {}
    facts = [dict(section="Key figures", fact=f"{v} {n}") for n, v in kpis]
    for s in sections:
        facts += [dict(section=_plain(s["title"]), fact=_plain(b)) for b in s.get("bullets", [])]
    sheets["Reading_facts"] = pd.DataFrame(facts, columns=["section", "fact"])
    if prog is not None:
        ptab = programmes_table(run, design, prog, groups)
        if fact and fact.get("S") is not None:
            V = fact["S"].drop(columns="sum").round(3).add_prefix("var_")
            ptab = ptab.merge(V, left_on="protein", right_index=True, how="left")
        sheets["Programmes"] = ptab
        sheets["Programme_summary"] = programme_summary(prog, ptab)
    if groups is not None:
        g, gt, gx = themes_tables(run, groups, prog, labels)
        sheets["Groups"], sheets["Group_terms"], sheets["Group_x_contrast"] = g, gt, gx
    if fact and fact.get("summary") is not None:
        vp = fact["summary"].copy()
        vp.insert(1, "factors", ", ".join(fact["F"].columns) + f" ({fact['source']})")
        sheets["Variance_partition"] = vp
        sheets["Design_factors"] = fact["F"].rename_axis("condition").reset_index()
    sheets["Reading_columns"] = pd.DataFrame(COLUMNS, columns=["sheet", "column", "meaning"])
    return sheets
