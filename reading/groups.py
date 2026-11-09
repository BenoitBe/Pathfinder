"""Functional groups of the biological reading = the Meta-GO redundancy groups.

The significant GO / pathway terms (q < 0.05, all contrasts) that rest on the same
proteins are already grouped by the Meta-GO axis (Jaccard distance on the
intersections, average linkage, --jaccard cut, 0.5 by default). The reading takes
exactly these groups, named by the Meta-GO representative term ("Non-redundant"
view): a group carries the same name in every tab.

Returned object: the interface of InSight's former "themes" (names, genes(), top(),
group x contrast table), used by the programmes, the volcanos and the tables.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

MIN_GENES = 3          # a term linked to fewer than 3 proteins joins no group


@dataclass
class Groups:
    terms: pd.DataFrame          # term, name, source, p, genes (set), theme, cluster, rep
    hits: pd.DataFrame           # contrast, term, p, z, genes, theme
    names: dict                  # theme -> name (Meta-GO representative term)
    cluster_of: dict             # theme -> Meta-GO redundancy_cluster
    table: pd.DataFrame = field(default_factory=pd.DataFrame)  # theme x contrast
    n_unmapped: int = 0          # intersection accessions matching no protein
    origin: str = "ORA"          # 'ORA' (GO sheets) or 'GSEA' (leading edges)

    def top(self, n=16, min_terms=1) -> list:
        best = self.table.groupby("theme").p.min().sort_values()
        size = self.terms.groupby("theme").term.nunique()
        return [t for t in best.index if size.get(t, 0) >= min_terms][:n]

    def genes(self, theme) -> set:
        sub = self.terms[self.terms.theme == theme]
        return set().union(*sub.genes) if len(sub) else set()


def _key(df: pd.DataFrame) -> pd.Series:
    """Same term key as axis_meta_go (term_id, else name|source)."""
    return pd.Series(np.where(df["term_id"].astype(str).str.len() > 0,
                              df["term_id"].astype(str),
                              df["term_name"].astype(str) + "|" + df["source"].astype(str)),
                     index=df.index)


def build_groups(go_res: dict, go_long: pd.DataFrame, lookup: dict,
                 q_thr: float = 0.05, origin: str = "ORA") -> Groups | None:
    """go_res: output of axis_meta_go.run; go_long: GOData.long; lookup: token ->
    protein identifier (reading.io.id_lookup)."""
    if go_res is None or go_res["matrix"].empty or go_long is None or not len(go_long):
        return None
    meta = go_res["meta"]
    df = go_long[go_long["p_value"] < q_thr].copy()
    df["term"] = _key(df)
    df = df[df.term.isin(meta.index)]
    if df.empty:
        return None
    seen, unmapped = set(), set()

    def _map(lst):
        out = set()
        for a in lst:
            a = str(a).strip()
            pid = lookup.get(a) or lookup.get(a.split(";")[0].strip())
            if pid:
                out.add(pid)
            elif a:
                unmapped.add(a)
            seen.add(a)
        return out

    df["genes"] = df["intersection"].map(_map)
    hits = df[["contrast", "term", "p_value", "z", "genes"]].rename(columns={"p_value": "p"})
    hits["z"] = pd.to_numeric(hits["z"], errors="coerce").fillna(0.0)
    terms = (hits.groupby("term")
             .agg(p=("p", "min"), genes=("genes", lambda s: set().union(*s)))
             .reset_index())
    terms["name"] = terms.term.map(meta["term_name"]).astype(str)
    terms["source"] = terms.term.map(meta["source"]).astype(str)
    terms["cluster"] = terms.term.map(meta["redundancy_cluster"]).fillna(-1).astype(int)
    terms["rep"] = terms.term.map(meta["is_representative"]).fillna(True).astype(bool)
    terms = terms[terms.genes.map(len) >= MIN_GENES].reset_index(drop=True)
    if terms.empty:
        return None
    # a term without a group (cluster -1: fewer than two terms) forms its own group
    solo = terms.cluster < 0
    terms.loc[solo, "cluster"] = -(np.arange(int(solo.sum())) + 1) * 1000
    # groups numbered 0..n-1 in the order of their best term
    best = terms.groupby("cluster").p.min().sort_values(kind="mergesort")
    theme_of = {cl: i for i, cl in enumerate(best.index)}
    terms["theme"] = terms.cluster.map(theme_of).astype(int)
    names = {}
    for cl, t in theme_of.items():
        sub = terms[terms.cluster == cl]
        rep = sub[sub.rep]
        pick = rep if len(rep) else sub
        names[t] = pick.sort_values(["p", "name"]).name.iloc[0]
    hits = hits[hits.term.isin(set(terms.term))].copy()
    hits["theme"] = hits.term.map(dict(zip(terms.term, terms.theme))).astype(int)
    table = (hits.groupby(["theme", "contrast"])
             .agg(p=("p", "min"), z=("z", "median"), n=("term", "nunique"))
             .reset_index())
    return Groups(terms=terms, hits=hits.reset_index(drop=True), names=names,
                  cluster_of={t: int(cl) for cl, t in theme_of.items()}, table=table,
                  n_unmapped=len(unmapped), origin=origin)


def jaccard(a: set, b: set) -> float:
    u = len(a | b)
    return len(a & b) / u if u else 0.0
