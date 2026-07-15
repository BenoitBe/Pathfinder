"""
Axe A — Publications (online).

Pour chaque contraste, prend les protéines significatives et interroge
l'enrichissement PMID de STRING : « quels articles sont sur-représentés dans
cette signature ». Résultats mis en cache (clé = accessions+taxid+endpoint) avec
date de fetch (snapshot), car le corpus littérature évolue dans le temps.

Opt-in (--online). Nécessite un taxid (résolu depuis --org gprofiler).
Biais connu : sur-représentation des protéines très étudiées → à lire avec le
LFC/Pi en tête (on classe par FDR mais on trie les gènes d'entrée par |LFC|).
"""
from __future__ import annotations
import pandas as pd

import cache

STRING_URL = "https://string-db.org/api/tsv/enrichment"


def _string_request(identifiers, taxid, timeout=30):
    """Appel réel STRING (POST). Levé si réseau indisponible → axe A ignoré proprement."""
    import requests
    data = {"identifiers": "\r".join(map(str, identifiers)),
            "species": int(taxid), "caller_identity": "pathfinder"}
    r = requests.post(STRING_URL, data=data, timeout=timeout)
    r.raise_for_status()
    return r.text


def _parse_pmid(tsv, top=15):
    lines = [l for l in tsv.splitlines() if l.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    ix = {h: i for i, h in enumerate(header)}
    def col(row, *names):
        for n in names:
            if n in ix and ix[n] < len(row):
                return row[ix[n]]
        return ""
    rows = []
    for line in lines[1:]:
        r = line.split("\t")
        if col(r, "category") != "PMID":
            continue
        term = col(r, "term")
        pmid = term.split(":")[-1] if term else ""
        rows.append({
            "pmid": pmid,
            "titre": col(r, "description"),
            "n_genes": col(r, "number_of_genes"),
            "genes": col(r, "preferredNames", "inputGenes"),
            "fdr": pd.to_numeric(col(r, "fdr"), errors="coerce"),
            "p_value": pd.to_numeric(col(r, "p_value"), errors="coerce"),
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
        })
    rows.sort(key=lambda d: (d["fdr"] if pd.notna(d["fdr"]) else 1.0))
    return rows[:top]


def _sig_identifiers(stats, contrast, p_thr, max_ids):
    """Gènes (défaut) ou accessions des protéines sig d'un contraste, triés par |LFC|."""
    comp = stats.comp
    padj, diff = f"{contrast}_p.adj", f"{contrast}_diff"
    if padj not in comp.columns:
        return []
    sub = comp[pd.to_numeric(comp[padj], errors="coerce") < p_thr].copy()
    if diff in comp.columns:
        sub["_abs"] = pd.to_numeric(sub[diff], errors="coerce").abs()
        sub = sub.sort_values("_abs", ascending=False)
    idcol = stats.gene_col if (stats.gene_col in comp.columns) else stats.id_col
    ids = [str(x) for x in sub[idcol].dropna().tolist() if str(x) not in ("", "nan")]
    # dédoublonne en gardant l'ordre
    seen, out = set(), []
    for x in ids:
        if x not in seen:
            seen.add(x); out.append(x)
        if len(out) >= max_ids:
            break
    return out


def run(stats, taxid, p_thr=0.05, max_ids=400, top_pmid=15, refresh=False, fetcher=None):
    fetch = fetcher or _string_request
    results, prov = {}, {}
    for c in stats.contrasts:
        ids = _sig_identifiers(stats, c, p_thr, max_ids)
        if not ids:
            continue
        cached = cache.get(ids, taxid, "string_pmid", refresh=refresh)
        if cached:
            rows, when = cached["value"], cached.get("_fetch_date", "cache")
        else:
            try:
                tsv = fetch(ids, taxid)
            except Exception as e:
                prov[c] = f"échec réseau ({type(e).__name__})"
                continue
            rows = _parse_pmid(tsv, top_pmid)
            put = cache.put(ids, taxid, "string_pmid", rows)
            when = put["_fetch_date"]
        results[c] = rows
        prov[c] = f"{len(ids)} gènes → {len(rows)} PMID (fetch {when})"
    return {"results": results, "provenance": prov, "taxid": taxid}


def to_dataframe(pub):
    recs = []
    for c, rows in pub["results"].items():
        for r in rows:
            recs.append({"contraste": c, **r})
    return pd.DataFrame(recs)
