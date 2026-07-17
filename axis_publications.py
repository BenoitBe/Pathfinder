"""
Publications axis (online).

For each contrast, takes the significant proteins and queries STRING's PMID
enrichment: "which papers are over-represented in this signature". Results are
cached (key = accessions+taxid+endpoint) with a fetch date (snapshot), because the
literature corpus evolves over time.

Opt-in (--online). Requires a taxid (resolved from --org, a gProfiler organism).
Known bias: heavily studied proteins are over-represented -> read alongside the
LFC/Pi score (input genes are ranked by |LFC| to attenuate it).
"""
from __future__ import annotations
import pandas as pd

import cache

STRING_URL = "https://string-db.org/api/tsv/enrichment"

# Enrichment is a statistical statement about a SET: below this many genes STRING
# merely returns papers citing them, which is not an enrichment. Such contrasts are
# reported as skipped rather than dressed up as a result.
MIN_GENES = 5


def _string_request(identifiers, taxid, timeout=30):
    """Live STRING call (POST). Raises if the network is unavailable -> axis skipped cleanly."""
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
            "title": col(r, "description"),
            "n_genes": col(r, "number_of_genes"),
            "genes": col(r, "preferredNames", "inputGenes"),
            "fdr": pd.to_numeric(col(r, "fdr"), errors="coerce"),
            "p_value": pd.to_numeric(col(r, "p_value"), errors="coerce"),
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
        })
    rows.sort(key=lambda d: (d["fdr"] if pd.notna(d["fdr"]) else 1.0))
    return rows[:top]


def _sig_identifiers(stats, contrast, p_thr, max_ids):
    """Genes (default) or accessions of a contrast's significant proteins, ranked by |LFC|.

    Returns (ids, n_total) so the caller can tell a small signature from a capped one.
    """
    comp = stats.comp
    padj, diff = f"{contrast}_p.adj", f"{contrast}_diff"
    if padj not in comp.columns:
        return [], 0
    sub = comp[pd.to_numeric(comp[padj], errors="coerce") < p_thr].copy()
    if diff in comp.columns:
        sub["_abs"] = pd.to_numeric(sub[diff], errors="coerce").abs()
        sub = sub.sort_values("_abs", ascending=False)
    idcol = stats.gene_col if (stats.gene_col in comp.columns) else stats.id_col
    ids = [str(x) for x in sub[idcol].dropna().tolist() if str(x) not in ("", "nan")]
    # de-duplicate while preserving order
    seen, out = set(), []
    for x in ids:
        if x not in seen:
            seen.add(x); out.append(x)
    n_total = len(out)
    return out[:max_ids], n_total


def run(stats, taxid, p_thr=0.05, max_ids=400, top_pmid=15, refresh=False, fetcher=None):
    fetch = fetcher or _string_request
    results, prov = {}, {}
    prov["_criterion"] = (f"significant proteins: p.adj < {p_thr} | max {max_ids} genes per query "
                          f"(ranked by |LFC|) | top {top_pmid} PMIDs kept per contrast")
    n_skipped = 0
    for c in stats.contrasts:
        ids, n_total = _sig_identifiers(stats, c, p_thr, max_ids)
        if n_total == 0:
            prov[c] = f"no protein with p.adj < {p_thr} - skipped"
            n_skipped += 1
            continue
        if n_total < MIN_GENES:
            prov[c] = (f"only {n_total} significant protein(s) (< {MIN_GENES}) - skipped: "
                       "enrichment is not meaningful on so few genes")
            n_skipped += 1
            continue
        capped = " (capped)" if n_total > max_ids else ""
        cached = cache.get(ids, taxid, "string_pmid", refresh=refresh)
        if cached:
            rows, when = cached["value"], cached.get("_fetch_date", "cache")
        else:
            try:
                tsv = fetch(ids, taxid)
            except Exception as e:
                prov[c] = f"network failure ({type(e).__name__})"
                continue
            rows = _parse_pmid(tsv, top_pmid)
            put = cache.put(ids, taxid, "string_pmid", rows)
            when = put["_fetch_date"]
        results[c] = rows
        trunc = f", truncated to top {top_pmid}" if len(rows) >= top_pmid else ""
        prov[c] = (f"{len(ids)}/{n_total} genes{capped} -> {len(rows)} PMIDs{trunc} "
                   f"(fetched {when})")
    prov["_skipped"] = f"{n_skipped} contrast(s) skipped (too few significant proteins)"
    return {"results": results, "provenance": prov, "taxid": taxid}


def to_dataframe(pub):
    recs = []
    for c, rows in pub["results"].items():
        for r in rows:
            recs.append({"contrast": c, **r})
    return pd.DataFrame(recs)
