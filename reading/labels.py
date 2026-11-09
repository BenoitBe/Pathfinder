"""Readable protein names: the gene symbol when informative, else the orthologue name
(Perseverance + reference FASTA, or STRING), else a short description, else the
accession."""
from __future__ import annotations

import glob
import os
import re

import pandas as pd

# identifiers that say nothing: locus tags, NCBI LOC, accessions
_UNIPROT = r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
_UNINFORMATIVE = re.compile(r"^(LOC\d+|[A-Z]{1,6}\d*_\d{3,}|" + _UNIPROT + r"|"
                            r"[A-Z]{2}_\d+(\.\d+)?|At\dg\d+|nan|None|)$", re.I)


def informative(name) -> bool:
    return isinstance(name, str) and not _UNINFORMATIVE.match(name.strip())


def _fasta_symbols(path: str, wanted: set) -> dict:
    out = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            tok = line[1:].split()[0]
            acc = tok.split("|")[1] if tok.count("|") >= 2 else tok
            if acc in wanted:
                m = re.search(r"GN=(\S+)", line)
                out[acc] = m.group(1) if m else None
    return out


def _reference_fasta_candidates(run, explicit):
    if explicit:
        return [explicit]
    cfg = run.config or {}
    c = [cfg.get(k) for k in ("perseverance_fasta_a", "perseverance_fasta_b") if cfg.get(k)]
    dirs = {run.out_dir, os.path.dirname(os.path.abspath(run.out_dir))}
    if cfg.get("perseverance_fasta_dir"):
        dirs.add(str(cfg["perseverance_fasta_dir"]))
    if cfg.get("tsv_path"):
        dirs.add(os.path.dirname(os.path.abspath(str(cfg["tsv_path"]))) or ".")
    for d in dirs:
        for ext in ("*.fasta", "*.fa", "*.faa"):
            c += glob.glob(os.path.join(d, ext))
    return [p for p in dict.fromkeys(c) if p and os.path.exists(p)
            and os.path.getsize(p) < 400e6]


def best_labels(run, xl, reference_fasta=None):
    de = run.de
    genes = (de["Genes"].astype(str).str.split(";").str[0].str.strip()
             if "Genes" in de.columns else pd.Series("", index=de.index))
    desc_col = next((c for c in ("First.Protein.Description", "Description")
                     if c in de.columns), None)
    desc = de[desc_col].fillna("").astype(str) if desc_col else pd.Series("", index=de.index)

    orth_names, src = {}, ""
    for sheet in ("Orthologues_Perseverance", "Orthologues_STRING"):
        if sheet not in xl.sheet_names:
            continue
        o = pd.read_excel(xl, sheet_name=sheet)
        pcol = next((c for c in ("protein_id", "query_id", "protein") if c in o.columns), None)
        name_col = next((c for c in ("preferredName", "preferred_name", "reference_name",
                                     "ortholog_name", "gene_name") if c in o.columns), None)
        if pcol is None:
            continue
        if name_col:
            orth_names = dict(zip(o[pcol].astype(str), o[name_col]))
            src = f"orthologue names ({sheet})"
        elif "ortholog_id" in o.columns:
            ids = dict(zip(o[pcol].astype(str), o["ortholog_id"].astype(str)))
            wanted = set(ids.values())
            best, best_n = {}, 0
            for f in _reference_fasta_candidates(run, reference_fasta):
                try:
                    sym = _fasta_symbols(f, wanted)
                except OSError:
                    continue
                n = sum(1 for v in sym.values() if v)
                if n > best_n:
                    best, best_n, src = sym, n, f"orthologues ({os.path.basename(f)})"
            orth_names = {p: best.get(a) for p, a in ids.items() if best.get(a)}
        if orth_names:
            break

    labels, n_gene, n_orth = {}, 0, 0
    for pid in de.index:
        g = genes.get(pid, "")
        if informative(g):
            labels[pid] = g; n_gene += 1; continue
        o = orth_names.get(pid)
        if informative(o):
            labels[pid] = o; n_orth += 1; continue
        d = re.sub(r"\s*(OS=|\[).*$", "", desc.get(pid, "")).strip()
        d = re.sub(r"^(PREDICTED:|uncharacterized protein)\s*", "", d, flags=re.I)
        labels[pid] = (d[:28] + "…") if len(d) > 29 else (d or pid)
    parts = []
    if n_gene:
        parts.append(f"{n_gene} gene symbols")
    if n_orth:
        parts.append(f"{n_orth} from {src}")
    return pd.Series(labels), ", ".join(parts) or "descriptions / accessions"
