"""GSEA of Deimos (GSEA_* sheets) as a source of functional groups.

When the ORA is empty (non-model species, H. glaber case: few proteins pass the
selection rule, or g:Profiler / STRING return nothing), the pre-ranked GSEA of Deimos
still reads coordinated shifts over the whole ranking. Its pathways are turned into the
same long table as the GO sheets (contrast, source, term_id, term_name, p_value, z,
intersection), with:
  - p_value = FDR q-value of gseapy (the sheets hold the pathways at FDR <= gsea_fdr);
  - z = NES (sign = direction, like the GO z-score);
  - intersection = the leading-edge proteins (Lead_proteins column, Deimos >= 2026-11-07;
    older sheets: leading-edge genes matched to the Genes column).
Redundant pathways (same leading edge) are then grouped exactly like the Meta-GO terms.
"""
from __future__ import annotations

import re

import pandas as pd

_ID = re.compile(r"^((?:GO|KEGG|REAC|WP|HP|HPA|MIRNA|CORUM|TF):[^\s]+)\s+(.*)$")
_ID_TAIL = re.compile(r"^(.*?)\s*\(((?:GO|R-[A-Z]{3})[:\-][\w\-.]+)\)\s*$")


def term_parts(term):
    """'GO:BP__GO:0009414 response to water deprivation' -> (source, id, name);
    'GO_Biological_Process_2023__translation (GO:0006412)' -> (library, GO:0006412, name)."""
    t = str(term)
    src, rest = (t.split("__", 1) if "__" in t else ("", t))
    m = _ID.match(rest)
    if m:
        return src, m.group(1), m.group(2).strip()
    m = _ID_TAIL.match(rest)
    if m:
        return src, m.group(2), m.group(1).strip()
    return src, "", rest.strip()


def _sheet_contrast(sheet, contrasts, declared):
    if declared.get(sheet) in contrasts:
        return declared[sheet]
    for c in contrasts:
        if ("GSEA_" + c.replace("_vs_", "v").replace("_", ""))[:31] == sheet:
            return c
    return None


def read_sheet(xl, sheet):
    """(table, provenance): the gseapy table starts at its 'Term' column (the NES plot
    fills the left of the sheet), the provenance sits in A1."""
    raw = pd.read_excel(xl, sheet_name=sheet, header=None)
    if raw.empty:
        return pd.DataFrame(), ""
    head = [str(v).strip() if pd.notna(v) else "" for v in raw.iloc[0]]
    prov = head[0] if head and head[0] and head[0] != "Term" else ""
    if "Term" not in head:
        return pd.DataFrame(), prov
    j = head.index("Term")
    cols = [h for h in head[j:] if h]
    tab = raw.iloc[1:, j:j + len(cols)]
    tab.columns = cols
    return tab.dropna(subset=["Term"]).reset_index(drop=True), prov


def load_gsea(stats_path, contrasts, manifest, de, lookup):
    """Long table of the GSEA pathways (one row per contrast x pathway), or None.

    de: reading.io Run.de (index = protein id); lookup: token -> protein id."""
    xl = pd.ExcelFile(stats_path)
    sheets = [s for s in xl.sheet_names if str(s).startswith("GSEA_")]
    if not sheets:
        return None
    declared = (manifest or {}).get("gsea_sheet", {})
    gene2prot = {}
    if "Genes" in de.columns:
        for pid, g in zip(de.index, de["Genes"].fillna("").astype(str)):
            k = g.split(";")[0].strip().upper()
            if k:
                gene2prot.setdefault(k, []).append(pid)
    rows, provs = [], {}
    for sh in sheets:
        c = _sheet_contrast(sh, contrasts, declared)
        if c is None:
            continue
        tab, prov = read_sheet(xl, sh)
        if tab.empty or "NES" not in tab.columns:
            continue
        provs[c] = prov
        for _, r in tab.iterrows():
            src, tid, name = term_parts(r["Term"])
            if "Lead_proteins" in tab.columns and pd.notna(r.get("Lead_proteins")):
                acc = [lookup.get(a.strip()) for a in str(r["Lead_proteins"]).split(";")]
                lead = [a for a in acc if a]
            else:
                lead = []
                for g in str(r.get("Lead_genes") or "").split(";"):
                    lead += gene2prot.get(g.strip().upper(), [])
            rows.append({"contrast": c, "sheet": sh, "source": src or "GSEA",
                         "term_id": tid, "term_name": name,
                         "p_value": pd.to_numeric(r.get("FDR q-val"), errors="coerce"),
                         "z": pd.to_numeric(r.get("NES"), errors="coerce"),
                         "intersection": list(dict.fromkeys(lead))})
    if not rows:
        return None
    out = pd.DataFrame(rows)
    out = out.dropna(subset=["z"])
    out["p_value"] = out["p_value"].fillna(1.0)
    out.attrs["provenance"] = provs
    return out.reset_index(drop=True)
