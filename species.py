"""
Species resolution for the Publications axis (STRING/UniProt require an NCBI taxid).

Resolution order: explicit argument (Deimos chaining) > organism detected in the
workbook (Deimos records its gProfiler organism) > interactive prompt.
"""
from __future__ import annotations
import re
import sys

# gProfiler organism -> NCBI taxid
GPROFILER_TO_TAXID = {
    "hsapiens": 9606, "mmusculus": 10090, "rnorvegicus": 10116, "drerio": 7955,
    "scerevisiae": 4932, "dmelanogaster": 7227, "celegans": 6239, "athaliana": 3702,
    "btaurus": 9913, "sscrofa": 9823, "ggallus": 9031, "xtropicalis": 8364,
    "oaries": 9940, "ocuniculus": 9986, "cfamiliaris": 9615, "mmulatta": 9544,
    "bnapus": 3708, "osativa": 4530, "zmays": 4577, "slycopersicum": 4081,
    "ecoli": 511145, "dmelanogaster_dm6": 7227,
}

COMMON = {
    "human": "hsapiens", "mouse": "mmusculus", "rat": "rnorvegicus",
    "zebrafish": "drerio", "yeast": "scerevisiae", "fly": "dmelanogaster",
    "worm": "celegans", "arabidopsis": "athaliana", "cow": "btaurus",
    "pig": "sscrofa", "chicken": "ggallus", "dog": "cfamiliaris",
}


def normalise_org(value):
    """Accept 'hsapiens', 'Homo sapiens', 'human', or a numeric taxid -> (taxid, org)."""
    if value is None:
        return None, None
    v = str(value).strip()
    if not v:
        return None, None
    if v.isdigit():
        return int(v), None
    key = re.sub(r"[^a-z]", "", v.lower())
    if key in GPROFILER_TO_TAXID:
        return GPROFILER_TO_TAXID[key], key
    if key in COMMON:
        org = COMMON[key]
        return GPROFILER_TO_TAXID[org], org
    # 'Homo sapiens' / 'homo_sapiens' -> hsapiens
    parts = re.split(r"[\s_]+", v.strip().lower())
    if len(parts) >= 2 and parts[0] and parts[1]:
        cand = parts[0][0] + re.sub(r"[^a-z]", "", parts[1])
        if cand in GPROFILER_TO_TAXID:
            return GPROFILER_TO_TAXID[cand], cand
    return None, None


def detect_org_in_workbook(path):
    """Scan a workbook for the organism recorded upstream (e.g. Deimos 'go_organism')."""
    try:
        import pandas as pd
        sheets = pd.read_excel(path, sheet_name=None)
    except Exception:
        return None
    keys = ("go_organism", "organism", "species", "espece", "taxid", "gprofiler")
    for name, df in sheets.items():
        if df.empty:
            continue
        # key/value sheets (Methods) : look for the key in any column, value nearby
        for col in df.columns:
            ser = df[col].astype(str)
            hit = ser.str.lower().str.replace(r"[^a-z_]", "", regex=True).isin(keys)
            if hit.any():
                for idx in df.index[hit]:
                    row = df.loc[idx]
                    for v in row.values[::-1]:
                        taxid, org = normalise_org(v)
                        if taxid:
                            return org or str(taxid)
        # or a column literally named like the key
        for col in df.columns:
            if re.sub(r"[^a-z_]", "", str(col).lower()) in keys:
                for v in df[col].dropna().head(5):
                    taxid, org = normalise_org(v)
                    if taxid:
                        return org or str(taxid)
    return None


def prompt_species(default=None):
    """Ask the user for the species. Returns the raw answer, '' to skip, None if not a TTY."""
    if not (sys.stdin and sys.stdin.isatty()):
        return None
    examples = "hsapiens, mmusculus, rnorvegicus, drerio... or an NCBI taxid"
    if default:
        msg = (f"\n[Pathfinder] Publications axis (STRING) — species detected: {default}\n"
               f"             Press Enter to accept, type another ({examples}), "
               f"or 'n' to skip: ")
    else:
        msg = (f"\n[Pathfinder] Publications axis (STRING) — species? ({examples})\n"
               f"             Leave empty to skip: ")
    try:
        ans = input(msg).strip()
    except (EOFError, KeyboardInterrupt):
        return ""
    if ans.lower() in ("n", "no", "skip"):
        return ""
    return ans or (default or "")


def resolve_species(org_arg=None, excel_org=None):
    """Return (taxid, organism, source) or (None, None, reason)."""
    for org, src in ((org_arg, "argument"), (excel_org, "workbook")):
        if org:
            taxid, name = normalise_org(org)
            if taxid:
                return taxid, name, src
            return None, None, f"unrecognised organism '{org}' ({src})"
    return None, None, "not provided"
