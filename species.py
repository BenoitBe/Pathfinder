"""
Species resolution for the Publications axis (STRING/UniProt require an NCBI taxid).
Resolution order: explicit argument (Deimos chaining) > Excel field > prompt.
"""
from __future__ import annotations

# gProfiler organism -> NCBI taxid
GPROFILER_TO_TAXID = {
    "hsapiens": 9606, "mmusculus": 10090, "rnorvegicus": 10116, "drerio": 7955,
    "scerevisiae": 4932, "dmelanogaster": 7227, "celegans": 6239,
}


def resolve_species(org_arg=None, excel_org=None):
    """Return (taxid, organism, source) or (None, None, 'to be provided')."""
    for org, src in ((org_arg, "argument (Deimos chaining)"), (excel_org, "Excel field")):
        if org:
            key = str(org).strip().lower()
            if key in GPROFILER_TO_TAXID:
                return GPROFILER_TO_TAXID[key], key, src
            if str(org).strip().isdigit():
                return int(org), None, src
    return None, None, "to be provided (missing / invalid)"
