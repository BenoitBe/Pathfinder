"""
Résolution de l'espèce pour l'axe A (STRING/UniProt exigent un taxid NCBI).
Échelle de résolution : argument explicite (chaînage Deimos) > champ Excel > re-demande.
"""
from __future__ import annotations

# org gprofiler -> taxid NCBI
GPROFILER_TO_TAXID = {
    "hsapiens": 9606, "mmusculus": 10090, "rnorvegicus": 10116, "drerio": 7955,
    "scerevisiae": 4932, "dmelanogaster": 7227, "celegans": 6239,
}


def resolve_species(org_arg=None, excel_org=None):
    """Retourne (taxid, org, source) ou (None, None, 'à demander')."""
    for org, src in ((org_arg, "argument (chaînage Deimos)"), (excel_org, "champ Excel")):
        if org:
            key = str(org).strip().lower()
            if key in GPROFILER_TO_TAXID:
                return GPROFILER_TO_TAXID[key], key, src
            if str(org).strip().isdigit():
                return int(org), None, src
    return None, None, "à demander (non renseigné / invalide)"
