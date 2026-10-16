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
    """Accept 'hsapiens', 'Homo sapiens', 'human', a numeric taxid, or a labelled taxid
    ('NCBI taxid 1304', 'taxid:1304', as Deimos writes it in STRING mode) -> (taxid, org).

    The labelled form must be the whole value: a taxid quoted inside free text (e.g. a
    strain in a Methods detail) is not the organism the analysis ran on."""
    if value is None:
        return None, None
    v = str(value).strip()
    if not v:
        return None, None
    if v.isdigit():
        return int(v), None
    m = re.fullmatch(r"(?:ncbi[\s_-]*)?tax(?:on)?[\s_-]*id[\W_]{0,3}(\d+)", v, re.I)
    if m:
        return int(m.group(1)), None
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
    # 1) the upstream manifest states the organism the analysis actually ran on
    #    (STRING taxid in STRING mode, g:Profiler code otherwise): it wins.
    org = sheets.get("_manifest")
    if org is not None and {"section", "key", "value"} <= set(map(str, org.columns)):
        rows = org[org["section"].astype(str) == "organism"]
        o = dict(zip(rows["key"].astype(str), rows["value"]))
        for k in ("taxid", "gprofiler"):
            v = o.get(k)
            if v is not None and pd.notna(v):
                if isinstance(v, float) and v.is_integer():
                    v = int(v)
                taxid, name = normalise_org(str(v).strip())
                if taxid:
                    return name or str(taxid)

    # 2) otherwise scan key/value sheets (Methods...) for an organism entry
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
    examples = "hsapiens, 'Brassica napus', 'camelina'... or an NCBI taxid"
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


STRING_SPECIES_URL = "https://stringdb-downloads.org/download/species.v12.0.txt"


def _species_cache_path():
    import os
    from pathlib import Path
    d = Path(os.environ.get("PATHFINDER_CACHE", Path.home() / ".pathfinder_cache"))
    return d / "string_species.tsv"


def fetch_string_species(timeout=30):
    """STRING's own species table (cached on disk). Returns [(taxid, name), ...].

    Resolving against STRING rather than a hard-coded list has two virtues: it covers
    every species STRING knows, and it proves upfront that the Publications axis can
    actually run for that organism.
    """
    p = _species_cache_path()
    if p.exists():
        try:
            text = p.read_text(encoding="utf-8")
        except Exception:
            text = ""
    else:
        text = ""
    if not text:
        import requests
        r = requests.get(STRING_SPECIES_URL, timeout=timeout)
        r.raise_for_status()
        text = r.text
        try:
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(text, encoding="utf-8")
        except Exception:
            pass
    out = []
    for line in text.splitlines()[1:]:
        f = line.split("\t")
        if len(f) >= 4 and f[0].strip().isdigit():
            out.append((int(f[0]), f[3].strip()))       # official NCBI name
    return out


def lookup_species(query, species=None):
    """Find candidate species for a free-text query. Returns [(taxid, name), ...].

    Accepts 'Cannabis sativa', 'cannabis', or the gProfiler-style 'csativa' — which is
    ambiguous on purpose (Cannabis / Camelina / Castanea sativa), so every candidate is
    returned and the caller asks rather than guessing.
    """
    if species is None:
        species = fetch_string_species()
    q = str(query).strip().lower()
    nq = re.sub(r"[^a-z]", "", q)
    exact, partial = [], []
    for taxid, name in species:
        nl = name.lower()
        nn = re.sub(r"[^a-z]", "", nl)
        if nq == nn or q == nl:
            exact.append((taxid, name))
            continue
        parts = nl.split()
        # gProfiler style: first letter of the genus + full species epithet
        if len(parts) >= 2 and nq == (parts[0][0] + re.sub(r"[^a-z]", "", parts[1])):
            partial.append((taxid, name))
        elif nq and nq in nn:
            partial.append((taxid, name))
    return exact or partial


def choose_candidate(candidates):
    """Pick among species matching an ambiguous query. Returns a taxid or None.

    'csativa' fits Cannabis, Camelina and Castanea sativa: ask rather than guess.
    """
    if not candidates:
        return None
    if len(candidates) == 1:
        taxid, name = candidates[0]
        print(f"[Pathfinder] species resolved: {name} (taxid {taxid})")
        return taxid
    if not (sys.stdin and sys.stdin.isatty()):
        print(f"[Pathfinder] \u26a0 ambiguous species ({len(candidates)} matches) and no "
              "terminal to ask: pass an NCBI taxid with --org.")
        return None
    print("\n[Pathfinder] Several species match - which one?")
    for i, (taxid, name) in enumerate(candidates[:10], 1):
        print(f"    {i}. {name}  (taxid {taxid})")
    try:
        ans = input("             Number, or Enter to skip: ").strip()
    except (EOFError, KeyboardInterrupt):
        return None
    if ans.isdigit() and 1 <= int(ans) <= min(len(candidates), 10):
        taxid, name = candidates[int(ans) - 1]
        print(f"[Pathfinder] species resolved: {name} (taxid {taxid})")
        return taxid
    return None


def resolve_species(org_arg=None, excel_org=None):
    """Return (taxid, organism, source) or (None, None, reason)."""
    for org, src in ((org_arg, "argument"), (excel_org, "workbook")):
        if org:
            taxid, name = normalise_org(org)
            if taxid:
                return taxid, name, src
            return None, None, f"unrecognised organism '{org}' ({src})"
    return None, None, "not provided"
