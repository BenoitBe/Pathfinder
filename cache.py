"""
Disk cache for online calls (Publications axis).
Key = sha1(sorted accessions) + taxid + endpoint. Stores the value + fetch date.
"""
from __future__ import annotations
import hashlib
import json
import os
import time
from pathlib import Path

CACHE_DIR = Path(os.environ.get("PATHFINDER_CACHE", Path.home() / ".pathfinder_cache"))


def _key(accessions, taxid, endpoint):
    h = hashlib.sha1((",".join(sorted(accessions))).encode()).hexdigest()[:16]
    return f"{endpoint}__{taxid}__{h}.json"


def get(accessions, taxid, endpoint, ttl_days=30, refresh=False):
    if refresh:
        return None
    p = CACHE_DIR / _key(accessions, taxid, endpoint)
    if not p.exists():
        return None
    try:
        obj = json.loads(p.read_text())
        if ttl_days and (time.time() - obj.get("_fetch_ts", 0)) > ttl_days * 86400:
            return None
        return obj
    except Exception:
        return None


def put(accessions, taxid, endpoint, value):
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    p = CACHE_DIR / _key(accessions, taxid, endpoint)
    obj = {"value": value, "_fetch_ts": time.time(),
           "_fetch_date": time.strftime("%Y-%m-%d %H:%M:%S"),
           "taxid": taxid, "endpoint": endpoint, "n_acc": len(accessions)}
    p.write_text(json.dumps(obj))
    return obj
