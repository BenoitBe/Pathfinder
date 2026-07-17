"""Shared fixtures: synthetic workbooks reproducing the real-world shapes."""
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# Condition names whose tokens NEST ('HPS7' is a substring of 'EPHPS7'/'IPHPS7'):
# this is what broke the GO sheet mapping on real data.
NESTED_CONDITIONS = ["EPHPS7", "IPHPS7", "LPHPS7", "EPS7", "IPS7", "LPS7", "HPS7", "CS8"]


def contrasts_of(conditions):
    return [f"{a}_vs_{b}" for a, b in itertools.combinations(conditions, 2)]


def gene_symbols(n):
    """Real proteomes repeat gene symbols (isoforms, protein groups) and leave gaps:
    a unique-index assumption breaks on real data."""
    out = []
    for i in range(n):
        if i % 17 == 0:
            out.append("")                    # unannotated protein group
        else:
            out.append(f"G{i // 3:04d}")      # ~3 protein groups share one symbol
    return out


def _stats_frame(rng, conditions, n_proteins, contrasts, sig_per_contrast=None):
    comp = pd.DataFrame({
        "Protein.Group": [f"P{i:04d}" for i in range(n_proteins)],
        "Genes": gene_symbols(n_proteins),
    })
    for c in contrasts:
        p = rng.uniform(0, 1, n_proteins)
        if sig_per_contrast is not None:
            k = sig_per_contrast.get(c, 0)
            p = np.ones(n_proteins)
            p[:k] = 1e-6
        comp[f"{c}_diff"] = rng.normal(0, 1.5, n_proteins)
        comp[f"{c}_p.val"] = p
        comp[f"{c}_p.adj"] = p
    comp["significant"] = True
    return comp


@pytest.fixture
def workbook_full(tmp_path):
    """Workbook with THREE competing per-sample matrices, like the real one:
    raw LFQ (not log2, extra runs), full log2 (the right one), z-scored subset."""
    rng = np.random.default_rng(5)
    conds = NESTED_CONDITIONS
    contrasts = contrasts_of(conds)
    samples = [f"Sple-{i + 37}" for i in range(24)]
    n = 400

    eff = np.zeros((n, 8))
    for j in range(1, 8):
        eff[:120, j] = rng.normal(0, 1.5, 120)
    base = rng.normal(24, 4, (n, 1))              # real dynamic range across proteins
    X = base + eff[:, np.repeat(np.arange(8), 3)] + rng.normal(0, 0.4, (n, 24))

    log2 = pd.DataFrame(X, columns=samples)
    log2.insert(0, "Protein.Group", [f"P{i:04d}" for i in range(n)])
    log2.insert(1, "Genes", gene_symbols(n))

    idx = np.argsort(-X.var(axis=1))[:120]        # SELECTED subset, then z-scored
    Z = (X[idx] - X[idx].mean(1, keepdims=True)) / (X[idx].std(1, keepdims=True) + 1e-9)
    z = pd.DataFrame(Z, columns=samples)
    z.insert(0, "Protein.Group", [f"P{i:04d}" for i in idx])
    z.insert(1, "Genes", [gene_symbols(n)[i] for i in idx])

    raw = pd.DataFrame(2 ** X, columns=[f"LFQ.intensity.{s}" for s in samples])
    raw.insert(0, "Protein.Group", [f"P{i:04d}" for i in range(n)])

    comp = _stats_frame(rng, conds, n, contrasts)

    path = tmp_path / "workbook_full.xlsx"
    with pd.ExcelWriter(path) as xw:
        raw.to_excel(xw, sheet_name="raw_data", index=False)
        log2.to_excel(xw, sheet_name="Log2_Impute", index=False)
        comp.to_excel(xw, sheet_name="Differential_Expression", index=False)
        z.to_excel(xw, sheet_name="Zscore_Heatmap", index=False)

    design = tmp_path / "ExperimentalDesign.csv"
    pd.DataFrame({"sample": samples,
                  "condition": [c for c in conds for _ in range(3)],
                  "replicate": [i % 3 + 1 for i in range(24)]}).to_csv(design, index=False)
    return {"path": str(path), "design": str(design), "samples": samples,
            "conditions": conds, "contrasts": contrasts}


@pytest.fixture
def workbook_zscore_only(tmp_path):
    """Only a z-scored matrix is available (the shape that yields inflated effect sizes)."""
    rng = np.random.default_rng(9)
    conds = ["A", "B", "C"]
    contrasts = contrasts_of(conds)
    samples = [f"Sple-{i + 1}" for i in range(9)]
    n = 200
    X = rng.normal(22, 3, (n, 1)) + rng.normal(0, 0.5, (n, 9))
    Z = (X - X.mean(1, keepdims=True)) / (X.std(1, keepdims=True) + 1e-9)
    z = pd.DataFrame(Z, columns=samples)
    z.insert(0, "Protein.Group", [f"P{i:04d}" for i in range(n)])
    z.insert(1, "Genes", gene_symbols(n))
    comp = _stats_frame(rng, conds, n, contrasts)
    path = tmp_path / "zscore_only.xlsx"
    with pd.ExcelWriter(path) as xw:
        comp.to_excel(xw, sheet_name="Comparaison", index=False)
        z.to_excel(xw, sheet_name="Zscore_Heatmap", index=False)
    design = tmp_path / "design_z.csv"
    pd.DataFrame({"sample": samples,
                  "condition": [c for c in conds for _ in range(3)]}).to_csv(design, index=False)
    return {"path": str(path), "design": str(design), "conditions": conds,
            "contrasts": contrasts}


@pytest.fixture
def go_workbook(tmp_path, workbook_full):
    """GO sheets named 'GO_<A>v<B>' with a DISTINCT term count per contrast, so a
    mis-assignment is detectable by counting."""
    rng = np.random.default_rng(11)
    conds = workbook_full["conditions"]
    path = tmp_path / "go_only.xlsx"
    expected = {}
    with pd.ExcelWriter(path) as xw:
        for k, (a, b) in enumerate(itertools.combinations(conds, 2)):
            nterms = 5 + k                      # unique count -> traceable
            expected[f"{a}_vs_{b}"] = nterms
            pd.DataFrame({
                "source": ["GO:BP"] * nterms,
                "term_id": [f"GO:{9000000 + k * 100 + i}" for i in range(nterms)],
                "term_name": [f"term_{a}_{b}_{i}" for i in range(nterms)],
                "p_value": rng.uniform(1e-6, 0.01, nterms),
                "intersection": [",".join(f"P{rng.integers(0, 400):04d}" for _ in range(5))
                                 for _ in range(nterms)],
                "z_score": rng.normal(0, 1.5, nterms),
            }).to_excel(xw, sheet_name=f"GO_{a}v{b}"[:31], index=False)
    return {"path": str(path), "expected_terms": expected}
