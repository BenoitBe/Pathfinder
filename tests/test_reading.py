"""Biological reading (ex-InSight, merged in 1.3.0): end to end on varied designs
(categorical, many levels, time course, crossed design, two conditions, minimal
workbook), same proteins as the Proteins tab, groups = Meta-GO groups, tables of the
synthesis workbook, valid and interactive SVG."""
import glob
import os
import re
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import pytest

import axis_meta_go
import axis_proteins
import io_layer
import pathfinder
from reading.factors import detect_design
from reading.groups import build_groups
from reading.io import id_lookup, load_run
from reading.programmes import build_programmes
from reading.run import run_reading
from reading.svg import fig_programmes
from reading.tables import TSV_NAME
from reading_fixture import make_deimos_run


def _xlsx(root):
    return glob.glob(os.path.join(root, "ProteomicAnalysis_Results*.xlsx"))[0]


def _load(root, pthr=None, ratio=None):
    x = _xlsx(root)
    stats = io_layer.load_stats(x)
    man = io_layer.read_manifest(x)
    rule = io_layer.selection_rule(man, p_thr=pthr, min_ratio=ratio)
    go_res, go_long = None, None
    try:
        go = io_layer.load_go(x, stats.contrasts, conditions=stats.conditions, verbose=False)
        if len(go.long):
            go_long = go.long
            go_res = axis_meta_go.run(go.long, jaccard_cut=0.5)
    except ValueError:
        pass
    return x, stats, man, rule, go_res, go_long


def _reading(root, **kw):
    x, stats, man, rule, go_res, go_long = _load(root)
    return run_reading(stats, x, rule, man, go_res, go_long, **kw)


def _prog(root):
    x, stats, man, rule, go_res, go_long = _load(root)
    run = load_run(stats, x, rule, man, verbose=False)
    d = detect_design(run.conditions, run.samples)
    groups = build_groups(go_res, go_long, id_lookup(run.de)) if go_res else None
    return run, d, build_programmes(run, d, groups)


def _main(root, tmp_path, *extra, monkeypatch=None):
    out = tmp_path / "out"
    pathfinder.main(["--stats", _xlsx(root), "--outdir", str(out), "--no-online", "--no-umap",
                     *extra])
    return out


def _ids(res):
    return [s["id"] for s in res["sections"]]


def _svgs(html):
    return re.findall(r"<svg class='rd'.*?</svg>", html, re.S)


@pytest.fixture(autouse=True)
def _isolated_cwd(tmp_path, monkeypatch):
    """Pathfinder looks for a design file in the current folder: isolate it."""
    monkeypatch.chdir(tmp_path)


def test_categorical_project_end_to_end(tmp_path, capsys):
    root = make_deimos_run(str(tmp_path / "run"), ["Ctrl", "DrugA", "DrugB", "Heat", "Cold"])
    out = _main(root, tmp_path)
    log = capsys.readouterr().out
    assert "reading:" in log and "skipped" not in log
    html = (out / "pathfinder_report.html").read_text(encoding="utf-8")
    assert ">Biological reading</button>" in html
    # the tab comes right after Summary
    tabs = re.findall(r"<button class='tab[^']*' onclick=\"show\(\d+\)\">([^<]+)</button>", html)
    assert tabs[:2] == ["Summary", "Biological reading"]
    for sid in ("landscape", "map", "programmes", "groups", "volcano1", "spotlight"):
        assert f"id='rd-{sid}'" in html, sid
    assert "no order" in html and "Largest step" not in html
    assert not re.search(r"""(src|href)=["']https?://""", html)       # offline report
    svgs = _svgs(html)
    assert len(svgs) >= 7
    for s in svgs:                                                   # well-formed SVG
        ET.fromstring(s.replace("&nbsp;", " "))
    assert (out / TSV_NAME).exists()


def test_categorical_layout_is_bars(tmp_path):
    root = make_deimos_run(str(tmp_path / "run"), ["Ctrl", "DrugA", "DrugB", "Heat"])
    run, d, prog = _prog(root)
    assert prog is not None and 3 <= prog.k <= 8
    _, facts = fig_programmes(run, d, prog)
    assert facts["layout"] == "bars"


def test_many_conditions_use_heatmap(tmp_path):
    tissues = [f"Tissue{chr(65 + i)}" for i in range(18)]
    root = make_deimos_run(str(tmp_path / "run"), tissues, n_rep=2, n_prot=500)
    run, d, prog = _prog(root)
    assert not d.ordered
    _, facts = fig_programmes(run, d, prog)
    assert facts["layout"] == "heatmap"
    assert "programmes" in _ids(_reading(root))


def test_time_course_uses_lines(tmp_path):
    conds = [f"{g}_D{t}" for g in ("WT", "KO") for t in (0, 3, 7, 14)]
    root = make_deimos_run(str(tmp_path / "run"), conds,
                           contrasts=[(f"{g}_D{t}", f"{g}_D0") for g in ("WT", "KO")
                                      for t in (3, 7, 14)] + [("KO_D7", "WT_D7")])
    run, d, prog = _prog(root)
    assert d.ordered and d.groups == ["WT", "KO"] and d.level_name == "time (d)"
    _, facts = fig_programmes(run, d, prog)
    assert facts["layout"] == "lines"
    res = _reading(root)
    mp = next(s for s in res["sections"] if s["id"] == "map")
    assert any("Largest step" in b for b in mp["bullets"])
    assert "marker-end" in mp["svgs"][0]                              # arrowed trajectories


def test_crossed_design_with_many_groups(tmp_path, capsys):
    """Hordeum-like: 4 genotypes x light x time = 8 ordered groups. Too many groups for
    trajectories: centroids coloured by level, heatmap of the programmes."""
    gens = ["Tron", "M1", "M2", "M3"]
    conds = [f"{g}_T0_Light" for g in gens] + [f"{g}_T{t}_{l}" for g in gens for t in (3, 7)
                                               for l in ("Light", "Dark")]
    contr = [(f"{g}_T0_Light", "Tron_T0_Light") for g in gens[1:]] + \
        [(f"{g}_T{t}_Light", f"{g}_T{t}_Dark") for g in gens for t in (3, 7)]
    root = make_deimos_run(str(tmp_path / "run"), conds, n_rep=2, n_prot=500, contrasts=contr)
    run, d, prog = _prog(root)
    assert d.ordered and len(d.groups) == 8 and d.level_name == "time"
    _, facts = fig_programmes(run, d, prog)
    assert facts["layout"] == "heatmap"
    capsys.readouterr()
    res = _reading(root)
    assert "skipped" not in capsys.readouterr().out
    mp = next(s for s in res["sections"] if s["id"] == "map")
    assert not any("Largest step" in b for b in mp["bullets"])


def test_two_conditions(tmp_path, capsys):
    root = make_deimos_run(str(tmp_path / "run"), ["Grp_Male4", "Grp_Female4"], n_rep=4)
    res = _reading(root)
    assert "skipped" not in capsys.readouterr().out
    ids = _ids(res)
    assert "programmes" not in ids                     # < 3 conditions: no programmes
    assert "landscape" in ids and "spotlight" in ids and res["prog"] is None


def test_minimal_workbook(tmp_path):
    root = make_deimos_run(str(tmp_path / "run"), ["A", "B", "C"], go=False, matrices=False,
                           pca=False, manifest=False)
    out = _main(root, tmp_path)
    html = (out / "pathfinder_report.html").read_text(encoding="utf-8")
    assert "id='rd-landscape'" in html and "id='rd-volcano1'" in html
    assert "id='rd-groups'" not in html and "id='rd-map'" not in html


def test_no_reading_flag(tmp_path):
    root = make_deimos_run(str(tmp_path / "run"), ["Ctrl", "DrugA", "DrugB", "Heat"])
    out = _main(root, tmp_path, "--no-reading")
    html = (out / "pathfinder_report.html").read_text(encoding="utf-8")
    assert "Biological reading" not in html and "<svg class='rd'" not in html
    assert "Programmes" not in pd.ExcelFile(out / "pathfinder_synthesis.xlsx").sheet_names
    assert not (out / TSV_NAME).exists()


def test_selection_is_the_proteins_tab_rule(tmp_path):
    """Same proteins as the Proteins tab and as Deimos' volcano_status, including when
    --pthr / --min-ratio override the upstream rule."""
    root = make_deimos_run(str(tmp_path / "run"), ["A", "B", "C", "D"])
    for pthr, ratio in ((None, None), (0.01, 1.0)):
        x, stats, man, rule, go_res, go_long = _load(root, pthr, ratio)
        run = load_run(stats, x, rule, man, verbose=False)
        prots = axis_proteins.run(stats, p_thr=rule.p_thr, rule=rule)
        union = np.zeros(len(run.de), bool)
        for c in run.contrasts:
            col = c + ("_p.adj" if rule.p_column == "p.adj" else "_p.val")
            p, d = run.de[col].to_numpy(), run.de[c + "_diff"].to_numpy()
            exp = (p < rule.p_thr) & (np.abs(d) >= rule.lfc_min)
            assert int(run.sig(c).sum()) == int(exp.sum())
            union |= run.sig(c).to_numpy()
        assert int(union.sum()) == prots["n_significant"]
    assert run.ratio_min == 1.0 and "overridden" in run.rule_source


def test_groups_are_the_metago_redundancy_groups(tmp_path):
    root = make_deimos_run(str(tmp_path / "run"), ["A", "B", "C", "D"], terms_per_programme=7)
    x, stats, man, rule, go_res, go_long = _load(root)
    run = load_run(stats, x, rule, man, verbose=False)
    g = build_groups(go_res, go_long, id_lookup(run.de))
    meta = go_res["meta"]
    for t in g.names:
        sub = g.terms[g.terms.theme == t]
        assert sub.term.map(meta["redundancy_cluster"]).nunique() == 1   # one group = one cluster
        reps = sub[sub.term.map(meta["is_representative"])]
        if len(reps):                                   # name = Meta-GO representative
            assert g.names[t] in set(reps.name)
    assert len(g.names) == g.terms.cluster.nunique() < len(g.terms)
    best = g.table.groupby("theme").p.min()
    assert list(best.sort_index()) == sorted(best)      # G1 = best term
    # every intersection points to a protein of the workbook
    assert g.n_unmapped == 0 and set().union(*g.terms.genes) <= set(run.de.index)


def test_intersections_reach_secondary_accessions(tmp_path):
    de = pd.DataFrame({"protein": ["P1", "P2"], "Protein.Group": ["P1;Q9", "P2"]}).set_index(
        "protein", drop=False)
    lk = id_lookup(de)
    assert lk["Q9"] == "P1" and lk["P2"] == "P2" and "X" not in lk


def test_volcano_tooltips_and_toggles(tmp_path):
    root = make_deimos_run(str(tmp_path / "run"), ["Ctrl", "DrugA", "DrugB", "Heat"])
    res = _reading(root)
    v = next(s for s in res["sections"] if s["id"] == "volcano1")["svgs"][0]
    tip = re.search(r"<title>([^<]*log2FC[^<]*)</title>", v).group(1)
    assert "p.adj" in tip and "programme P" in v
    assert "class='lg' data-t=" in v and "data-g='ns'" in v
    sp = next(s for s in res["sections"] if s["id"] == "spotlight")["svgs"][0]
    assert "log2 " in sp and "<title>" in sp


def test_tables_in_the_synthesis_workbook(tmp_path, capsys):
    root = make_deimos_run(str(tmp_path / "run"), ["Ctrl", "DrugA", "DrugB", "Heat", "Cold"])
    md = os.path.join(root, "matrices")                # DPC: 20 proteins never seen in Ctrl
    A = pd.read_csv(os.path.join(md, "log2_imputed.tsv"), sep="\t")
    s = pd.read_csv(os.path.join(md, "samples.tsv"), sep="\t")
    A.to_csv(os.path.join(md, "dpc_abundance.tsv"), sep="\t", index=False)
    O = A.copy()
    labs = list(s.label)
    O[labs] = (np.random.default_rng(1).random((len(O), len(labs))) > 0.05).astype(int)
    O.loc[:19, list(s.label[s.condition == "Ctrl"])] = 0
    O.to_csv(os.path.join(md, "dpc_observed.tsv"), sep="\t", index=False)
    gone = list(O.protein[:20])
    out = _main(root, tmp_path)
    x = pd.ExcelFile(out / "pathfinder_synthesis.xlsx")
    for sh in ("Reading_facts", "Programmes", "Programme_summary", "Groups", "Group_terms",
               "Group_x_contrast", "Reading_columns"):
        assert sh in x.sheet_names, sh
    pr = pd.read_excel(x, "Programmes")
    summ = pd.read_excel(x, "Programme_summary")
    assert pr.protein.is_unique
    assert pr.programme.value_counts().to_dict() == dict(zip(summ.programme, summ.n_proteins))
    assert pr.r_to_programme.between(-1, 1).all() and pr.silhouette.between(-1, 1).all()
    assert (pr.n_significant >= 1).all()
    for c in ["Ctrl", "DrugA", "DrugB", "Heat", "Cold"]:
        assert f"log2_{c}" in pr.columns and f"z_{c}" in pr.columns
    z = pr[[c for c in pr.columns if c.startswith("z_")]]
    assert np.allclose(z.mean(axis=1), 0, atol=1e-2)
    sel = pr[pr.protein.isin(gone)]
    assert len(sel) and (sel.min_observed_fraction == 0).all()
    assert (sel.least_observed_in == "Ctrl").all()
    g, gt, gx = (pd.read_excel(x, s) for s in ("Groups", "Group_terms", "Group_x_contrast"))
    assert set(gt.group_id) == set(g.group_id) == set(gx.group_id)
    assert g.best_fdr.is_monotonic_increasing and (gt.best_fdr < 0.05).all()
    assert (g.n_terms == gt.groupby("group_id").size().reindex(g.group_id).values).all()
    mg = pd.read_excel(x, "MetaGO_matrix", index_col=0)               # link to the Meta-GO
    assert set(g.metago_cluster) <= set(mg.redundancy_cluster)
    facts = pd.read_excel(x, "Reading_facts")
    assert "Response programmes" in set(facts.section)
    assert not facts.fact.str.contains("<").any()
    prov = pd.read_excel(x, "provenance").set_index("key")["value"]
    assert "reading_programmes" in prov.index and "DPC" in prov["reading_abundance"]
    tsv = pd.read_csv(out / TSV_NAME, sep="\t")
    assert list(tsv.protein) == list(pr.protein)
    summary = (out / "pathfinder_summary.md").read_text(encoding="utf-8")
    assert "## Response programmes (Biological reading)" in summary
    html = (out / "pathfinder_report.html").read_text(encoding="utf-8")
    assert "hollow = not detected" in html                          # DPC: non-detections


def test_cli_options_reach_the_reading(tmp_path):
    root = make_deimos_run(str(tmp_path / "run"), [f"Line{i}" for i in range(1, 6)])
    res = _reading(root, k=3, volcanos=["Line3_vs_Line1"],
                   design_overrides={"order_regex": r"Line(\d+)"})
    assert res["prog"].k == 3 and res["design"].ordered
    assert [s["title"] for s in res["sections"] if s["id"].startswith("volcano")] == \
        ["Volcano · Line3 vs Line1"]
    out = _main(root, tmp_path, "--k", "4", "--categorical", "--max-volcanos", "2")
    prov = pd.read_excel(out / "pathfinder_synthesis.xlsx", "provenance").set_index("key")["value"]
    assert prov["reading_programmes"].startswith("k = 4 (set by --k)")
    assert "forced categorical" in prov["reading_design"]
    html = (out / "pathfinder_report.html").read_text(encoding="utf-8")
    assert html.count("id='rd-volcano") == 2


def test_a_failing_figure_does_not_break_the_report(tmp_path, monkeypatch, capsys):
    import reading.run as rr
    root = make_deimos_run(str(tmp_path / "run"), ["Ctrl", "DrugA", "DrugB", "Heat"])

    def boom(*a, **k):
        raise RuntimeError("synthetic")
    monkeypatch.setattr(rr, "fig_map", boom)
    out = _main(root, tmp_path)
    assert "condition map skipped (RuntimeError: synthetic)" in capsys.readouterr().out
    html = (out / "pathfinder_report.html").read_text(encoding="utf-8")
    assert "id='rd-map'" not in html and "id='rd-programmes'" in html


# --------------------------------------------------------------------------- 1.3.0 additions

def test_gsea_groups_when_the_ora_is_empty(tmp_path, capsys):
    """H. glaber case: no usable ORA term, the GSEA of Deimos gives the functional groups
    (leading-edge proteins, NES as direction), grouped like the Meta-GO terms."""
    root = make_deimos_run(str(tmp_path / "run"), ["Ctrl", "DrugA", "DrugB", "Heat"],
                           go=False, gsea=True)
    res = _reading(root)
    g = res["groups"]
    assert g is not None and g.origin == "GSEA"
    assert "GSEA leading edges" in res["meta"]["reading_groups"]
    # 3 redundant pathways per programme -> grouped together
    assert len(g.names) < g.terms.term.nunique()
    assert set().union(*g.terms.genes) <= set(res["run"].de.index)
    sec = next(s for s in res["sections"] if s["id"] == "groups")
    assert "GSEA" in sec["bullets"][0] and "colour = NES" in sec["svgs"][0]
    assert "functional groups from the GSEA leading edges" in capsys.readouterr().out
    # the programmes are annotated with the GSEA groups
    assert any(res["prog"].themes.get(p) for p in range(1, res["prog"].k + 1))
    # --groups-from ora: no group without ORA; gsea forced even when an ORA exists
    assert _reading(root, groups_from="ora")["groups"] is None
    root2 = make_deimos_run(str(tmp_path / "run2"), ["Ctrl", "DrugA", "DrugB", "Heat"], gsea=True)
    assert _reading(root2)["groups"].origin == "ORA"
    assert _reading(root2, groups_from="gsea")["groups"].origin == "GSEA"


def test_gsea_sheets_without_lead_proteins_use_the_genes(tmp_path):
    from reading.gsea import load_gsea, term_parts
    assert term_parts("GO:BP__GO:0009414 response to water deprivation") == \
        ("GO:BP", "GO:0009414", "response to water deprivation")
    assert term_parts("GO_Biological_Process_2023__translation (GO:0006412)")[1:] == \
        ("GO:0006412", "translation")
    root = make_deimos_run(str(tmp_path / "run"), ["A", "B", "C"], go=False, gsea=True)
    x = _xlsx(root)
    sheets = pd.read_excel(x, sheet_name=None, header=None)
    with pd.ExcelWriter(x, engine="openpyxl") as w:               # drop Lead_proteins
        for n, d in sheets.items():
            if n.startswith("GSEA_"):
                d = d.loc[:, [j for j in d.columns if str(d.iat[0, j]) != "Lead_proteins"]]
            d.to_excel(w, sheet_name=n, index=False, header=False)
    _, stats, man, rule, _, _ = _load(root)
    run = load_run(stats, x, rule, man, verbose=False)
    gl = load_gsea(x, run.contrasts, man, run.de, id_lookup(run.de))
    assert gl is not None and gl.z.abs().min() > 1
    named = run.de["Genes"].astype(str).str.len() > 0
    for lst in gl.intersection:                                   # only proteins with a symbol
        assert lst and set(lst) <= set(run.de.index[named])


def test_programme_stability(tmp_path):
    """Resampling the replicates: well-separated programmes are stable (Jaccard >= 0.75);
    asking for twice the true number of programmes splits some of them arbitrarily."""
    from reading.programmes import programme_stability
    root = make_deimos_run(str(tmp_path / "run"), ["Ctrl", "DrugA", "DrugB", "Heat", "Cold"],
                           n_rep=4, k_true=4)
    run, d, prog = _prog(root)
    prog = programme_stability(run, prog, n_boot=20)
    assert prog.n_boot == 20 and set(prog.stability) == set(range(1, prog.k + 1))
    assert min(prog.stability.values()) >= 0.75
    assert prog.protein_stability.between(0, 1).all()
    x, stats, man, rule, go_res, go_long = _load(root)
    over = build_programmes(run, d, None, k=8)
    over = programme_stability(run, over, n_boot=20)
    assert min(over.stability.values()) < min(prog.stability.values())
    res = _reading(root, n_boot=10)
    sec = next(s for s in res["sections"] if s["id"] == "programmes")
    assert "Stability under replicate resampling" in sec["bullets"][0] and "stability" in sec["svgs"][0]
    assert "stability_jaccard" in res["tables"]["Programme_summary"]
    assert "stability" in res["tables"]["Programmes"]
    assert _reading(root, n_boot=0)["prog"].stability is None


HV_GENS = ["Tron", "M1", "M2", "M3"]


def _hordeum(tmp_path, samples_columns=False):
    gens = HV_GENS
    conds = [f"{g}_T0_Light" for g in gens] + [f"{g}_T{t}_{l}" for g in gens for t in (3, 7)
                                               for l in ("Light", "Dark")]

    def f(c):
        g, t, l = c.split("_")
        return gens.index(g), int(t[1:]), l
    eff = {1: [2.0 * (f(c)[0] in (1, 2)) for c in conds],                       # genotype
           2: [(1.5 if f(c)[2] == "Dark" else -0.5) * (f(c)[1] > 0) for c in conds],  # light
           3: [1.8 * (f(c)[0] == 3) + (1.2 if f(c)[2] == "Dark" else 0) for c in conds]}
    contr = [(f"{g}_T0_Light", "Tron_T0_Light") for g in gens[1:]] + \
        [(f"{g}_T{t}_Dark", f"{g}_T{t}_Light") for g in gens for t in (3, 7)]
    root = make_deimos_run(str(tmp_path / "run"), conds, n_rep=3, n_prot=700, k_true=3,
                           contrasts=contr, cond_effects=eff)
    if samples_columns:
        p = os.path.join(root, "matrices", "samples.tsv")
        s = pd.read_csv(p, sep="\t")
        s["genotype"] = s.condition.str.split("_").str[0]
        s["time"] = s.condition.str.split("_").str[1]
        s["light"] = s.condition.str.split("_").str[2]
        s.to_csv(p, sep="\t", index=False)
    return root


def test_factorial_design_variance_partition(tmp_path):
    from reading.factorial import design_factors, variance_partition
    root = _hordeum(tmp_path)
    run, d, prog = _prog(root)
    F, src = design_factors(run, d.conditions)
    assert src == "condition names" and set(F.columns) == {"Tron/M1/M2/M3", "time", "Light/Dark"}
    S, terms = variance_partition(run, list(prog.z.index), F)
    assert "Tron/M1/M2/M3 × Light/Dark" in terms and {"higher-order", "replicates"} <= set(S)
    assert S.drop(columns="sum").ge(0).all().all()
    lab = prog.labels.reindex(S.index)
    med = S.groupby(lab).median()
    geno = med["Tron/M1/M2/M3"].idxmax()                       # a genotype-driven programme
    light = med["Light/Dark"].idxmax()                         # a light-driven programme
    assert med.loc[geno, "Tron/M1/M2/M3"] > 0.7 and med.loc[light, "Light/Dark"] > 0.7
    assert geno != light
    res = _reading(root, n_boot=0)
    sec = next(s for s in res["sections"] if s["id"] == "factors")
    assert "Design factors (from condition names)" in sec["bullets"][0]
    assert "variance share per design factor" in sec["svgs"][0]
    t = res["tables"]
    assert {"Variance_partition", "Design_factors"} <= set(t)
    assert "var_Light/Dark" in t["Programmes"] and "reading_factors" in res["meta"]


def test_factorial_factors_from_samples_columns(tmp_path):
    from reading.factorial import design_factors
    root = _hordeum(tmp_path, samples_columns=True)
    run, d, prog = _prog(root)
    F, src = design_factors(run, d.conditions)
    assert src == "samples.tsv columns" and list(F.columns) == ["genotype", "time", "light"]
    # a single-factor design is not factorial
    root1 = make_deimos_run(str(tmp_path / "one"), ["Ctrl", "DrugA", "DrugB", "Heat"])
    run1, d1, _ = _prog(root1)
    assert design_factors(run1, d1.conditions)[0] is None
    assert "factors" not in _ids(_reading(root1, n_boot=0))


HV45 = ["Tron"] + [f"M{i}" for i in range(1, 9)]


def _hordeum45(tmp_path, n_prot=700):
    """The real Hordeum layout: 9 genotypes x (T0, T3 / T7 x Light / Dark) = 45 groups, T0
    named without the light (Tron_T0); mutants vs Tron at T0, Dark vs Light at T3 and T7."""
    conds = [f"{g}_T0" for g in HV45] + [f"{g}_T{t}_{l}" for g in HV45 for t in (3, 7)
                                         for l in ("Light", "Dark")]

    def f(c):
        p = c.split("_")
        return HV45.index(p[0]), int(p[1][1:]), (p[2] if len(p) > 2 else "Light")
    eff = {1: [2.0 * (f(c)[0] in (2, 5)) for c in conds],                     # genotype
           2: [(1.6 if f(c)[2] == "Dark" else -0.4) * (f(c)[1] / 7) for c in conds],  # light x time
           3: [1.8 * (f(c)[0] == 7) + (1.2 if f(c)[2] == "Dark" and f(c)[1] == 7 else 0)
               for c in conds]}
    contr = [(f"{g}_T0", "Tron_T0") for g in HV45[1:]] + \
        [(f"{g}_T{t}_Dark", f"{g}_T{t}_Light") for g in HV45 for t in (3, 7)]
    return make_deimos_run(str(tmp_path / "run"), conds, n_rep=4, n_prot=n_prot, k_true=3,
                           contrasts=contr, cond_effects=eff), conds


def test_irregular_names_45_groups(tmp_path):
    from reading.factorial import (MISSING, _from_names, condition_order, contrast_factors,
                                   factor_blocks, levels, roles)
    _, conds = _hordeum45(tmp_path, n_prot=300)
    F = pd.DataFrame(_from_names(conds))
    assert set(F.columns) == {"time", "factor 1", "Light/Dark"}
    assert F.loc["Tron_T0", "Light/Dark"] == MISSING and F.loc["M3_T7_Dark", "factor 1"] == "M3"
    r = roles(F)
    assert r == {"x": "time", "lines": "factor 1", "facets": ["Light/Dark"]}
    assert levels(F, "time") == ["T0", "T3", "T7"]
    order = condition_order(F, r)
    assert order[:5] == ["Tron_T0", "Tron_T3_Light", "Tron_T7_Light", "Tron_T3_Dark", "Tron_T7_Dark"]
    pairs = {f"{g}_T0_vs_Tron_T0": (f"{g}_T0", "Tron_T0") for g in HV45[1:]}
    pairs.update({f"{g}_T{t}_Dark_vs_{g}_T{t}_Light": (f"{g}_T{t}_Dark", f"{g}_T{t}_Light")
                  for g in HV45 for t in (3, 7)})
    cf = contrast_factors(pairs, F, list(pairs))
    assert cf["M3_T0_vs_Tron_T0"][0] == ["factor 1"]
    assert cf["M3_T7_Dark_vs_M3_T7_Light"][:2] == (["Light/Dark"], {"time": "T7", "factor 1": "M3"})
    assert [(k, len(v)) for k, v in factor_blocks(cf, F)] == [("Light/Dark", 18), ("factor 1", 8)]
    # a short name with a token no position knows is not guessed
    assert _from_names(["A_T0", "A_T3_x", "B_T3_y", "Q_T7"]) == {}


def test_factor_regex_names_the_factors(tmp_path, capsys):
    from reading.factorial import design_factors
    root, conds = _hordeum45(tmp_path, n_prot=300)
    run, d, prog = _prog(root)
    rx = r"(?P<genotype>[^_]+)_T(?P<time>\d+)(?:_(?P<light>Light|Dark))?"
    F, src = design_factors(run, d.conditions, regex=rx)
    assert src == "--factor-regex" and list(F.columns) == ["genotype", "time", "light"]
    assert F.loc["M1_T0", "light"] == "—" and F.loc["M1_T3_Dark", "time"] == "3"
    with pytest.raises(ValueError, match="does not match"):
        design_factors(run, d.conditions, regex=r"(?P<g>X\d)")
    res = _reading(root, n_boot=0, factor_regex=r"(?P<g>X\d)")           # bad regex: fallback
    assert "factors read automatically instead" in capsys.readouterr().out
    assert "factor-effects" in _ids(res)


def test_factorial_views_end_to_end_45_groups(tmp_path):
    root, conds = _hordeum45(tmp_path)
    rx = r"(?P<genotype>[^_]+)_T(?P<time>\d+)(?:_(?P<light>Light|Dark))?"
    out = _main(root, tmp_path, "--factor-regex", rx, "--prog-boot", "5")
    html = (out / "pathfinder_report.html").read_text(encoding="utf-8")
    for sid in ("factor-effects", "map", "programmes", "factor-profiles", "factors"):
        assert f"id='rd-{sid}'" in html, sid
    for s in _svgs(html):
        ET.fromstring(s)
    eff = re.search(r"id='rd-factor-effects'.*?</svg>", html, re.S).group(0)
    assert "light varied — 18 comparison(s)" in eff and "genotype varied — 8 comparison(s)" in eff
    assert "up = higher in Dark than in Light" in eff
    prof = re.search(r"id='rd-factor-profiles'.*?</svg>", html, re.S).group(0)
    assert "one line per genotype" in prof and "Tron (the reference" in prof
    assert "shared start" in prof and "data-g='l0'" in prof
    mp = re.search(r"id='rd-map'.*?</svg>", html, re.S).group(0)
    assert "colour = genotype · shape = light · size = time" in mp
    # programme heatmap columns follow genotype > light > time
    hm = re.search(r"id='rd-programmes'.*?</svg>", html, re.S).group(0)
    hm = hm[hm.index("<svg"):]
    assert hm.index("Tron_T0") < hm.index("Tron_T3_Light") < hm.index("Tron_T3_Dark") < hm.index("M1_T0")
    x = pd.ExcelFile(out / "pathfinder_synthesis.xlsx")
    vp = pd.read_excel(x, "Variance_partition")
    assert {"genotype", "light", "time", "genotype × light"} <= set(vp.columns)
    prov = pd.read_excel(x, "provenance").set_index("key")["value"]
    assert prov["reading_factor_roles"] == "x = time, lines = genotype, panels = light"


def test_glued_sex_age_names_are_a_2x2_design(tmp_path):
    """H. glaber run: Grp_Male4 / Grp_Female4 / Grp_Male14 / Grp_Female14 = sex x age; the
    conditions are written Grp-Male4 in samples.tsv while the contrasts say Grp_Male4."""
    from reading.factorial import _from_names, contrast_factors, factor_blocks, roles
    conds = ["Grp-Male4", "Grp-Female4", "Grp-Male14", "Grp-Female14"]
    F = pd.DataFrame(_from_names(conds))
    assert list(F.columns) == ["Male/Female", "4/14"]               # 4 / 14 is not 'time'
    assert roles(F) == {"x": "4/14", "lines": "Male/Female", "facets": []}
    pairs = {"Grp_Male4_vs_Grp_Female4": ("Grp_Male4", "Grp_Female4"),
             "Grp_Male4_vs_Grp_Male14": ("Grp_Male4", "Grp_Male14"),
             "Grp_Female4_vs_Grp_Male14": ("Grp_Female4", "Grp_Male14")}
    cf = contrast_factors(pairs, F, list(pairs))
    assert cf["Grp_Male4_vs_Grp_Female4"][0] == ["Male/Female"] and len(cf) == 3
    assert [k for k, _ in factor_blocks(cf, F)] == ["Male/Female", "4/14", "several"]
    # no glued split when the separators already give the factors, nor for plain series
    assert _from_names([f"Line{i}" for i in range(1, 6)]) == {}
    assert "time" in _from_names(["WT_0h", "WT_6h", "WT_24h", "KO_0h", "KO_6h", "KO_24h"])
    emb = pd.DataFrame(_from_names([f"emb{c}{t}" for c in ("Aviso", "Major") for t in (0, 2, 8, 24)]))
    assert roles(emb)["x"] == "0/2/8/24" and roles(emb)["lines"] == "embAviso/embMajor"
