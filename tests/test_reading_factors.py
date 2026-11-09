"""Design reading (reading.factors, ex-InSight): an order is read only with evidence."""
import pandas as pd
import pytest

from reading.factors import detect_design

ORDERED = {
    "crossed (P5)": ([f"emb{g}{t}" for g in ("Aviso", "Major") for t in (0, 2, 8, 24)],
                     ["embAviso", "embMajor"], "level", ["0", "2", "8", "24"]),
    "days prefix": (["D0", "D3", "D7", "D14"], [], "time (d)", ["0 d", "3 d", "7 d", "14 d"]),
    "groups x days": (["Ctrl_D0", "Ctrl_D7", "Ctrl_D14", "KO_D0", "KO_D7", "KO_D14"],
                      ["Ctrl", "KO"], "time (d)", ["0 d", "7 d", "14 d"]),
    "hours suffix": (["WT_0h", "WT_6h", "WT_24h", "mut_0h", "mut_6h", "mut_24h"],
                     ["WT", "mut"], "time (h)", ["0 h", "6 h", "24 h"]),
    "T inside": (["Tron_T0_Light", "Tron_T3_Light", "Tron_T7_Light", "Tron_T3_Dark",
                  "Tron_T7_Dark"], ["Tron · Light", "Tron · Dark"], "time", ["T0", "T3", "T7"]),
    "unit inside": (["WT_24h_Light", "WT_48h_Light", "WT_72h_Light", "WT_24h_Dark",
                     "WT_48h_Dark", "WT_72h_Dark"], ["WT · Light", "WT · Dark"], "time (h)",
                    ["24 h", "48 h", "72 h"]),
    "DAP first": (["15DAP", "20DAP", "30DAP"], [], "time (d)", ["15 d", "20 d", "30 d"]),
    "hai": (["Seed_0hai", "Seed_12hai", "Seed_24hai"], [], "time (h)", ["0 h", "12 h", "24 h"]),
    "dose": (["Dose0", "Dose10", "Dose100"], [], "dose", ["0", "10", "100"]),
}
CATEGORICAL = {
    "two groups": ["Grp_Male4", "Grp_Female4"],
    "lines (IDs, not an order)": ["Line1", "Line2", "Line3", "Line4"],
    "treatments": ["Ctrl", "DrugA", "DrugB", "Heat"],
    "too few levels": ["WT1", "WT2", "KO1", "KO2"],
    "sites": ["S1", "S2", "S3", "S4"],
    "clones": ["WT_clone1", "WT_clone2", "WT_clone3"],
    "tissues": ["Liver", "Brain", "Heart", "Kidney", "Lung"],
}


@pytest.mark.parametrize("name", ORDERED)
def test_ordered(name):
    conds, groups, lname, labels = ORDERED[name]
    d = detect_design(conds)
    assert d.ordered, name
    assert d.groups == groups
    assert d.level_name == lname
    assert [d.label(lv) for lv in d.levels] == labels
    assert sorted(d.conditions) == sorted(conds)


@pytest.mark.parametrize("name", CATEGORICAL)
def test_categorical(name):
    d = detect_design(CATEGORICAL[name])
    assert not d.ordered, name
    assert d.conditions == CATEGORICAL[name]          # original order kept


def test_overrides():
    conds = ["Line1", "Line2", "Line3"]
    d = detect_design(conds, overrides={"order_regex": r"Line(\d+)"})
    assert d.ordered and d.source == "order_regex"
    d = detect_design([f"emb{g}{t}" for g in "AB" for t in (0, 2, 8)],
                      overrides={"categorical": True})
    assert not d.ordered
    s = pd.DataFrame({"label": ["a", "b", "c", "d", "e", "f"],
                      "condition": ["x", "y", "z", "u", "v", "w"],
                      "stage": [1, 2, 3, 1, 2, 3], "line": ["L1"] * 3 + ["L2"] * 3})
    d = detect_design(list("xyzuvw"), s, {"order_column": "stage", "group_column": "line"})
    assert d.ordered and d.groups == ["L1", "L2"] and d.levels == [1, 2, 3]


def test_display_names_strip_common_parts():
    d = detect_design(["Tissue_Liver", "Tissue_Brain", "Tissue_Heart"])
    assert [d.display(c) for c in d.conditions] == ["Liver", "Brain", "Heart"]
    d = detect_design(["A_rep", "B_rep", "C_rep"])
    assert d.display("A_rep") == "A"
    d = detect_design(["Ctrl", "DrugA"])                 # nothing in common: unchanged
    assert d.display("DrugA") == "DrugA"
