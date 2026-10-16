"""Rebuild tests/fixtures/deimos_contract.xlsx from a real Deimos workbook.

Keeps only what Pathfinder reads (values only, no images, hidden _manifest kept
hidden), so the fixture stays a few hundred KB. Regenerate it whenever Deimos bumps
its manifest schema version:

    python tests/fixtures/make_deimos_fixture.py path/to/ProteomicAnalysis_Results.xlsx

The committed fixture comes from Deimos on its example_data (16 samples, 4
conditions renamed in the strain style 'Live_K12' / 'Heat_killed_M18', WGCNA on),
with the gProfiler HTTP call replaced by a deterministic offline stand-in: the GO
sheets have Deimos' real layout and names, the terms themselves are synthetic.
"""
import sys
from pathlib import Path

import openpyxl

KEEP = ("Methods", "Log2_Impute", "Differential_Expression", "Zscore_Heatmap",
        "WGCNA_Results", "WGCNA_Modules", "WGCNA_kME", "WGCNA_Eigengenes", "_manifest")


def build(src, dst):
    wb_in = openpyxl.load_workbook(src, read_only=True)
    wb_out = openpyxl.Workbook()
    wb_out.remove(wb_out.active)
    for name in wb_in.sheetnames:
        if name not in KEEP and not name.startswith("GO_"):
            continue
        ws_in = wb_in[name]
        ws_out = wb_out.create_sheet(name)
        for row in ws_in.iter_rows(values_only=True):
            ws_out.append(list(row))
        if name == "_manifest":
            ws_out.sheet_state = "hidden"
    wb_out.save(dst)
    return dst


if __name__ == "__main__":
    out = Path(__file__).with_name("deimos_contract.xlsx")
    build(sys.argv[1], out)
    print(f"{out} ({out.stat().st_size // 1024} KB)")
