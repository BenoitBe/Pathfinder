"""Species resolution: built-in table, STRING species list, ambiguity handling."""
import builtins

import pytest

import species

FAKE_SPECIES = [
    (3483, "Cannabis sativa"),
    (90675, "Camelina sativa"),
    (60357, "Castanea sativa"),
    (3708, "Brassica napus"),
    (9606, "Homo sapiens"),
    (4530, "Oryza sativa"),
    (10090, "Mus musculus"),
]


class TestBuiltInTable:
    @pytest.mark.parametrize("query,taxid", [
        ("hsapiens", 9606), ("Homo sapiens", 9606), ("homo_sapiens", 9606),
        ("human", 9606), ("rnorvegicus", 10116), ("bnapus", 3708),
        ("10116", 10116), ("3708", 3708),
    ])
    def test_known_spellings(self, query, taxid):
        assert species.normalise_org(query)[0] == taxid

    @pytest.mark.parametrize("value,taxid", [
        ("NCBI taxid 1304", 1304), ("ncbi taxid 1304", 1304), ("taxid:1304", 1304),
        ("taxid 1304", 1304), ("NCBI_taxid_511145", 511145), ("taxon id = 9606", 9606),
    ])
    def test_labelled_taxid_as_written_by_deimos(self, value, taxid):
        """Deimos records 'Species = NCBI taxid 1304' in Methods in STRING mode."""
        assert species.normalise_org(value) == (taxid, None)

    def test_taxid_inside_free_text_is_not_taken(self):
        """A strain taxid quoted in a Methods detail is not the analysed organism."""
        assert species.normalise_org("S. salivarius K12 (taxid 1200793), via UniProt")[0] is None

    def test_unknown_returns_nothing_rather_than_guessing(self):
        assert species.normalise_org("csativa")[0] is None
        assert species.normalise_org("Camelina sativa")[0] is None
        assert species.normalise_org("not-a-species")[0] is None


class TestStringLookup:
    """REGRESSION: a hard-coded table of 20 organisms cannot serve a lab working on
    varied biological contexts; 'csativa' was simply rejected."""

    def test_full_scientific_name(self):
        assert species.lookup_species("Camelina sativa", FAKE_SPECIES) == [(90675, "Camelina sativa")]

    def test_genus_alone(self):
        assert species.lookup_species("cannabis", FAKE_SPECIES) == [(3483, "Cannabis sativa")]

    def test_gprofiler_style_code(self):
        assert species.lookup_species("bnapus", FAKE_SPECIES) == [(3708, "Brassica napus")]

    def test_ambiguous_code_returns_every_candidate(self):
        """'csativa' fits Cannabis, Camelina and Castanea sativa: never pick one."""
        got = species.lookup_species("csativa", FAKE_SPECIES)
        assert len(got) == 3
        assert {t for t, _ in got} == {3483, 90675, 60357}

    def test_unknown_returns_empty(self):
        assert species.lookup_species("xyzzy", FAKE_SPECIES) == []


class TestCandidateChoice:
    def test_single_candidate_is_taken(self, capsys):
        assert species.choose_candidate([(90675, "Camelina sativa")]) == 90675
        assert "Camelina sativa" in capsys.readouterr().out

    def test_no_candidate(self):
        assert species.choose_candidate([]) is None

    def test_ambiguity_without_a_terminal_is_refused(self, monkeypatch, capsys):
        """Non-interactive runs must not silently pick the first match."""
        monkeypatch.setattr(species.sys.stdin, "isatty", lambda: False, raising=False)
        got = species.choose_candidate(FAKE_SPECIES[:3])
        assert got is None
        assert "ambiguous" in capsys.readouterr().out

    def test_user_pick(self, monkeypatch):
        monkeypatch.setattr(species.sys.stdin, "isatty", lambda: True, raising=False)
        monkeypatch.setattr(builtins, "input", lambda *a: "2")
        assert species.choose_candidate(FAKE_SPECIES[:3]) == 90675

    def test_user_skips(self, monkeypatch):
        monkeypatch.setattr(species.sys.stdin, "isatty", lambda: True, raising=False)
        monkeypatch.setattr(builtins, "input", lambda *a: "")
        assert species.choose_candidate(FAKE_SPECIES[:3]) is None

    def test_out_of_range_pick_is_refused(self, monkeypatch):
        monkeypatch.setattr(species.sys.stdin, "isatty", lambda: True, raising=False)
        monkeypatch.setattr(builtins, "input", lambda *a: "99")
        assert species.choose_candidate(FAKE_SPECIES[:3]) is None


class TestWorkbookDetection:
    def test_deimos_string_methods_row(self, tmp_path):
        import pandas as pd
        p = tmp_path / "wb.xlsx"
        methods = pd.DataFrame({"Item": ["Tool", "Species"],
                                "Choice / Value": ["STRING (API 'enrichment')", "NCBI taxid 1304"],
                                "Method / Detail": ["", "Streptococcus salivarius - strain K12 "
                                                        "(taxid 1200793)"]})
        with pd.ExcelWriter(p) as xw:
            methods.to_excel(xw, sheet_name="Methods", index=False)
        assert species.detect_org_in_workbook(str(p)) == "1304"
