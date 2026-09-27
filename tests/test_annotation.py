"""Annotation: library search, lipid rules and harmonization, on the synthetic study."""

import numpy as np
import pandas as pd
import pytest

from atlas_ms.annotation.harmonize import harmonize
from atlas_ms.annotation.lipids import (
    annotate_feature, combine_chains, load_rules, parse_lipid, rules_path, species_table,
)
from atlas_ms.annotation.schema import candidate_table
from atlas_ms.config import HarmonizationSettings, LipidSettings
from test_workflow import feature_of


def annotations(project, name: str) -> pd.DataFrame:
    return pd.read_parquet(project.results_dir / "annotations" / f"{name}.parquet")


# ---- Lipid rules -------------------------------------------------------------

def test_species_formulas_agree_with_goslin():
    """Every class's species formulas, derived from its reference, match Goslin's."""
    from pygoslin.parser.Parser import LipidParser

    parser = LipidParser()
    species = species_table(load_rules(rules_path()))
    # Per class: the smallest, the largest and the most unsaturated species.
    for _, group in species[species["adduct"] == species.groupby("lipid_class")["adduct"].transform("first")].groupby("lipid_class"):
        for _, row in group.iloc[[0, -1, group["double_bonds"].argmax()]].iterrows():
            try:
                expected = parser.parse(row["name"]).get_sum_formula()
            except Exception:  # classes Goslin does not know (OL, DGTS)
                continue
            assert row["formula"] == expected, row["name"]
    # Ion m/z: PC 34:1 [M+H]+ and TG 52:2 [M+NH4]+ of the synthetic study.
    pc = species[(species["name"] == "PC 34:1") & (species["adduct"] == "[M+H]+")]
    tg = species[(species["name"] == "TG 52:2") & (species["adduct"] == "[M+NH4]+")]
    assert pc["mz"].iloc[0] == pytest.approx(760.5851, abs=1e-4)
    assert tg["mz"].iloc[0] == pytest.approx(876.8014, abs=1e-4)


def test_chain_combinations():
    # TG 52:2 with losses of 16:0 and 18:1: only 16:0_18:1_18:1 uses both.
    assert combine_chains([(16, 0), (18, 1)], 3, (52, 2)) == [(16, 0), (18, 1), (18, 1)]
    # TG 54:5 with losses of 16:0 and 18:1: the third chain (unseen) is 20:4.
    assert combine_chains([(16, 0), (18, 1)], 3, (54, 5)) == [(16, 0), (18, 1), (20, 4)]
    # TG 52:2 with four losses: several compositions fit, so none is chosen.
    assert combine_chains([(16, 0), (18, 0), (18, 1), (18, 2)], 3, (52, 2)) is None


def test_isobars_are_separated_by_their_headgroup():
    """PC 34:1 and PE 37:1 [M+H]+ have the same formula: the fragments decide."""
    rules = load_rules(rules_path())
    species = species_table(rules)

    def best(fragments):
        mz = np.array(sorted(fragments))
        rows = annotate_feature(1, 760.5851, (mz, np.ones_like(mz)), species, rules, LipidSettings())
        return rows[0]["name"], rows[0]["proposed_level"]

    assert best([184.0733, 104.1070]) == ("PC 34:1", "3")
    assert best([760.5851 - 141.0191]) == ("PE 37:1", "3")
    name, level = best([300.0])  # no class evidence: m/z suggestions only
    assert level == "5"


def test_lipid_rules_on_the_study(processed_project):
    features = pd.read_parquet(processed_project.results_dir / "features.parquet")
    table = annotations(processed_project, "lipid_rules")
    best = table[table["rank"] == 1].set_index("feature_id")
    expected = {
        "PC 34:1 [M+H]+": ("PC 34:1", "[M+H]+", "species", "3"),
        "PC 34:1 [M+Na]+": ("PC 34:1", "[M+Na]+", "species", "3"),
        "SM 34:1;O2 [M+H]+": ("SM 18:1;O2/16:0", "[M+H]+", "sn-position", "3"),
        "LPC 16:0 [M+H]+": ("LPC 16:0", "[M+H]+", "molecular species", "3"),
        "Cer 34:1;O2 [M+H]+": ("Cer 18:1;O2/16:0", "[M+H]+", "sn-position", "3"),
        "TG 52:2 [M+NH4]+": ("TG 16:0_18:1_18:1", "[M+NH4]+", "molecular species", "3"),
        "CE 18:1 [M+NH4]+": ("CE 18:1", "[M+NH4]+", "molecular species", "3"),
        # No MS2 spectrum: m/z suggestion only.
        "DG 34:1 [M+NH4]+": ("DG 34:1", "[M+NH4]+", "species", "5"),
    }
    for compound, values in expected.items():
        row = best.loc[feature_of(features, compound)["feature_id"]]
        assert (row["name"], row["adduct"], row["lipid_level"], row["proposed_level"]) == values, compound


def test_lipid_names_from_any_source():
    assert parse_lipid("PC(16:0/18:1(9Z))") | {} == {
        "name": "PC 16:0/18:1(9Z)", "species": "PC 34:1", "lipid_class": "PC", "level": "full structure"}
    assert parse_lipid("CE 18:1")["lipid_class"] == "CE"
    assert parse_lipid("PC O-34:1")["lipid_class"] == "PC-O"
    assert parse_lipid("OL 34:1")["lipid_class"] == "OL"  # unknown to Goslin, read directly
    assert parse_lipid("Caffeine") is None


# ---- Library search -----------------------------------------------------------

def test_library_search_levels(processed_project):
    features = pd.read_parquet(processed_project.results_dir / "features.parquet")
    table = annotations(processed_project, "library").set_index("feature_id")

    def hit(compound):
        return table.loc[feature_of(features, compound)["feature_id"]]

    # Reference standard, same spectrum and retention time: level 1.
    pc = hit("PC 34:1 [M+H]+")
    assert (pc["source"], pc["name"], pc["proposed_level"]) == ("library:standards", "PC 16:0_18:1", "1")
    assert pc["score"] == pytest.approx(1.0, abs=1e-3) and pc["matched_peaks"] == 3
    assert abs(pc["rt_error_s"]) <= 10
    # Standard without retention time, or with the wrong one: level 2a.
    assert hit("SM 34:1;O2 [M+H]+")["proposed_level"] == "2a"
    cer = hit("Cer 34:1;O2 [M+H]+")
    assert cer["proposed_level"] == "2a" and cer["rt_error_s"] == pytest.approx(200 - 300, abs=5)
    # In-silico spectrum: level 3.
    tg = hit("TG 52:2 [M+NH4]+")
    assert (tg["source"], tg["library_kind"], tg["proposed_level"]) == ("library:insilico", "in_silico", "3")
    # The negative-mode spectrum of the standards library was skipped.
    prepared = pd.read_parquet(processed_project.root / "work" / "annotation" / "libraries" / "standards.parquet")
    assert len(prepared) == 3
    # Library spectra are kept for mirror plots.
    assert 184.0733 == pytest.approx(max(pc["reference_mz"], key=lambda mz: -abs(mz - 184.07)), abs=1e-4)


# ---- Harmonization ------------------------------------------------------------

def test_best_annotations_of_the_study(processed_project):
    features = pd.read_parquet(processed_project.results_dir / "features.parquet")
    best = annotations(processed_project, "best").set_index("feature_id")
    assert len(best) == len(features)  # one row per feature, annotated or not

    def label(compound):
        return best.loc[feature_of(features, compound)["feature_id"], "label"]

    assert label("PC 34:1 [M+H]+") == "L1 · PC 16:0_18:1 (molecular species)"
    assert label("Cer 34:1;O2 [M+H]+") == "L2a · Cer 18:1;O2/16:0 (sn-position)"
    assert label("TG 52:2 [M+NH4]+") == "L3 · TG 16:0_18:1_18:1 (molecular species)"
    assert label("DG 34:1 [M+NH4]+") == "L5 · DG 34:1 (species)"
    # The phosphocholine family (PC, SM, LPC...): consensus class PC.
    pc = best.loc[feature_of(features, "PC 34:1 [M+H]+")["feature_id"]]
    assert pc["family_class"] == "PC" and 0 < pc["family_class_score"] <= 1


def candidate(feature_id, source, name, level, **values):
    return {"feature_id": feature_id, "source": source, "rank": 1, "name": name, "proposed_level": level,
            "score": 0.9, "library_kind": "experimental" if source.startswith("library") else ""} | values


def test_levels_are_capped_never_raised():
    candidates = candidate_table([
        candidate(1, "library:lipidblast", "PC 34:1", "2a", library_kind="in_silico"),  # in silico: 3
        candidate(2, "library:standards", "SM 34:1;O2", "1"),  # level 1 without retention time: 2a
        candidate(3, "sirius:csi", "PE 34:1", "2a"),  # CSI:FingerID: at most 3
        candidate(4, "lipid_rules", "PG 34:1", "5"),  # a level is never raised
    ])
    features = pd.DataFrame({"feature_id": [1, 2, 3, 4], "rt": [60.0] * 4})
    nodes = pd.DataFrame({"feature_id": [1, 2, 3, 4], "family": [-1] * 4, "qc": [""] * 4})
    result, _ = harmonize(candidates, features, nodes, HarmonizationSettings())
    assert result.set_index("feature_id")["level"].to_dict() == {1: "3", 2: "2a", 3: "3", 4: "5"}


def test_flags_class_disagreement_and_rt_trend():
    # Six PC species on a clean RT trend (RT rises 0.5 min per carbon, falls
    # 0.4 min per double bond), one of them 2 min late. Feature 7: two
    # sources disagree on its class. Feature 8 is a blank feature.
    compositions = [(32, 0), (34, 1), (34, 2), (36, 1), (36, 2), (38, 4)]
    rows, rt = [], {}
    for feature_id, (c, d) in enumerate(compositions, start=1):
        rows.append(candidate(feature_id, "lipid_rules", f"PC {c}:{d}", "3", lipid_class="PC"))
        rt[feature_id] = 60 * (0.5 * c - 0.4 * d + (2.0 if feature_id == 4 else 0.0))
    rows += [candidate(7, "lipid_rules", "PC 34:1", "3", lipid_class="PC"),
             candidate(7, "library:std", "PE 37:1", "2a"),
             candidate(8, "lipid_rules", "PG 34:1", "3", lipid_class="PG")]
    rt |= {7: 1000.0, 8: 1000.0}
    features = pd.DataFrame({"feature_id": list(rt), "rt": list(rt.values())})
    nodes = pd.DataFrame({"feature_id": list(rt), "family": [1] * 6 + [-1, -1],
                          "qc": [""] * 7 + ["present in blanks"]})
    _, best = harmonize(candidate_table(rows), features, nodes, HarmonizationSettings())
    flags = best.set_index("feature_id")["flags"]
    assert flags[4] == "RT off class trend"
    assert all(flags[i] == "" for i in (1, 2, 3, 5, 6))
    assert "lipid class disagreement" in flags[7]
    assert flags[8] == "blank"
    # Family 1: six PC annotations, all agreeing.
    family = best.set_index("feature_id").loc[1]
    assert (family["family_class"], family["family_class_score"]) == ("PC", 1.0)


# ---- Library harmonization ------------------------------------------------------

POPC = "CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCC=CCCCCCCCC"
CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"


def test_library_harmonization(tmp_path):
    """FragHub-like clean-up of two messy libraries, and duplicates between them."""
    from atlas_ms.annotation.libraries import combine_libraries, prepare_library, write_library_mgf
    from synthetic import LibraryEntry, write_msp

    pc_peaks = [(104.107, 5.0), (184.0733, 100.0), (577.519, 3.0)]
    messy = write_msp(tmp_path / "messy.msp", [
        # Adduct written without brackets, no ion mode, a lipid abbreviation as name, a correct structure.
        LibraryEntry("POPC", 760.5851, "M+H", "", pc_peaks, ionmode="", smiles=POPC),
        # A structure that cannot be this precursor (caffeine at m/z 760.59): removed.
        LibraryEntry("Wrong", 760.5851, "[M+H]+", "", pc_peaks, smiles=CAFFEINE),
        # A single fragment: removed (min_library_peaks = 2).
        LibraryEntry("One peak", 496.3398, "[M+H]+", "", [(184.0733, 100.0)]),
        # Predicted, although the library is declared experimental.
        LibraryEntry("LPC 16:0", 496.3398, "[M+H]+", "", [(104.107, 30.0), (184.0733, 100.0)],
                     comment="in-silico MSMS by LipidBlast"),
        # Negative mode: removed.
        LibraryEntry("PC 34:1", 804.5760, "[M+HCOO]-", "", [(255.233, 100.0), (281.2486, 80.0)], ionmode="Negative"),
    ])
    standards = write_msp(tmp_path / "standards.msp", [LibraryEntry("PC 16:0_18:1", 760.5851, "[M+H]+", "", pc_peaks, rt_min=3.0)])

    table, report = prepare_library(messy, "positive", 17.0)
    assert list(table["name"]) == ["LPC 16:0", "POPC"]  # sorted by precursor m/z
    popc = table.set_index("name").loc["POPC"]
    assert popc["adduct"] == "[M+H]+" and popc["formula"] == "C42H82NO8P"  # formula derived from the SMILES
    assert popc["inchikey"].startswith("WTJKGGKOPKCXLL")
    assert table.set_index("name")["in_silico"].to_dict() == {"LPC 16:0": True, "POPC": False}
    removed = report.set_index("step")["removed spectra"]
    assert removed["keep_ion_mode"] == 1 and removed["repair_structure_annotation"] == 1
    assert removed["require_minimum_number_of_peaks"] == 1 and (report["spectra read"] == 5).all()

    # The standard's copy of PC is kept, the messy library's duplicate removed.
    std_table, _ = prepare_library(standards, "positive", 17.0)
    entries = [{"name": "messy", "kind": "experimental", "reference_standards": False},
               {"name": "standards", "kind": "experimental", "reference_standards": True}]
    combined, summary = combine_libraries([table, std_table], entries)
    assert sorted(zip(combined["library"], combined["name"], combined["kind"])) == [
        ("messy", "LPC 16:0", "in_silico"), ("standards", "PC 16:0_18:1", "experimental")]
    assert summary.set_index("library")["duplicates removed"].to_dict() == {"messy": 1, "standards": 0}

    # The harmonized library as MGF, read back as a library (RTINSECONDS: rt_unit "s").
    mgf = tmp_path / "harmonized.mgf"
    write_library_mgf(combined, mgf)
    text = mgf.read_text()
    assert text.count("BEGIN IONS") == 2 and "LIBRARY=standards" in text and "KIND=in_silico" in text
    again, _ = prepare_library(mgf, "positive", 17.0, rt_unit="s")
    assert sorted(again["name"]) == sorted(combined["name"])
    assert again.set_index("name").loc["PC 16:0_18:1", "rt_s"] == pytest.approx(180.0)


def test_json_libraries(tmp_path):
    """MoNA and MassBank JSON (matchms reads neither), GNPS JSON, streamed in small chunks."""
    from atlas_ms.annotation.library_files import json_records
    from atlas_ms.annotation.libraries import prepare_library
    from synthetic import LibraryEntry, massbank_record, mona_record, write_json

    pc_peaks = [(104.107, 5.0), (184.0733, 100.0), (577.519, 3.0)]
    pc = LibraryEntry("PC 16:0_18:1", 760.5851, "[M+H]+", "", pc_peaks, rt_min=3.0, smiles=POPC)
    lpc = LibraryEntry("LPC 16:0", 496.3398, "[M+H]+", "C24H50NO7P", [(104.107, 30.0), (184.0733, 100.0)])
    mona = write_json(tmp_path / "mona.json", [mona_record(pc), mona_record(lpc, tags=("LipidBlast", "In-Silico"))])

    # The records come back whole, whatever the chunk size.
    assert [r["id"] for r in json_records(mona, chunk_chars=50)] == [r["id"] for r in json_records(mona)]

    table, report = prepare_library(mona, "positive", 17.0)
    rows = table.set_index("name")
    assert rows["in_silico"].to_dict() == {"LPC 16:0": True, "PC 16:0_18:1": False}  # from the MoNA tags
    assert rows.loc["PC 16:0_18:1", "formula"] == "C42H82NO8P"  # derived from the SMILES
    assert rows.loc["PC 16:0_18:1", "rt_s"] == pytest.approx(180.0)  # "3.0 min"
    assert rows.loc["PC 16:0_18:1", "adduct"] == "[M+H]+"
    assert report.set_index("step").loc["mark_in_silico", "changed spectra"] == 1

    massbank = write_json(tmp_path / "MassBank.json", [
        massbank_record(pc, "MSBNK-test-00001"),
        massbank_record(lpc, "MSBNK-test-00002", ms_type="MS"),  # an MS1 spectrum: removed
        {"ACCESSION": "MSBNK-test-00003", "DEPRECATED": True, "DEPRECATED_CONTENT": "..."},  # skipped
        # A GNPS record in the same file (its own JSON format).
        {"Compound_Name": "LPC 16:0", "Precursor_MZ": "496.3398", "Adduct": "M+H", "Ion_Mode": "Positive",
         "peaks_json": "[[104.107, 30.0], [184.0733, 100.0]]", "spectrum_id": "CCMSLIB00000000001"},
    ])
    table, report = prepare_library(massbank, "positive", 17.0)
    assert list(table["library_id"]) == ["CCMSLIB00000000001", "MSBNK-test-00001"]
    assert table.loc[1, "rt_s"] == pytest.approx(180.0)  # "180.0 sec"
    assert report.set_index("step").loc["keep_ms2", "removed spectra"] == 1
    assert (report["spectra read"] == 3).all()

    # A file with no readable spectrum is an error, not an empty library.
    with pytest.raises(ValueError, match="No spectrum could be read"):
        prepare_library(write_json(tmp_path / "other.json", [{"something": "else"}]), "positive", 17.0)
