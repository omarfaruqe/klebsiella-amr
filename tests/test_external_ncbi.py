import json

import pandas as pd
import pytest

from kamr.external.ncbi_ast import (
    PHENOTYPE_MAP,
    cohort_identifiers,
    harmonise,
    read_export,
    strip_version,
)

BP = """antibiotic,s_le_mg_per_l,r_gt_mg_per_l,applies_to_kpneumoniae,standard,version,verified
colistin,2,2,yes,EUCAST,v16.1,yes
"""


@pytest.fixture
def breakpoints(tmp_path):
    p = tmp_path / "bp.csv"
    p.write_text(BP)
    return p


def frame(**over):
    base = {
        "biosample": ["SAMN1"],
        "assembly": ["GCA_000000001.1"],
        "bioproject": ["PRJNA1"],
        "collection_date": ["2020"],
        "location": ["Canada"],
        "submitted_phenotype": ["resistant"],
        "mic": ["8"],
        "sign": ["=="],
        "unit": ["mg/L"],
        "standard": ["CLSI"],
    }
    base.update({k: [v] for k, v in over.items()})
    return pd.DataFrame(base)


def test_strip_version_ignores_assembly_version():
    assert strip_version("GCA_000240185.2") == strip_version("GCA_000240185.1")
    assert strip_version("GCA_000240185.2") == "GCA_000240185"


def test_mic_overrides_submitted_phenotype(breakpoints):
    # submitted as susceptible under CLSI, but 8 mg/L is resistant under EUCAST v16.1
    out = harmonise(frame(submitted_phenotype="susceptible"), breakpoints)
    assert out["label"].iloc[0] == "R"
    assert out["label_source"].iloc[0] == "mic_eucast"
    assert out["breakpoint"].iloc[0] == "EUCAST v16.1"


def test_submitted_phenotype_used_only_without_usable_mic(breakpoints):
    out = harmonise(frame(mic="", unit="", submitted_phenotype="resistant"), breakpoints)
    assert out["label"].iloc[0] == "R"
    assert out["label_source"].iloc[0] == "submitted"


def test_censored_mic_that_cannot_decide_falls_back(breakpoints):
    # "<=4" is consistent with both S and R when s_le is 2, so the submitted call is used
    out = harmonise(frame(mic="4", sign="<=", submitted_phenotype="susceptible"), breakpoints)
    assert out["label_source"].iloc[0] == "submitted"
    assert out["label"].iloc[0] == "S"


def test_open_ended_mic_above_the_breakpoint_is_still_decisive(breakpoints):
    # ">2" establishes R outright, because R is defined as MIC > 2; it must not fall back
    out = harmonise(frame(mic="2", sign=">", submitted_phenotype="susceptible"), breakpoints)
    assert out["label_source"].iloc[0] == "mic_eucast"
    assert out["label"].iloc[0] == "R"


def test_nonsusceptible_is_not_silently_called_resistant(breakpoints):
    assert "nonsusceptible" not in PHENOTYPE_MAP
    out = harmonise(frame(mic="", unit="", submitted_phenotype="nonsusceptible"), breakpoints)
    assert out["label"].iloc[0] == ""
    assert out["label_source"].iloc[0] == "unusable"


def test_unrecognised_unit_does_not_reach_the_breakpoint(breakpoints):
    # a disk diameter in mm must never be read as an MIC in mg/L
    out = harmonise(frame(mic="14", unit="mm", submitted_phenotype="susceptible"), breakpoints)
    assert out["label_source"].iloc[0] == "submitted"


def test_read_export_normalises_column_spelling(tmp_path):
    p = tmp_path / "e.tsv"
    p.write_text(
        "BioSample\tasm_acc\tMIC (mg/L)\tMeasurement sign\tTesting standard\n"
        "SAMN9\tGCA_9.1\t4\t==\tEUCAST\n"
    )
    out = read_export(p)
    assert out["biosample"].iloc[0] == "SAMN9"
    assert out["assembly"].iloc[0] == "GCA_9.1"
    assert out["mic"].iloc[0] == "4"
    assert out["standard"].iloc[0] == "EUCAST"


def test_read_export_rejects_an_export_without_assembly_accessions(tmp_path):
    p = tmp_path / "e.csv"
    p.write_text("biosample_acc,asm_acc,mic\nSAMN9,,4\n")
    with pytest.raises(SystemExit, match="assembly accessions"):
        read_export(p)


def test_cohort_identifiers_reads_the_newest_snapshot(tmp_path):
    for day, sample in (("2026-01-01", "SAMN_OLD"), ("2026-02-01", "SAMN_NEW")):
        d = tmp_path / day
        d.mkdir()
        (d / "genome.jsonl").write_text(
            json.dumps({"biosample_accession": sample, "assembly_accession": "GCA_1.1"}) + "\n"
        )
    biosamples, assemblies = cohort_identifiers(tmp_path)
    assert biosamples == {"SAMN_NEW"}
    assert assemblies == {"GCA_1.1"}


def test_cohort_identifiers_requires_a_snapshot(tmp_path):
    with pytest.raises(SystemExit, match="no BV-BRC genome snapshot"):
        cohort_identifiers(tmp_path / "absent")


def test_export_without_units_column_treats_mic_as_mg_per_litre(tmp_path, breakpoints):
    # NCBI's ast table has no units column: mic is mg/L by definition. Without this, every MIC
    # would fail the unit guard and fall back to the submitter's interpretation, silently
    # defeating the breakpoint harmonisation.
    p = tmp_path / "e.csv"
    p.write_text(
        "biosample_acc,asm_acc,mic,measurement_sign,phenotype,standard\n"
        "SAMN1,GCA_1.1,8,==,susceptible,CLSI\n"
    )
    out = read_export(p)
    assert out["unit"].iloc[0] == "mg/L"
    harmonised = harmonise(out, breakpoints)
    assert harmonised["label_source"].iloc[0] == "mic_eucast"
    assert harmonised["label"].iloc[0] == "R"


def test_explicit_units_are_not_overridden(tmp_path, breakpoints):
    # a real units column must win, so a disk diameter is never read as an MIC
    p = tmp_path / "e.csv"
    p.write_text(
        "biosample_acc,asm_acc,mic,measurement_units,measurement_sign,phenotype\n"
        "SAMN1,GCA_1.1,14,mm,==,susceptible\n"
    )
    out = read_export(p)
    assert out["unit"].iloc[0] == "mm"
    assert harmonise(out, breakpoints)["label_source"].iloc[0] == "submitted"
