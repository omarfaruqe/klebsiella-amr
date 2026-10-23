import json

import pandas as pd

from kamr.curation.breakpoints import Breakpoint
from kamr.curation.cohort import aggregate_labels, annotate_records, read_lab_records

BPS = {
    "meropenem": Breakpoint("meropenem", 2, 8, True, "test"),
    "tigecycline": Breakpoint("tigecycline", 0.5, 0.5, False, "test"),
}


def rec(gid, drug, pheno=None, sign="=", value=None, unit="mg/L", evidence="Laboratory Method"):
    return {
        "id": f"{gid}-{drug}-{pheno}-{value}",
        "genome_id": gid,
        "antibiotic": drug,
        "resistant_phenotype": pheno,
        "measurement_sign": sign,
        "measurement_value": value,
        "measurement_unit": unit,
        "laboratory_typing_method": "Broth dilution",
        "testing_standard": "CLSI",
        "testing_standard_year": None,
        "evidence": evidence,
    }


def test_read_lab_records_filters(tmp_path):
    path = tmp_path / "snap.jsonl"
    rows = [
        rec("1", "Meropenem", "Resistant", value="16"),
        rec("2", "meropenem", "Resistant", value="16", evidence="Computational Method"),
        rec("3", "ampicillin", "Resistant", value="32"),
    ]
    path.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    df = read_lab_records(path, drugs=("meropenem",))
    assert list(df["genome_id"]) == ["1"]
    assert df.loc[0, "antibiotic"] == "meropenem"


def labels_for(rows):
    df = pd.DataFrame(rows).drop(columns=["evidence"])
    out = aggregate_labels(annotate_records(df, BPS))
    return out.set_index(["genome_id", "antibiotic"])


def test_unknown_with_mic_is_recovered():
    out = labels_for([rec("1", "meropenem", None, value="16")])
    row = out.loc[("1", "meropenem")]
    assert row["label"] == "R"
    assert bool(row["recovered_from_unknown"])


def test_reinterpretation_overrides_reported_call():
    out = labels_for([rec("1", "meropenem", "Resistant", value="4")])
    assert out.loc[("1", "meropenem"), "label"] == "I"


def test_conflict_and_intermediate_rules():
    out = labels_for(
        [
            rec("1", "meropenem", None, value="0.5"),
            rec("1", "meropenem", None, value="16"),
            rec("2", "meropenem", None, value="4"),
            rec("2", "meropenem", None, value="16"),
        ]
    )
    assert out.loc[("1", "meropenem"), "label"] == "conflict"
    assert out.loc[("2", "meropenem"), "label"] == "R"


def test_breakpoint_not_applicable_uses_reported_call():
    out = labels_for([rec("1", "tigecycline", "Susceptible", value="4")])
    row = out.loc[("1", "tigecycline")]
    assert row["label"] == "S"
    assert row["label_sources"] == "reported"


def test_disk_diffusion_without_call_is_unlabelled():
    out = labels_for([rec("1", "meropenem", None, value="22", unit="mm")])
    assert out.loc[("1", "meropenem"), "label"] == "unlabelled"
