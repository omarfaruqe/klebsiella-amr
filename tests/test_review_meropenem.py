import pandas as pd

from kamr.analysis.review_meropenem import (
    build_table,
    carbapenemase_family,
    family_summary,
    genome_mechanisms,
    md_table,
)


def row(gid, sym, sub, subtype="AMR", etype="AMR"):
    return {
        "genome_id": gid,
        "Element symbol": sym,
        "Type": etype,
        "Subtype": subtype,
        "Class": "BETA-LACTAM",
        "Subclass": sub,
    }


AMR = pd.DataFrame(
    [
        row("1", "blaKPC-2", "CARBAPENEM"),
        row("2", "blaOXA-48", "CARBAPENEM"),
        row("3", "blaOXA-232", "CARBAPENEM"),
        row("3", "blaNDM-1", "CARBAPENEM"),
        row("4", "blaCTX-M-15", "CEPHALOSPORIN"),
        row("4", "ompK36_D135DGD", "CARBAPENEM", subtype="POINT"),
        row("5", "blaSHV-11", "BETA-LACTAM"),
    ]
)


def test_family_names():
    assert carbapenemase_family("blaKPC-3") == "KPC"
    assert carbapenemase_family("blaOXA-181") == "OXA-48-like"
    assert carbapenemase_family("blaOXA-23") == "other"
    assert carbapenemase_family("blaOXA-4811") == "other"
    assert carbapenemase_family("blaVIM-1") == "VIM"


def test_mechanisms():
    m = genome_mechanisms(AMR)
    assert m.loc["3", "carbapenemase_family"] == "NDM;OXA-48-like"
    assert m.loc["4", "carbapenemase_family"] == ""
    assert m.loc["4", "esbl_ampc"] == "blaCTX-M-15"
    assert m.loc["4", "porin_defects"] == "ompK36_D135DGD"


def test_table_and_summary():
    labels = pd.DataFrame(
        {"genome_id": list("12345"), "antibiotic": "meropenem", "label": list("RSRRS")}
    )
    records = pd.DataFrame(
        {
            "genome_id": ["2", "2", "4"],
            "antibiotic": "meropenem",
            "eucast_call": ["S", "S", "R"],
            "measurement_value": ["0.5", "2", "16"],
            "measurement_sign": ["=", "<=", ">"],
            "laboratory_typing_method": "Broth dilution",
            "testing_standard": "CLSI",
        }
    )
    qc = pd.DataFrame(
        {
            "genome_id": list("12345"),
            "st": "ST11",
            "bioproject_accession": "P",
            "collection_year": 2015,
            "isolation_country": "USA",
        }
    )
    t = build_table(AMR, labels, records, qc, set("12345")).set_index("genome_id")
    assert t.loc["2", "mic_text"] == "<=2"
    assert t.loc["4", "mic_text"] == ">16"
    assert t.loc["3", "family_group"] == "multiple"
    assert not bool(t.loc["4", "has_carbapenemase"])
    s = family_summary(t.reset_index())
    assert s.loc["OXA-48-like", "S"] == 1
    assert s.loc["KPC", "R"] == 1


def test_md_table_formats_counts_and_shares():
    df = pd.DataFrame({"genomes": [58, 29], "S_share": [0.517, 0.69]}, index=["OXA-48-like", "KPC"])
    lines = md_table(df, "Family")
    assert lines[2] == "| OXA-48-like | 58 | 51.7% |"
    assert lines[3] == "| KPC | 29 | 69.0% |"
