import pandas as pd
import pytest

from kamr.benchmark.features import (
    amr_features,
    borderline,
    feature_matrix,
    gene_family,
    label_table,
    mgrb_features,
)
from kamr.curation.breakpoints import Breakpoint


@pytest.mark.parametrize(
    ("symbol", "family"),
    [
        ("blaKPC-2", "blaKPC"),
        ("blaCTX-M-15", "blaCTX-M"),
        ("rmtB1", "rmtB"),
        ("qnrS1", "qnrS"),
        ("aac(6')-Ib-cr5", "aac(6')-Ib-cr"),
        ("aph(3')-Ia", "aph(3')-Ia"),
        ("tet(X4)", "tet(X4)"),
        ("armA", "armA"),
    ],
)
def test_gene_family(symbol, family):
    assert gene_family(symbol) == family


def test_features_and_matrix():
    amr = pd.DataFrame(
        {
            "genome_id": ["1", "1", "1", "2", "2"],
            "Element symbol": ["blaKPC-2", "gyrA_S83I", "oqxA", "blaKPC-3", "fieF"],
            "Type": ["AMR", "AMR", "AMR", "AMR", "STRESS"],
            "Subtype": ["AMR", "POINT", "AMR", "AMR", "STRESS"],
        }
    )
    status = pd.DataFrame(
        {
            "genome_id": ["1", "2"],
            "mgrb_status": ["split_contigs", "missense"],
            "mgrb_disrupted": ["True", "False"],
        }
    )
    long = pd.concat([amr_features(amr), mgrb_features(status)])
    x = feature_matrix(long, ["1", "2", "3"], min_count=1)
    assert list(x.columns) == ["blaKPC", "gyrA_S83I", "mgrB_disrupted", "mgrB_missense"]
    assert x.loc["2"].tolist() == [1, 0, 0, 1]
    assert x.loc["3"].sum() == 0
    assert list(feature_matrix(long, ["1", "2"], min_count=2).columns) == ["blaKPC"]


def test_borderline():
    assert borderline(2, True, 2, 2) and borderline(4, True, 2, 2)
    assert not borderline(8, True, 2, 2)
    assert not borderline(4, False, 2, 2)
    assert borderline(16, True, 2, 8) and not borderline(4, True, 2, 8)


def test_label_table_flags_and_borderline():
    labels = pd.DataFrame(
        {"genome_id": ["1", "2", "3"], "antibiotic": "gentamicin", "label": ["R", "S", "R"]}
    )
    records = pd.DataFrame(
        {
            "genome_id": ["1", "2", "3"],
            "antibiotic": "gentamicin",
            "eucast_call": ["R", "S", "R"],
            "record_label": ["R", "S", "R"],
            "measurement_value": ["4", "1", "32"],
            "measurement_sign": ["=", "=", ">"],
            "laboratory_typing_method": "Broth dilution",
        }
    )
    qc = pd.DataFrame(
        {
            "genome_id": ["1", "2", "3"],
            "bioproject_accession": ["P1", "P2", "P1"],
            "st": "ST11",
            "collection_year": [2012, 2015, None],
        }
    )
    flags = pd.DataFrame(
        {"bioproject_accession": ["P2"], "antibiotic": ["gentamicin"], "action": ["exclude"]}
    )
    bps = {"gentamicin": Breakpoint("gentamicin", 2, 2, True, "t")}
    t = label_table(labels, records, qc, flags, bps, ["1", "2", "3"]).set_index("genome_id")
    assert bool(t.loc["1", "borderline"]) and not bool(t.loc["3", "borderline"])
    assert bool(t.loc["2", "flagged"]) and not bool(t.loc["1", "flagged"])
    assert t["y"].tolist() == [1, 0, 1]
