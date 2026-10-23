import json

import pandas as pd
import pytest

from kamr.curation.qc import apply_qc, load_metadata, parse_st, select_sample
from kamr.ingestion.bvbrc_genome import batches, build_query

CFG = {
    "assembly": {
        "species_pattern": "Klebsiella pneumoniae",
        "min_length": 5_000_000,
        "max_length": 6_500_000,
        "max_contigs": 500,
        "min_n50": 20_000,
        "min_checkm_completeness": 95.0,
        "max_checkm_contamination": 5.0,
    }
}


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("MLST.Klebsiella_pneumoniae.11", "ST11"),
        (["MLST.Klebsiella_pneumoniae.258"], "ST258"),
        ("ST147", "ST147"),
        (None, None),
        ("novel", None),
    ],
)
def test_parse_st(value, expected):
    assert parse_st(value) == expected


def test_build_query_and_batches():
    q = build_query(["573.1", "573.2"])
    assert q.startswith("in(genome_id,(573.1,573.2))")
    assert q.endswith("limit(12)")
    assert [len(b) for b in batches(list(range(450)), 200)] == [200, 200, 50]


def genome(gid, **kw):
    base = {
        "genome_id": gid,
        "genome_name": "Klebsiella pneumoniae strain X",
        "species": "Klebsiella pneumoniae",
        "contigs": 80,
        "contig_n50": 200_000,
        "genome_length": 5_500_000,
        "checkm_completeness": 99.5,
        "checkm_contamination": 0.5,
        "mlst": "MLST.Klebsiella_pneumoniae.11",
        "assembly_accession": f"GCA_{gid}",
    }
    base.update(kw)
    return base


def meta_frame(tmp_path, rows):
    path = tmp_path / "genome.jsonl"
    path.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    return load_metadata(path)


def test_qc_filters_and_funnel(tmp_path):
    meta = meta_frame(
        tmp_path,
        [
            genome("1"),
            genome("2", species="Klebsiella variicola", genome_name="Klebsiella variicola Y"),
            genome("3", genome_length=4_000_000),
            genome("4", contigs=900),
            genome("5", checkm_completeness=None, checkm_contamination=None),
            genome("6", checkm_contamination=9.0),
            genome("1"),
        ],
    )
    qc, funnel = apply_qc(meta, CFG)
    reasons = qc.set_index("genome_id")["qc_fail_reason"].to_dict()
    assert len(qc) == 6
    assert reasons["2"] == "species"
    assert reasons["3"] == "min_length"
    assert reasons["4"] == "max_contigs"
    assert reasons["6"] == "checkm_contamination"
    assert set(qc.loc[qc["qc_pass"], "genome_id"]) == {"1", "5"}
    assert funnel.iloc[-1]["remaining"] == 2
    assert int(funnel["removed"].sum()) == 4


def test_select_sample_covers_strata_and_spreads_sts():
    qc = pd.DataFrame(
        {
            "genome_id": [str(i) for i in range(40)],
            "st": ["ST258"] * 30 + [f"ST{i}" for i in range(10)],
            "qc_pass": [True] * 39 + [False],
        }
    )
    labels = pd.DataFrame(
        {
            "genome_id": [str(i) for i in range(40)],
            "antibiotic": "meropenem",
            "label": ["R"] * 20 + ["S"] * 20,
        }
    )
    a = select_sample(qc, labels, ["meropenem"], per_stratum=5, seed=1)
    b = select_sample(qc, labels, ["meropenem"], per_stratum=5, seed=1)
    pd.testing.assert_frame_equal(a, b)
    assert (a["meropenem"] == "R").sum() == 5
    assert (a["meropenem"] == "S").sum() == 5
    assert "39" not in set(a["genome_id"])
    s_rows = a[a["meropenem"] == "S"]
    assert s_rows["st"].nunique() == 5
