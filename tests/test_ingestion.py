import pandas as pd

from kamr.ingestion.bvbrc_amr import build_query
from kamr.ingestion.summarize import drug_inventory

LAB = "Laboratory Method"


def _rec(genome_id, antibiotic, phenotype, evidence=LAB):
    return {
        "genome_id": genome_id,
        "antibiotic": antibiotic,
        "resistant_phenotype": phenotype,
        "evidence": evidence,
    }


def test_build_query_pages():
    q = build_query(573, 100, 200)
    assert q.startswith("eq(taxon_id,573)")
    assert q.endswith("limit(100,200)")


def test_drug_inventory_lab_only_and_conflicts():
    df = pd.DataFrame(
        [
            _rec("1", "Meropenem", "Resistant"),
            _rec("2", "meropenem ", "Susceptible"),
            _rec("2", "meropenem", "Resistant"),
            _rec("3", "meropenem", "Resistant", evidence="Computational Method"),
        ]
    )
    row = drug_inventory(df).loc["meropenem"]
    assert row["genomes"] == 2
    assert row["Resistant"] == 2
    assert row["Susceptible"] == 1
    assert row["conflicting"] == 1
