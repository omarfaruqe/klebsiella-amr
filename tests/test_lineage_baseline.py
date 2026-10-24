import numpy as np
import pandas as pd

from kamr.benchmark.lineage_baseline import lineage_only_predictions, overlap_fraction


def frame(n=200):
    # resistance perfectly determined by ST: ST0/ST1 resistant, ST2/ST3 susceptible
    st = [f"ST{i % 4}" for i in range(n)]
    y = [1 if int(s[2:]) < 2 else 0 for s in st]
    return pd.DataFrame(
        {
            "genome_id": [f"g{i}" for i in range(n)],
            "y": y,
            "st": st,
            "bioproject_accession": [f"P{i % 5}" for i in range(n)],
            "collection_year": [2010 + i % 10 for i in range(n)],
        }
    )


def test_lineage_only_is_perfect_under_random_and_chance_under_lineage_split():
    lab = frame()
    pred_r, seen_r = lineage_only_predictions(lab, "random")
    y = lab["y"].to_numpy()
    assert seen_r.all()  # every test ST appears in training
    assert (pred_r == y).mean() == 1.0  # ST fully determines the label here

    pred_l, seen_l = lineage_only_predictions(lab, "lineage")
    assert not seen_l.any()  # no test ST appears in training, by construction
    assert set(np.unique(pred_l)) <= {0, 1}  # falls back to the training majority


def test_overlap_fraction():
    lab = frame()
    assert overlap_fraction(lab, "random", "st") == 1.0
    assert overlap_fraction(lab, "lineage", "st") == 0.0
    assert overlap_fraction(lab, "study", "bioproject_accession") == 0.0
