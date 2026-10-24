import numpy as np
import pandas as pd

from kamr.benchmark.run import run_drug, score
from kamr.benchmark.splits import folds, groups_for, temporal_split


def lab_frame(n=60):
    rng = np.random.default_rng(0)
    return pd.DataFrame(
        {
            "genome_id": [f"g{i}" for i in range(n)],
            "y": rng.integers(0, 2, n),
            "bioproject_accession": [f"P{i % 6}" if i % 10 else None for i in range(n)],
            "st": [f"ST{i % 4}" for i in range(n)],
            "collection_year": [2005 + i % 15 if i % 7 else None for i in range(n)],
            "borderline": [i % 9 == 0 for i in range(n)],
        }
    )


def test_missing_groups_are_singletons():
    lab = lab_frame()
    g = groups_for(lab, "bioproject_accession")
    assert g[0] == "missing:g0" and g[10] == "missing:g10" and g[1] == "P1"


def test_grouped_folds_do_not_share_groups():
    lab = lab_frame()
    g = groups_for(lab, "bioproject_accession")
    for train, test in folds(lab, "study"):
        assert not set(g[train]) & set(g[test])
    covered = np.concatenate([t for _, t in folds(lab, "lineage")])
    assert sorted(covered) == list(range(len(lab)))


def test_temporal_split_is_ordered_and_dated_only():
    lab = lab_frame()
    train, test, cutoff = temporal_split(lab)
    years = lab["collection_year"]
    assert (years.iloc[train] <= cutoff).all() and (years.iloc[test] > cutoff).all()
    assert years.iloc[np.concatenate([train, test])].notna().all()


def test_score():
    y = np.array([1, 1, 1, 0, 0])
    pred = np.array([1, 1, 0, 0, 1])
    s = score(y, pred, np.array([0.9, 0.8, 0.4, 0.1, 0.6]))
    assert (s["TP"], s["FN"], s["FP"], s["TN"]) == (2, 1, 1, 1)
    assert abs(s["balanced_accuracy"] - (2 / 3 + 1 / 2) / 2) < 1e-9
    assert 0.5 < s["auroc"] <= 1


def test_run_drug_learns_a_marker():
    lab = lab_frame(200)
    x = pd.DataFrame({"marker": lab["y"].to_numpy(), "noise": 0}, index=lab["genome_id"])
    rules = pd.DataFrame({"meropenem": lab["y"].to_numpy()}, index=lab["genome_id"])
    rows, oof = run_drug(lab, x, rules, "meropenem")
    res = pd.DataFrame(rows)
    main = res[(res.labels == "all") & (res.scheme == "random")].set_index("model")
    assert main.loc["rule", "balanced_accuracy"] == 1.0
    assert main.loc["logreg", "balanced_accuracy"] > 0.95
    assert set(res.scheme) == {"random", "study", "lineage", "temporal"}
    assert set(oof.model) == {"rule", "logreg", "gboost"}


def test_small_test_sets_are_not_evaluable():
    lab = lab_frame(200)
    lab.loc[lab["collection_year"] > 2016, "y"] = 0
    x = pd.DataFrame({"marker": lab["y"].to_numpy()}, index=lab["genome_id"])
    rules = pd.DataFrame({"meropenem": lab["y"].to_numpy()}, index=lab["genome_id"])
    rows, _ = run_drug(lab, x, rules, "meropenem")
    res = pd.DataFrame(rows)
    temporal = res[(res.scheme == "temporal") & (res.labels == "all")]
    assert not temporal["evaluable"].any()
    assert temporal["balanced_accuracy"].isna().all()
    assert "not evaluable" in temporal["note"].iloc[0]
    assert res[res.scheme == "random"]["evaluable"].all()
