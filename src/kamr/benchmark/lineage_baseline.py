"""Measure how much of model performance is attributable to lineage rather than mechanism.

The main benchmark shows model performance falling when sequence types are held out, which
suggests the models were partly exploiting lineage. That is an inference from behaviour. This
module measures the quantity directly in two ways, for each drug and split scheme:

1. Overlap. What fraction of test genomes share a sequence type, or a BioProject, with at least
   one training genome? This is the leakage itself, before any model is involved.
2. A lineage-only predictor. For each test genome, predict the majority label of its sequence
   type among training genomes, falling back to the overall training majority when the type is
   unseen. It uses no resistance determinant of any kind. Whatever accuracy it achieves is
   accuracy obtainable from lineage alone.

If the models' advantage under random splitting came from lineage, the lineage-only predictor
should score well there and collapse to near chance under lineage-aware splitting, where every
test type is unseen by construction.

Outputs: reports/benchmark/lineage_baseline.md, reports/benchmark/lineage_baseline.csv
"""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import numpy as np
import pandas as pd

from kamr.benchmark.compare import balanced_accuracy
from kamr.benchmark.features import DRUGS
from kamr.benchmark.splits import SCHEMES, folds, groups_for

MIN_PER_CLASS = 10  # must match kamr.benchmark.run


def lineage_only_predictions(lab: pd.DataFrame, scheme: str) -> tuple[np.ndarray, np.ndarray]:
    """Out-of-fold predictions from sequence type alone, plus whether each ST was seen."""
    y = lab["y"].to_numpy()
    st = groups_for(lab, "st")
    pred = np.full(len(lab), -1)
    seen = np.zeros(len(lab), dtype=bool)
    for train, test in folds(lab, scheme):
        majority = int(round(y[train].mean())) if len(train) else 0
        by_st: dict[str, int] = {}
        for group in set(st[train]):
            vals = y[train][st[train] == group]
            by_st[group] = int(vals.mean() >= 0.5)
        for i in test:
            pred[i] = by_st.get(st[i], majority)
            seen[i] = st[i] in by_st
    return pred, seen


def overlap_fraction(lab: pd.DataFrame, scheme: str, column: str) -> float:
    """Fraction of test genomes sharing a value of `column` with some training genome."""
    g = groups_for(lab, column)
    shared, total = 0, 0
    for train, test in folds(lab, scheme):
        train_groups = set(g[train])
        shared += sum(1 for i in test if g[i] in train_groups)
        total += len(test)
    return shared / total if total else float("nan")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--features", type=Path, default=Path("data/features"))
    p.add_argument("--benchmark", type=Path, default=Path("reports/benchmark/results.csv"))
    p.add_argument("--out", type=Path, default=Path("reports/benchmark"))
    a = p.parse_args()

    labels = pd.read_parquet(a.features / "labels.parquet")
    labels = labels[~labels["flagged"]]
    bench = pd.read_csv(a.benchmark)
    bench = bench[(bench.labels == "all") & bench.evaluable]

    rows = []
    for drug in DRUGS:
        lab = labels[labels["drug"] == drug].reset_index(drop=True)
        if lab.empty:
            continue
        for scheme in SCHEMES:
            pred, seen = lineage_only_predictions(lab, scheme)
            done = pred >= 0
            y = lab["y"].to_numpy()[done]
            if min(int(y.sum()), int((1 - y).sum())) < MIN_PER_CLASS:
                continue
            best = bench[
                (bench.drug == drug)
                & (bench.scheme == scheme)
                & (bench.model.isin(["logreg", "gboost"]))
            ]["balanced_accuracy"]
            rows.append(
                {
                    "drug": drug,
                    "scheme": scheme,
                    "n_test": int(done.sum()),
                    "st_overlap": overlap_fraction(lab, scheme, "st"),
                    "bioproject_overlap": overlap_fraction(lab, scheme, "bioproject_accession"),
                    "lineage_only_ba": balanced_accuracy(y, pred[done]),
                    "best_model_ba": float(best.max()) if len(best) else np.nan,
                }
            )
    res = pd.DataFrame(rows)
    a.out.mkdir(parents=True, exist_ok=True)
    res.to_csv(a.out / "lineage_baseline.csv", index=False)

    def fmt(v: float) -> str:
        return "n/a" if pd.isna(v) else f"{v:.3f}"

    L = [
        "# How much of model performance is lineage?",
        "",
        f"- Generated (UTC): {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        "- ST overlap: share of test genomes whose sequence type also appears in training",
        "- Lineage-only: predicts the majority label of the genome's sequence type in training,",
        "  falling back to the overall training majority for unseen types. Uses no resistance",
        "  determinant of any kind.",
        "- Best model: the better of logistic regression and gradient boosting from the benchmark",
        "",
        "| Drug | Split | ST overlap | BioProject overlap | Lineage-only | Best model |",
        "|---|---|---|---|---|---|",
    ]
    for r in res.itertuples():
        L.append(
            f"| {r.drug} | {r.scheme} | {r.st_overlap:.0%} | {r.bioproject_overlap:.0%} | "
            f"{fmt(r.lineage_only_ba)} | {fmt(r.best_model_ba)} |"
        )
    L += [
        "",
        "A lineage-only predictor that scores well under random splitting and near 0.5 under",
        "lineage-aware splitting demonstrates that lineage alone carries predictive signal, and",
        "that conventional validation lets models use it.",
    ]
    (a.out / "lineage_baseline.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
