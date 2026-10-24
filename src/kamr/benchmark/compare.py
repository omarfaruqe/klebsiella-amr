"""Paired bootstrap comparison of predictors from held-out predictions.

The benchmark reports point estimates of balanced accuracy per drug, split scheme and model.
This module tests whether the differences between predictors are larger than sampling noise, by
resampling genomes with replacement and recomputing each pairwise difference on the same
resampled set (paired, so the comparison is not inflated by between-resample variation).

Input:  data/features/oof_predictions.parquet (written by kamr.benchmark.run)
Output: reports/benchmark/comparisons.md, reports/benchmark/comparisons.csv
"""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20261024
MIN_PER_CLASS = 10  # must match kamr.benchmark.run; degenerate splits are not comparable
PAIRS = [("rule", "logreg"), ("rule", "gboost"), ("logreg", "gboost")]


def balanced_accuracy(y: np.ndarray, pred: np.ndarray) -> float:
    pos, neg = y == 1, y == 0
    if not pos.any() or not neg.any():
        return np.nan
    return (pred[pos].mean() + (1 - pred[neg]).mean()) / 2


def paired_bootstrap(
    wide: pd.DataFrame, a: str, b: str, n_boot: int, rng: np.random.Generator
) -> dict:
    y = wide["y"].to_numpy()
    pa, pb = wide[a].to_numpy(), wide[b].to_numpy()
    obs = balanced_accuracy(y, pa) - balanced_accuracy(y, pb)
    idx = rng.integers(0, len(y), size=(n_boot, len(y)))
    diffs = np.array([balanced_accuracy(y[i], pa[i]) - balanced_accuracy(y[i], pb[i]) for i in idx])
    diffs = diffs[~np.isnan(diffs)]
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return {
        "model_a": a,
        "model_b": b,
        "diff": obs,
        "lo": lo,
        "hi": hi,
        "excludes_zero": bool(lo > 0 or hi < 0),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--oof", type=Path, default=Path("data/features/oof_predictions.parquet"))
    p.add_argument("--out", type=Path, default=Path("reports/benchmark"))
    p.add_argument("--n-boot", type=int, default=2000)
    a = p.parse_args()

    oof = pd.read_parquet(a.oof)
    rng = np.random.default_rng(SEED)
    rows = []
    for (drug, scheme), g in oof.groupby(["drug", "scheme"]):
        wide = g.pivot_table(
            index=["genome_id", "y"], columns="model", values="pred", aggfunc="first"
        ).reset_index()
        if wide[["rule", "logreg", "gboost"]].isna().any().any() or wide["y"].nunique() < 2:
            continue
        if min(int(wide.y.sum()), int((1 - wide.y).sum())) < MIN_PER_CLASS:
            continue  # same guard as the benchmark: too few of one class to compare
        for x, y_ in PAIRS:
            rows.append(
                {
                    "drug": drug,
                    "scheme": scheme,
                    "n": len(wide),
                    **paired_bootstrap(wide, x, y_, a.n_boot, rng),
                }
            )
    res = pd.DataFrame(rows)
    a.out.mkdir(parents=True, exist_ok=True)
    res.to_csv(a.out / "comparisons.csv", index=False, float_format="%.4f")

    L = [
        "# Paired comparisons between predictors",
        "",
        f"- Generated (UTC): {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        f"- {a.n_boot:,} paired bootstrap resamples of held-out genomes, seed {SEED}",
        "- Difference is balanced accuracy of A minus B; interval is the 2.5–97.5 percentile",
        "- 'Yes' means the 95% interval excludes zero",
        f"- Splits with fewer than {MIN_PER_CLASS} genomes of either class are omitted",
        "",
        "| Drug | Split | Comparison | Difference | 95% interval | Distinguishable |",
        "|---|---|---|---|---|---|",
    ]
    for r in res.itertuples():
        L.append(
            f"| {r.drug} | {r.scheme} | {r.model_a} vs {r.model_b} | {r.diff:+.3f} | "
            f"{r.lo:+.3f} to {r.hi:+.3f} | {'Yes' if r.excludes_zero else 'No'} |"
        )
    (a.out / "comparisons.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
