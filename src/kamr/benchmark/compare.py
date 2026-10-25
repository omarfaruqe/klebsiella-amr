"""Paired bootstrap comparison of predictors from held-out predictions.

The benchmark reports point estimates of balanced accuracy per drug, split scheme and model.
This module tests whether the differences between predictors are larger than sampling noise, by
resampling genomes with replacement and recomputing each pairwise difference on the same
resampled set (paired, so the comparison is not inflated by between-resample variation).

Every drug, split scheme and predictor pair is compared, which is a family of tests rather than
one, so some intervals will exclude zero by chance alone. Reporting only the per-comparison
intervals would overstate what the benchmark establishes. Each comparison therefore also carries
a two-sided bootstrap p-value and a Benjamini-Hochberg q-value computed across the whole family,
together with whether it would survive the far stricter Bonferroni threshold.
Benjamini-Hochberg is the primary correction: these comparisons are positively dependent (they
share predictors, genomes and folds), and controlling the false discovery rate is the appropriate
goal for a benchmark that is mapping where differences exist rather than testing one hypothesis.

The p-value resolution is bounded by the number of resamples: with B resamples no p-value below
1/B is distinguishable, and such values are reported at that floor rather than as zero.

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


def bootstrap_pvalue(diffs: np.ndarray) -> float:
    """Two-sided bootstrap p-value for H0: the difference is zero.

    Twice the smaller tail mass on either side of zero, floored at 1/B because no finer
    resolution exists with B resamples. Returns 1.0 for an empty distribution.
    """
    if not len(diffs):
        return 1.0
    tail = min(float(np.mean(diffs <= 0)), float(np.mean(diffs >= 0)))
    return float(min(1.0, max(2.0 * tail, 1.0 / len(diffs))))


def bh_qvalues(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values, order preserved."""
    p = np.asarray(pvalues, dtype=float)
    m = len(p)
    if not m:
        return p
    order = np.argsort(p)
    ranked = p[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]  # enforce monotonicity
    q = np.empty(m, dtype=float)
    q[order] = np.minimum(ranked, 1.0)
    return q


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
        "p": bootstrap_pvalue(diffs),
        "n_boot": len(diffs),
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
    res["q"] = bh_qvalues(res["p"].to_numpy())
    bonferroni = 0.05 / len(res) if len(res) else float("nan")
    res["survives_bh"] = res["q"] < 0.05
    res["survives_bonferroni"] = res["p"] < bonferroni
    a.out.mkdir(parents=True, exist_ok=True)
    res.to_csv(a.out / "comparisons.csv", index=False, float_format="%.4f")

    L = [
        "# Paired comparisons between predictors",
        "",
        f"- Generated (UTC): {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        f"- {a.n_boot:,} paired bootstrap resamples of held-out genomes, seed {SEED}",
        "- Difference is balanced accuracy of A minus B; interval is the 2.5–97.5 percentile",
        f"- Splits with fewer than {MIN_PER_CLASS} genomes of either class are omitted",
        f"- Family of {len(res)} comparisons; q is the Benjamini-Hochberg adjusted p-value",
        f"- Bonferroni threshold for this family is p < {bonferroni:.5f}",
        f"- p is floored at 1/{a.n_boot:,} = {1 / a.n_boot:.5f}, the resolution of the bootstrap",
        "",
        f"Uncorrected, {int(res['excludes_zero'].sum())} of {len(res)} intervals exclude zero; "
        f"{int(res['survives_bh'].sum())} survive Benjamini-Hochberg at 5% and "
        f"{int(res['survives_bonferroni'].sum())} survive Bonferroni.",
        "",
        "| Drug | Split | Comparison | Difference | 95% interval | p | BH q | Survives BH |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in res.itertuples():
        pstr = f"<{1 / a.n_boot:.4f}" if r.p <= 1 / a.n_boot else f"{r.p:.4f}"
        L.append(
            f"| {r.drug} | {r.scheme} | {r.model_a} vs {r.model_b} | {r.diff:+.3f} | "
            f"{r.lo:+.3f} to {r.hi:+.3f} | {pstr} | {r.q:.4f} | "
            f"{'Yes' if r.survives_bh else 'No'} |"
        )
    (a.out / "comparisons.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
