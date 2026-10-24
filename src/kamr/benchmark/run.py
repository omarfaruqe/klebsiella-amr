"""Leakage-aware benchmark: rules vs logistic regression vs gradient boosting, per drug and split.

For every drug and split scheme (random, study, lineage, temporal), out-of-fold predictions are
collected for three predictors:
  - rule: AMRFinderPlus core rule (rules.py); for colistin plus disrupted or missense mgrB;
  - logreg: L2 logistic regression, class_weight="balanced";
  - gboost: HistGradientBoostingClassifier, class_weight="balanced".
Models are trained on all non-flagged labels in the training fold. Metrics are computed on the
pooled held-out predictions, on all test labels and on non-borderline test labels (ADR-004).

Inputs:  data/features/X.parquet, data/features/labels.parquet, results/cohort/amrfinder.parquet,
         reports/review/mgrb_status.tsv
Outputs: reports/benchmark/results.csv, reports/benchmark/results.md,
         data/features/oof_predictions.parquet (not tracked)
"""

from __future__ import annotations

import argparse
import datetime as dt
import math
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from kamr.analysis.concordance import wilson
from kamr.analysis.rules import amrfinder_markers, predictions
from kamr.benchmark.features import DRUGS
from kamr.benchmark.splits import SCHEMES, folds, temporal_split

SEED = 20261024
MIN_PER_CLASS = 10  # a split is scored only if its test set has >= 10 R and >= 10 S


def make_model(name: str):
    if name == "logreg":
        return LogisticRegression(C=1.0, max_iter=2000, class_weight="balanced")
    if name == "gboost":
        return HistGradientBoostingClassifier(
            max_iter=200,
            learning_rate=0.1,
            max_leaf_nodes=15,
            l2_regularization=1.0,
            class_weight="balanced",
            random_state=SEED,
        )
    raise ValueError(name)


def rule_calls(amr: pd.DataFrame, mgrb: pd.DataFrame, genomes: list[str]) -> pd.DataFrame:
    """Genome x drug table of rule predictions (1 = R), colistin including the mgrB check.

    Only mechanism-specified mgrB disruption is used. The missense tier is excluded: a
    leave-one-study-out analysis (reports/review/mgrb_robustness.md) shows it adds 0.9 points of
    sensitivity out-of-sample against 19.1 points in-sample, i.e. it is almost entirely fitted to
    the cohort it was derived from.
    """
    p = predictions(amrfinder_markers(amr[amr.genome_id.isin(genomes)]), genomes, "rule")
    wide = p.assign(r=(p.predicted == "R").astype("uint8")).pivot(
        index="genome_id", columns="drug", values="r"
    )
    disrupted = set(mgrb.loc[mgrb["mgrb_disrupted"].astype(str).str.lower() == "true", "genome_id"])
    wide.loc[wide.index.isin(disrupted), "colistin"] = 1
    return wide


def score(y: np.ndarray, pred: np.ndarray, prob: np.ndarray | None) -> dict:
    tp = int(((y == 1) & (pred == 1)).sum())
    fn = int(((y == 1) & (pred == 0)).sum())
    fp = int(((y == 0) & (pred == 1)).sum())
    tn = int(((y == 0) & (pred == 0)).sum())
    sens = tp / (tp + fn) if tp + fn else math.nan
    spec = tn / (tn + fp) if tn + fp else math.nan
    s_lo, s_hi = wilson(tp, tp + fn)
    p_lo, p_hi = wilson(tn, tn + fp)
    auc = math.nan
    if prob is not None and len(set(y)) == 2:
        auc = roc_auc_score(y, prob)
    return {
        "n": len(y),
        "R": tp + fn,
        "S": tn + fp,
        "TP": tp,
        "FN": fn,
        "FP": fp,
        "TN": tn,
        "sensitivity": sens,
        "sens_lo": s_lo,
        "sens_hi": s_hi,
        "specificity": spec,
        "spec_lo": p_lo,
        "spec_hi": p_hi,
        "balanced_accuracy": (sens + spec) / 2,
        "auroc": auc,
    }


def run_drug(lab: pd.DataFrame, x: pd.DataFrame, rules: pd.DataFrame, drug: str) -> list[dict]:
    lab = lab.reset_index(drop=True)
    xd = x.loc[lab["genome_id"]].to_numpy(dtype=np.float32)
    y = lab["y"].to_numpy()
    rule = rules.reindex(lab["genome_id"])[drug].fillna(0).to_numpy().astype(int)
    rows, oof = [], []
    for scheme in SCHEMES:
        pred = {m: np.full(len(lab), -1) for m in ("rule", "logreg", "gboost")}
        prob = {m: np.full(len(lab), np.nan) for m in ("logreg", "gboost")}
        n_folds = 0
        for train, test in folds(lab, scheme):
            n_folds += 1
            pred["rule"][test] = rule[test]
            if len(set(y[train])) < 2:
                continue
            for m in ("logreg", "gboost"):
                model = make_model(m).fit(xd[train], y[train])
                p = model.predict_proba(xd[test])[:, 1]
                prob[m][test] = p
                pred[m][test] = (p >= 0.5).astype(int)
        extra = ""
        if scheme == "temporal":
            _, _, cutoff = temporal_split(lab)
            extra = f"train <= {cutoff}, test > {cutoff}"
        for m in ("rule", "logreg", "gboost"):
            done = pred[m] >= 0
            for view, keep in (
                ("all", done),
                ("no_borderline", done & ~lab["borderline"].to_numpy()),
            ):
                pr = prob[m][keep] if m != "rule" else None
                metrics = score(y[keep], pred[m][keep], pr)
                evaluable = min(metrics["R"], metrics["S"]) >= MIN_PER_CLASS
                note = extra
                if not evaluable:
                    note = f"{extra}; not evaluable ({metrics['R']} R / {metrics['S']} S in test)"
                    for k in ("sensitivity", "specificity", "balanced_accuracy", "auroc"):
                        metrics[k] = math.nan
                rows.append(
                    {
                        "drug": drug,
                        "scheme": scheme,
                        "model": m,
                        "labels": view,
                        "folds": n_folds,
                        "evaluable": evaluable,
                        "note": note.strip("; "),
                        **metrics,
                    }
                )
            oof.append(
                pd.DataFrame(
                    {
                        "genome_id": lab["genome_id"][done],
                        "drug": drug,
                        "scheme": scheme,
                        "model": m,
                        "y": y[done],
                        "pred": pred[m][done],
                        "prob": prob[m][done] if m != "rule" else np.nan,
                    }
                )
            )
    return rows, pd.concat(oof, ignore_index=True)


def fmt(v: float) -> str:
    return "n/a" if pd.isna(v) else f"{v:.3f}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, default=Path("data/features"))
    parser.add_argument("--amrfinder", type=Path, default=Path("results/cohort/amrfinder.parquet"))
    parser.add_argument("--mgrb", type=Path, default=Path("reports/review/mgrb_status.tsv"))
    parser.add_argument("--out", type=Path, default=Path("reports/benchmark"))
    parser.add_argument("--drugs", nargs="*", default=list(DRUGS))
    args = parser.parse_args()

    x = pd.read_parquet(args.features / "X.parquet")
    labels = pd.read_parquet(args.features / "labels.parquet")
    labels = labels[~labels["flagged"]]
    amr = pd.read_parquet(args.amrfinder)
    mgrb = pd.read_csv(args.mgrb, sep="\t", dtype={"genome_id": str, "mgrb_disrupted": str})
    rules = rule_calls(amr, mgrb, sorted(x.index))

    rows, oofs = [], []
    for drug in args.drugs:
        lab = labels[labels["drug"] == drug]
        print(f"{drug}: {len(lab):,} labels ({int(lab.y.sum()):,} R)", flush=True)
        r, o = run_drug(lab, x, rules, drug)
        rows += r
        oofs.append(o)
    res = pd.DataFrame(rows)
    args.out.mkdir(parents=True, exist_ok=True)
    res.to_csv(args.out / "results.csv", index=False, float_format="%.4f")
    pd.concat(oofs).to_parquet(args.features / "oof_predictions.parquet", index=False)

    main_view = res[res["labels"] == "all"]
    ba = main_view.pivot_table(
        index=["drug", "model"], columns="scheme", values="balanced_accuracy"
    ).reindex(columns=list(SCHEMES))
    auc = (
        main_view[main_view.model != "rule"]
        .pivot_table(index=["drug", "model"], columns="scheme", values="auroc")
        .reindex(columns=list(SCHEMES))
    )
    nb = (
        res[res["labels"] == "no_borderline"]
        .pivot_table(index=["drug", "model"], columns="scheme", values="balanced_accuracy")
        .reindex(columns=list(SCHEMES))
    )

    def table(df: pd.DataFrame) -> list[str]:
        out = [
            "| Drug | Model | " + " | ".join(df.columns) + " |",
            "|---|---|" + "---|" * len(df.columns),
        ]
        for (d, m), r in df.iterrows():
            out.append(f"| {d} | {m} | " + " | ".join(fmt(v) for v in r.tolist()) + " |")
        return out

    gap = (ba["random"] - ba[["study", "lineage"]].min(axis=1)).rename("gap")
    temporal_notes = main_view[
        (main_view.scheme == "temporal") & (main_view.model == "rule")
    ].drop_duplicates("drug")
    lines = [
        "# Leakage-aware benchmark",
        "",
        f"- Generated (UTC): {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        f"- Features: {x.shape[1]} (reports/benchmark/feature_summary.md); flagged labels excluded",
        "- Predictors: rule (AMRFinderPlus core; colistin + mgrB), logreg, gboost",
        "- Splits: random 5-fold; study = 5-fold grouped by BioProject; lineage = 5-fold grouped "
        "by ST; temporal = train on earlier, test on the latest ~20% of dated genomes",
        "- Temporal splits: "
        + "; ".join(f"{r.drug} {r.note}" for r in temporal_notes.itertuples()),
        f"- Splits whose test set has fewer than {MIN_PER_CLASS} R or {MIN_PER_CLASS} S genomes "
        "are not evaluable and shown as n/a.",
        "",
        "## Balanced accuracy, all test labels",
        "",
        *table(ba),
        "",
        "## AUROC, all test labels (models only)",
        "",
        *table(auc),
        "",
        "## Balanced accuracy without borderline-MIC labels (ADR-004)",
        "",
        *table(nb),
        "",
        "## Optimism of random splits",
        "",
        "Gap = balanced accuracy under random split minus the lower of study and lineage splits.",
        "",
        "| Drug | Model | Gap |",
        "|---|---|---|",
        *[f"| {d} | {m} | {fmt(v)} |" for (d, m), v in gap.items()],
        "",
        "Full metrics with sensitivity, specificity and 95% intervals: "
        "`reports/benchmark/results.csv`.",
    ]
    (args.out / "results.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
