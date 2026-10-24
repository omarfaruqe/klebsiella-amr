"""Pilot genotype-phenotype concordance: rule predictions vs laboratory labels, per drug and tool.

Inputs:  results/pilot/amrfinder.parquet, results/pilot/kleborate.parquet,
         data/curated/cohort_labels.parquet, data/curated/ast_records.parquet,
         reports/cohort/pilot_sample.tsv
Outputs: reports/pilot/concordance.csv, reports/pilot/concordance.md,
         reports/pilot/discordant.tsv (tracked)

Colistin follows ADR-003: the primary set uses broth-microdilution labels only; "colistin_all"
uses every colistin label for comparison.
"""

from __future__ import annotations

import argparse
import datetime as dt
import math
from pathlib import Path

import pandas as pd

from kamr.analysis.rules import (
    DRUGS,
    amrfinder_markers,
    kleborate_markers,
    predictions,
)


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (math.nan, math.nan)
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (max(0.0, centre - half), min(1.0, centre + half))


def metrics(df: pd.DataFrame) -> dict:
    """Confusion-matrix metrics for rows with columns label and predicted (R/S)."""
    tp = int(((df.label == "R") & (df.predicted == "R")).sum())
    fn = int(((df.label == "R") & (df.predicted == "S")).sum())
    fp = int(((df.label == "S") & (df.predicted == "R")).sum())
    tn = int(((df.label == "S") & (df.predicted == "S")).sum())
    sens_lo, sens_hi = wilson(tp, tp + fn)
    spec_lo, spec_hi = wilson(tn, tn + fp)
    n = tp + fn + fp + tn
    return {
        "n": n,
        "R": tp + fn,
        "S": tn + fp,
        "TP": tp,
        "FN": fn,
        "FP": fp,
        "TN": tn,
        "sensitivity": tp / (tp + fn) if tp + fn else math.nan,
        "sens_lo": sens_lo,
        "sens_hi": sens_hi,
        "specificity": tn / (tn + fp) if tn + fp else math.nan,
        "spec_lo": spec_lo,
        "spec_hi": spec_hi,
        "accuracy": (tp + tn) / n if n else math.nan,
    }


def evaluation_labels(
    labels: pd.DataFrame, records: pd.DataFrame, genomes: set[str]
) -> pd.DataFrame:
    """S/R labels for evaluated genomes; adds the ADR-003 broth-microdilution colistin set."""
    lab = labels[labels["label"].isin(["R", "S"]) & labels["genome_id"].isin(genomes)]
    lab = lab.rename(columns={"antibiotic": "drug"})[["genome_id", "drug", "label"]]
    col = records[
        (records["antibiotic"] == "colistin")
        & (records["laboratory_typing_method"].fillna("") == "Broth dilution")
        & records["record_label"].isin(["R", "S"])
    ]
    bmd = set(col["genome_id"])
    colistin = lab["drug"] == "colistin"
    extra = lab[colistin].assign(drug="colistin_all")
    primary = lab[~colistin | lab["genome_id"].isin(bmd)]
    return pd.concat([primary, extra], ignore_index=True)


def fmt(x: float) -> str:
    return "n/a" if pd.isna(x) else f"{x:.1%}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path("results/pilot"))
    parser.add_argument("--curated", type=Path, default=Path("data/curated"))
    parser.add_argument("--sample", type=Path, default=Path("reports/cohort/pilot_sample.tsv"))
    parser.add_argument("--reports", type=Path, default=Path("reports/pilot"))
    args = parser.parse_args()

    amr = pd.read_parquet(args.results / "amrfinder.parquet")
    kl = pd.read_parquet(args.results / "kleborate.parquet")
    labels = pd.read_parquet(args.curated / "cohort_labels.parquet")
    records = pd.read_parquet(args.curated / "ast_records.parquet")
    sample = pd.read_csv(args.sample, sep="\t", dtype=str)

    complete = set(amr["genome_id"]) & set(kl["genome_id"]) & set(sample["genome_id"])
    excluded = sorted(set(sample["genome_id"]) - complete)
    genomes = sorted(complete)

    preds = pd.concat(
        [
            predictions(
                amrfinder_markers(amr[amr.genome_id.isin(complete)]), genomes, "AMRFinderPlus"
            ),
            predictions(kleborate_markers(kl[kl.genome_id.isin(complete)]), genomes, "Kleborate"),
        ],
        ignore_index=True,
    )
    preds_all = pd.concat(
        [preds, preds[preds.drug == "colistin"].assign(drug="colistin_all")], ignore_index=True
    )
    evald = evaluation_labels(labels, records, complete).merge(
        preds_all, on=["genome_id", "drug"], how="inner"
    )

    order = [*DRUGS[:3], "colistin", "colistin_all", "tigecycline"]
    rows = []
    for drug in order:
        for tool in ("AMRFinderPlus", "Kleborate"):
            sub = evald[(evald.drug == drug) & (evald.tool == tool)]
            rows.append({"drug": drug, "tool": tool, **metrics(sub)})
    table = pd.DataFrame(rows)

    wide = evald.pivot_table(
        index=["genome_id", "drug", "label"],
        columns="tool",
        values=["predicted", "evidence"],
        aggfunc="first",
    )
    wide.columns = [f"{t}_{v}" for v, t in wide.columns]
    wide = wide.reset_index()
    agree = (
        wide.assign(same=wide["AMRFinderPlus_predicted"] == wide["Kleborate_predicted"])
        .groupby("drug")["same"]
        .mean()
    )
    st = sample.set_index("genome_id")["st"]
    wide["st"] = wide["genome_id"].map(st)
    disc = wide[
        (wide["AMRFinderPlus_predicted"] != wide["label"])
        | (wide["Kleborate_predicted"] != wide["label"])
    ]
    disc = disc[disc.drug != "colistin_all"].sort_values(["drug", "label", "genome_id"])

    args.reports.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.reports / "concordance.csv", index=False, float_format="%.4f")
    disc.to_csv(args.reports / "discordant.tsv", sep="\t", index=False)

    lines = [
        "# Pilot genotype-phenotype concordance",
        "",
        f"- Generated (UTC): {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        f"- Genomes evaluated: {len(genomes)} of {len(sample)} in the pilot sample"
        + (f" (excluded, incomplete tool output: {', '.join(excluded)})" if excluded else ""),
        "- Labels: EUCAST v16.1 re-interpreted laboratory S/R (I and conflicts excluded)",
        "- Colistin: broth microdilution only (ADR-003); colistin_all uses every method",
        "- Rules: see src/kamr/analysis/rules.py",
        "",
        "Sensitivity = resistant isolates the rule calls R; specificity = susceptible isolates the",
        "rule calls S. 95% Wilson intervals in brackets. The pilot over-samples rare STs, so these",
        "are tool-validation figures, not population estimates.",
        "",
        "| Drug | Tool | n (R/S) | Sensitivity | Specificity | FN | FP |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in table.itertuples():
        lines.append(
            f"| {r.drug} | {r.tool} | {r.n} ({r.R}/{r.S}) | {fmt(r.sensitivity)} "
            f"[{fmt(r.sens_lo)}-{fmt(r.sens_hi)}] | {fmt(r.specificity)} "
            f"[{fmt(r.spec_lo)}-{fmt(r.spec_hi)}] | {r.FN} | {r.FP} |"
        )
    lines += [
        "",
        "## AMRFinderPlus vs Kleborate prediction agreement",
        "",
        "| Drug | Agreement |",
        "|---|---|",
        *[f"| {d} | {fmt(agree.get(d, math.nan))} |" for d in order],
        "",
        f"Discordant genome-drug pairs (either tool wrong): {len(disc)}, listed in "
        "`reports/pilot/discordant.tsv` with each tool's evidence for review.",
    ]
    (args.reports / "concordance.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
