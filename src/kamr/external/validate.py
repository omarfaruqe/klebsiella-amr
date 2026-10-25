"""Apply the colistin rules, frozen, to the held-out external cohort.

This module deliberately exposes no tunable parameter. It imports `colistin_rules` from the
analysis package rather than reimplementing it, and imports `metrics` for the same reason, so the
rule evaluated here is the same object evaluated internally and cannot drift from it. There is no
threshold to adjust, no variant list to extend and no selection step: whatever the rule does on
these genomes is what it does. That is the only thing separating external validation from a
second round of tuning.

The external labels carry a provenance column, so results are reported overall and restricted to
isolates whose label was re-derived from an MIC under our pinned EUCAST table. If the two strata
disagree, the headline should be the MIC-derived one, because the other inherits each submitter's
own breakpoint choice.

Inputs:
  data/external/ncbi/colistin_external.csv     labels, from kamr.external.ncbi_ast
  data/external/ncbi/assemblies/<acc>.fna      genomes, from kamr.external.fetch_assemblies
  results/external/{amrfinder,kleborate}.parquet   from kamr.pipeline.collate

Outputs:
  reports/external/mgrb_status.tsv
  reports/external/colistin_external.md
  reports/external/colistin_external.csv
"""

from __future__ import annotations

import argparse
import datetime as dt
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from kamr.analysis.concordance import fmt, metrics
from kamr.analysis.mgrb import (
    DISRUPTED,
    MIN_HIT_AA,
    QUERY,
    classify,
    colistin_rules,
    read_hits,
    run_tblastn,
)


def mgrb_status(assemblies: Path, genomes: list[str], hit_dir: Path, jobs: int) -> pd.DataFrame:
    """tblastn MgrB against each external assembly and classify the result."""
    hit_dir.mkdir(parents=True, exist_ok=True)
    failed = []
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futs = {
            pool.submit(run_tblastn, assemblies / f"{g}.fna", hit_dir / f"{g}.tsv"): g
            for g in genomes
        }
        for n, fut in enumerate(as_completed(futs), start=1):
            if not fut.result():
                failed.append(futs[fut])
            if n % 50 == 0 or n == len(genomes):
                print(f"  tblastn {n}/{len(genomes)}, {len(failed)} failed", flush=True)
    rows = []
    for g in genomes:
        if g in failed:
            continue
        state, detail = classify(read_hits(hit_dir / f"{g}.tsv"))
        rows.append({"genome_id": g, "mgrb_status": state, "mgrb_detail": detail})
    status = pd.DataFrame(rows)
    status["mgrb_disrupted"] = status["mgrb_status"].isin(DISRUPTED)
    return status


def evaluate(labels: pd.DataFrame, tables: dict[str, pd.Series]) -> pd.DataFrame:
    """Sensitivity and specificity with Wilson intervals, for each rule, on these labels."""
    rows = []
    for name, pred in tables.items():
        df = labels[["genome_id", "label"]].assign(predicted=labels["genome_id"].map(pred))
        m = metrics(df.dropna())
        rows.append({"rule": name, "n": len(df.dropna()), **m})
    return pd.DataFrame(rows)


def table(res: pd.DataFrame) -> list[str]:
    lines = ["| Rule | n | Sensitivity | Specificity | FN | FP |", "|---|---|---|---|---|---|"]
    for r in res.itertuples(index=False):
        lines.append(
            f"| {r.rule} | {r.n:,} "
            f"| {fmt(r.sensitivity)} [{fmt(r.sens_lo)}-{fmt(r.sens_hi)}] "
            f"| {fmt(r.specificity)} [{fmt(r.spec_lo)}-{fmt(r.spec_hi)}] "
            f"| {r.FN} | {r.FP} |"
        )
    return lines


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--labels", type=Path, default=Path("data/external/ncbi/colistin_external.csv"))
    p.add_argument("--assemblies", type=Path, default=Path("data/external/ncbi/assemblies"))
    p.add_argument("--results", type=Path, default=Path("results/external"))
    p.add_argument("--out", type=Path, default=Path("reports/external"))
    p.add_argument("--jobs", type=int, default=16)
    a = p.parse_args()

    lab = pd.read_csv(a.labels, dtype=str)
    lab = lab.rename(columns={"assembly": "genome_id"})
    lab = lab[lab["label"].isin(["R", "S"])]  # the single I isolate cannot inform sens/spec
    amr = pd.read_parquet(a.results / "amrfinder.parquet")
    kl = pd.read_parquet(a.results / "kleborate.parquet")

    with_tools = set(amr["genome_id"]) & set(kl["genome_id"])
    missing = sorted(set(lab["genome_id"]) - with_tools)
    lab = lab[lab["genome_id"].isin(with_tools)]
    genomes = sorted(lab["genome_id"])
    print(f"{len(genomes):,} external genomes with labels and tool output")
    if missing:
        print(f"{len(missing):,} labelled genomes lack tool output and are excluded")

    status = mgrb_status(a.assemblies, genomes, a.results / "mgrb_hits", a.jobs)
    tables = colistin_rules(amr, kl, status, genomes)

    overall = evaluate(lab, tables)
    mic_only = lab[lab["label_source"] == "mic_eucast"]
    mic = evaluate(mic_only, tables)

    a.out.mkdir(parents=True, exist_ok=True)
    status.to_csv(a.out / "mgrb_status.tsv", sep="\t", index=False)
    overall.assign(stratum="all").to_csv(a.out / "colistin_external.csv", index=False)

    counts = lab["label"].value_counts().to_dict()
    lines = [
        "# External validation of the colistin rules (NCBI Pathogen Detection)",
        "",
        f"- Generated (UTC): {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        f"- Genomes: {len(genomes):,} ({counts.get('R', 0):,} R, {counts.get('S', 0):,} S)",
        f"- Query: MgrB {QUERY.name}, tblastn e-value <= 1e-5, hits >= {MIN_HIT_AA} aa",
        "- Rules imported unchanged from `kamr.analysis.mgrb.colistin_rules`; nothing was",
        "  refitted, reselected or thresholded for this cohort.",
        "",
        "These isolates are not in the BV-BRC cohort, matched on BioSample accession and on",
        "assembly accession ignoring the version suffix. Their phenotypes were submitted to NCBI",
        "as BioSample antibiograms and were not available to us when the rule was specified.",
        "",
        "## mgrB status",
        "",
        "| Status | Genomes |",
        "|---|---|",
        *[f"| {k} | {v:,} |" for k, v in status["mgrb_status"].value_counts().items()],
        "",
        "## All external isolates",
        "",
        *table(overall),
        "",
        f"## Restricted to MIC-derived labels (n = {len(mic_only):,})",
        "",
        "Labels re-derived from the submitted MIC under our pinned EUCAST table, excluding",
        "isolates that retain the submitter's own interpretation.",
        "",
        *table(mic),
    ]
    text = "\n".join(lines) + "\n"
    (a.out / "colistin_external.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
