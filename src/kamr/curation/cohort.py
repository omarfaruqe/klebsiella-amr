"""Build the labelled cohort: laboratory AST for target drugs, re-interpreted under one standard.

Outputs (under data/curated/, not tracked):
  ast_records.parquet   one row per laboratory AST record with reported and re-interpreted calls
  cohort_labels.parquet one row per genome x drug with the final label
Reports (under reports/cohort/, tracked):
  label_summary.csv     per-drug label counts and recovery / agreement statistics
  label_report.md       human-readable data-quality report
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path

import pandas as pd

from kamr.curation.breakpoints import (
    MIC_UNITS,
    Breakpoint,
    interpret_mic,
    load_breakpoints,
    normalize_unit,
)

TARGET_DRUGS = ("meropenem", "ciprofloxacin", "gentamicin", "colistin", "tigecycline")
LAB = "laboratory method"
PHENO_MAP = {
    "resistant": "R",
    "susceptible": "S",
    "intermediate": "I",
    "susceptible-dose dependent": "I",
    "nonsusceptible": "R",
}
KEEP = (
    "id",
    "genome_id",
    "antibiotic",
    "resistant_phenotype",
    "measurement_sign",
    "measurement_value",
    "measurement_unit",
    "laboratory_typing_method",
    "testing_standard",
    "testing_standard_year",
)


def read_lab_records(path: Path, drugs: tuple[str, ...] = TARGET_DRUGS) -> pd.DataFrame:
    """Stream the JSONL snapshot, keeping only laboratory records for the target drugs."""
    wanted = set(drugs)
    rows = []
    with path.open() as fh:
        for line in fh:
            if "aboratory" not in line:
                continue
            rec = json.loads(line)
            if str(rec.get("evidence") or "").strip().lower() != LAB:
                continue
            drug = str(rec.get("antibiotic") or "").strip().lower()
            if drug not in wanted:
                continue
            row = {k: rec.get(k) for k in KEEP}
            row["antibiotic"] = drug
            rows.append(row)
    df = pd.DataFrame(rows, columns=list(KEEP))
    df["genome_id"] = df["genome_id"].astype(str)
    return df


def annotate_records(df: pd.DataFrame, bps: dict[str, Breakpoint]) -> pd.DataFrame:
    """Add reported call, re-interpreted call and the record-level label used downstream."""
    out = df.copy()
    pheno = out["resistant_phenotype"].fillna("").str.strip().str.lower()
    out["reported_call"] = pheno.map(PHENO_MAP)
    out["unit_norm"] = out["measurement_unit"].map(normalize_unit)
    out["has_mic"] = (
        out["unit_norm"].isin(MIC_UNITS)
        & pd.to_numeric(out["measurement_value"], errors="coerce").notna()
    )

    def reinterpret(row: pd.Series) -> str | None:
        bp = bps.get(row["antibiotic"])
        if bp is None or not bp.applies or not row["has_mic"]:
            return None
        return interpret_mic(row["measurement_sign"], row["measurement_value"], bp)

    out["eucast_call"] = out.apply(reinterpret, axis=1)
    out["record_label"] = out["eucast_call"].where(out["eucast_call"].notna(), out["reported_call"])
    out["label_source"] = "none"
    out.loc[out["reported_call"].notna(), "label_source"] = "reported"
    out.loc[out["eucast_call"].notna(), "label_source"] = "eucast_reinterpreted"
    return out


def aggregate_labels(records: pd.DataFrame) -> pd.DataFrame:
    """One label per genome x drug. S and R together is a conflict; I only if nothing else."""

    def resolve(group: pd.DataFrame) -> pd.Series:
        labels = set(group["record_label"].dropna())
        decisive = labels & {"S", "R"}
        if len(decisive) == 2:
            label = "conflict"
        elif decisive:
            label = decisive.pop()
        elif "I" in labels:
            label = "I"
        else:
            label = "unlabelled"
        sources = sorted(set(group.loc[group["record_label"].notna(), "label_source"]))
        return pd.Series(
            {
                "label": label,
                "n_records": len(group),
                "label_sources": ";".join(sources),
                "recovered_from_unknown": bool(
                    group["reported_call"].isna().all() and group["eucast_call"].notna().any()
                ),
            }
        )

    keys = ["genome_id", "antibiotic"]
    cols = ["record_label", "label_source", "reported_call", "eucast_call"]
    return records.groupby(keys)[cols].apply(resolve).reset_index()


def summarize(records: pd.DataFrame, labels: pd.DataFrame) -> pd.DataFrame:
    counts = labels.pivot_table(
        index="antibiotic", columns="label", values="genome_id", aggfunc="count", fill_value=0
    )
    for col in ("S", "R", "I", "conflict", "unlabelled"):
        if col not in counts:
            counts[col] = 0
    counts = counts[["R", "S", "I", "conflict", "unlabelled"]]
    counts["genomes"] = labels.groupby("antibiotic")["genome_id"].nunique()
    counts["recovered_from_unknown"] = labels.groupby("antibiotic")["recovered_from_unknown"].sum()

    both = records[records["reported_call"].isin(["S", "R"]) & records["eucast_call"].notna()]
    agree = (both["reported_call"] == both["eucast_call"]).groupby(both["antibiotic"])
    counts["records_both_calls"] = agree.size()
    counts["reported_vs_eucast_agreement"] = agree.mean().round(3)
    return counts.fillna(0).sort_values("genomes", ascending=False)


def write_report(
    summary: pd.DataFrame, snapshot: Path, bp_path: Path, verified: bool, n_records: int, out: Path
) -> None:
    now = dt.datetime.now(dt.UTC).isoformat(timespec="seconds")
    lines = [
        "# Cohort label report",
        "",
        f"- Generated (UTC): {now}",
        f"- Snapshot: `{snapshot}`",
        f"- Breakpoints: `{bp_path}`" + ("" if verified else " **(UNVERIFIED)**"),
        f"- Laboratory records for target drugs: {n_records:,}",
        "",
        "Label policy (ADR-002): laboratory records only; MICs re-interpreted with one breakpoint",
        "table where it applies to K. pneumoniae, otherwise the reported call is used; genomes",
        "with both S and R for a drug are conflicts and excluded from modelling; I is kept",
        "separately.",
        "",
        "| Drug | Genomes | R | S | I | Conflict | Unlabelled | Recovered from Unknown |"
        " Reported vs EUCAST agreement |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for drug, r in summary.iterrows():
        agreement = r["reported_vs_eucast_agreement"]
        agreement_text = (
            f"{agreement:.1%} (n={int(r['records_both_calls']):,})" if agreement else "n/a"
        )
        lines.append(
            f"| {drug} | {int(r['genomes']):,} | {int(r['R']):,} | {int(r['S']):,} | "
            f"{int(r['I']):,} | {int(r['conflict']):,} | {int(r['unlabelled']):,} | "
            f"{int(r['recovered_from_unknown']):,} | {agreement_text} |"
        )
    out.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument(
        "--breakpoints", type=Path, default=Path("config/breakpoints_eucast_v16.1.csv")
    )
    parser.add_argument("--curated", type=Path, default=Path("data/curated"))
    parser.add_argument("--reports", type=Path, default=Path("reports/cohort"))
    parser.add_argument(
        "--allow-unverified",
        action="store_true",
        help="run even if the breakpoint table has rows not marked verified",
    )
    args = parser.parse_args()

    bps, verified = load_breakpoints(args.breakpoints)
    if not verified and not args.allow_unverified:
        raise SystemExit(
            f"Breakpoint table {args.breakpoints} has unverified rows. Check the values against "
            "the official EUCAST table, set verified=yes, or pass --allow-unverified."
        )

    print(f"reading {args.snapshot} ...")
    raw = read_lab_records(args.snapshot)
    print(f"  laboratory records for target drugs: {len(raw):,}")
    records = annotate_records(raw, bps)
    labels = aggregate_labels(records)
    summary = summarize(records, labels)

    args.curated.mkdir(parents=True, exist_ok=True)
    args.reports.mkdir(parents=True, exist_ok=True)
    records.to_parquet(args.curated / "ast_records.parquet", index=False)
    labels.to_parquet(args.curated / "cohort_labels.parquet", index=False)
    summary.to_csv(args.reports / "label_summary.csv")
    write_report(
        summary,
        args.snapshot,
        args.breakpoints,
        verified,
        len(raw),
        args.reports / "label_report.md",
    )
    print(summary.to_string())
    print(f"\nwrote {args.curated}/ast_records.parquet, {args.curated}/cohort_labels.parquet")
    print(f"wrote {args.reports}/label_summary.csv, {args.reports}/label_report.md")


if __name__ == "__main__":
    main()
