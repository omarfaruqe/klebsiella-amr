"""Collate per-genome AMRFinderPlus and Kleborate outputs into cohort tables and a run summary.

Outputs: results/pilot/amrfinder.parquet, results/pilot/kleborate.parquet (not tracked)
         reports/pilot/run_summary.md (tracked)
"""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import pandas as pd


def read_amrfinder(directory: Path) -> pd.DataFrame:
    frames = []
    for path in sorted(directory.glob("*.tsv")):
        df = pd.read_csv(path, sep="\t", dtype=str)
        df.insert(0, "genome_id", path.stem)
        frames.append(df.drop(columns=["Name"], errors="ignore"))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def kleborate_full_table(directory: Path) -> pd.DataFrame | None:
    """Kleborate writes a short summary table and a full table; return the full one."""
    for path in sorted(directory.glob("*.txt")):
        df = pd.read_csv(path, sep="\t", dtype=str)
        if any(c.endswith("mlst__ST") for c in df.columns):
            return df
    return None


def read_kleborate(directory: Path) -> pd.DataFrame:
    frames = []
    for sub in sorted(p for p in directory.iterdir() if p.is_dir() and not p.name.startswith(".")):
        if not (sub / ".done").exists():
            continue
        df = kleborate_full_table(sub)
        if df is not None:
            df.insert(0, "genome_id", sub.name)
            frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="snapshots/assemblies-pilot-<date>.tsv")
    parser.add_argument("--results", type=Path, default=Path("results/pilot"))
    parser.add_argument("--report", type=Path, default=Path("reports/pilot/run_summary.md"))
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t", dtype={"genome_id": str})
    amr = read_amrfinder(args.results / "amrfinder")
    kl = read_kleborate(args.results / "kleborate")
    amr.to_parquet(args.results / "amrfinder.parquet", index=False)
    kl.to_parquet(args.results / "kleborate.parquet", index=False)

    log = pd.read_csv(args.results / "run_log.tsv", sep="\t", dtype={"genome_id": str})
    timed = log.drop_duplicates("genome_id", keep="last")
    ok_ids = set(manifest.loc[manifest["status"] == "ok", "genome_id"])
    amr_ids = set(amr["genome_id"]) if len(amr) else set()
    kl_ids = set(kl["genome_id"]) if len(kl) else set()
    versions = (args.results / "versions.txt").read_text().strip().splitlines()

    def median(col: str) -> str:
        vals = timed.loc[timed[col] > 0, col]
        return f"{vals.median():.0f} s" if len(vals) else "n/a"

    lines = [
        "# Pilot run summary",
        "",
        f"- Generated (UTC): {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        f"- Assembly manifest: `{args.manifest}`",
        *[f"- {v}" for v in versions],
        "",
        "## Assemblies",
        "",
        "| Status | Genomes |",
        "|---|---|",
        *[f"| {k} | {v} |" for k, v in manifest["status"].value_counts().items()],
        "",
        "## Tool runs",
        "",
        "| Tool | Genomes with output | Missing (of downloaded) | Median time |",
        "|---|---|---|---|",
        f"| AMRFinderPlus | {len(amr_ids)} | {len(ok_ids - amr_ids)} | {median('amrfinder_s')} |",
        f"| Kleborate | {len(kl_ids)} | {len(ok_ids - kl_ids)} | {median('kleborate_s')} |",
        "",
        f"AMRFinderPlus rows: {len(amr):,}. Genomes complete in both tools: "
        f"{len(amr_ids & kl_ids)} of {len(ok_ids)}.",
    ]
    missing = sorted(ok_ids - (amr_ids & kl_ids))
    if missing:
        lines += ["", "Incomplete genomes: " + ", ".join(missing[:50])]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
