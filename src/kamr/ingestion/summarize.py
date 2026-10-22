"""Per-drug inventory of laboratory AST phenotypes from a BV-BRC snapshot."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

LAB_EVIDENCE = "laboratory method"


def load_snapshot(path: Path) -> pd.DataFrame:
    return pd.read_json(path, lines=True, dtype={"genome_id": str})


def drug_inventory(df: pd.DataFrame) -> pd.DataFrame:
    """Unique genomes per antibiotic x phenotype, laboratory evidence only."""
    lab = df.copy()
    if "evidence" not in lab:
        lab["evidence"] = ""
    lab["evidence"] = lab["evidence"].fillna("").str.strip().str.lower()
    lab = lab[lab["evidence"] == LAB_EVIDENCE].copy()
    lab["antibiotic"] = lab["antibiotic"].fillna("").str.strip().str.lower()
    lab["resistant_phenotype"] = (
        lab["resistant_phenotype"].fillna("Unknown").str.strip().str.title()
    )

    counts = (
        lab.groupby(["antibiotic", "resistant_phenotype"])["genome_id"]
        .nunique()
        .unstack(fill_value=0)
    )
    counts["genomes"] = lab.groupby("antibiotic")["genome_id"].nunique()
    per_genome = lab.groupby(["antibiotic", "genome_id"])["resistant_phenotype"].nunique()
    counts["conflicting"] = (per_genome > 1).groupby(level="antibiotic").sum()
    return counts.sort_values("genomes", ascending=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("--out", type=Path, default=Path("reports/spike/bvbrc_drug_inventory.csv"))
    parser.add_argument("--top", type=int, default=25)
    args = parser.parse_args()

    df = load_snapshot(args.snapshot)
    print(f"records: {len(df):,}   unique genomes: {df['genome_id'].nunique():,}")
    print("evidence types:")
    print(df.get("evidence", pd.Series(dtype=str)).value_counts(dropna=False).to_string())

    inv = drug_inventory(df)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    inv.to_csv(args.out)
    print(f"\nTop {args.top} drugs by genomes with laboratory AST:")
    print(inv.head(args.top).to_string())
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
