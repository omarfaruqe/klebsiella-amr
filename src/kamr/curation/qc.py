"""Assembly QC on labelled genomes, cohort data-quality report and a stratified pilot sample.

Inputs:  data/curated/cohort_labels.parquet, data/curated/ast_records.parquet,
         data/raw/bvbrc_genome/<date>/genome.jsonl, config/qc.toml
Outputs: data/curated/genome_qc.parquet (not tracked)
         reports/cohort/qc_report.md, reports/cohort/qc_funnel.csv,
         reports/cohort/pilot_sample.tsv (tracked)
"""

from __future__ import annotations

import argparse
import datetime as dt
import random
import re
import tomllib
from pathlib import Path

import pandas as pd

NUMERIC = (
    "contigs",
    "contig_n50",
    "genome_length",
    "gc_content",
    "checkm_completeness",
    "checkm_contamination",
    "collection_year",
)
ST_RE = re.compile(r"(?:^|[._\s])(?:ST)?(\d+)$", re.IGNORECASE)


def parse_st(value: object) -> str | None:
    """Extract the ST from BV-BRC mlst strings such as 'MLST.Klebsiella_pneumoniae.11'."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value[0] if isinstance(value, list) and value else value).strip()
    match = ST_RE.search(text)
    return f"ST{match.group(1)}" if match else None


def load_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_json(path, lines=True, dtype={"genome_id": str})
    for col in NUMERIC:
        meta[col] = pd.to_numeric(meta[col], errors="coerce") if col in meta else pd.NA
    for col in ("species", "genome_name", "assembly_accession", "bioproject_accession", "mlst"):
        if col not in meta:
            meta[col] = None
    meta["st"] = meta["mlst"].map(parse_st)
    return meta.drop_duplicates("genome_id").reset_index(drop=True)


def apply_qc(meta: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return per-genome QC results and a funnel table (filters applied in order)."""
    a = cfg["assembly"]
    pattern = a["species_pattern"]
    species_ok = (
        meta["species"].fillna("").astype(str).str.contains(pattern, case=False, regex=False)
    )
    # BV-BRC records can contradict themselves (e.g. species "Klebsiella pneumoniae" with
    # genome_name "Citrobacter freundii ..."), so the genome name must agree as well.
    name_ok = (
        meta["genome_name"].fillna("").astype(str).str.contains(pattern, case=False, regex=False)
    )
    comp = meta["checkm_completeness"]
    cont = meta["checkm_contamination"]
    steps = [
        ("species_field", species_ok),
        ("genome_name_species", name_ok),
        ("assembly_stats_present", meta[["genome_length", "contigs"]].notna().all(axis=1)),
        ("min_length", meta["genome_length"] >= a["min_length"]),
        ("max_length", meta["genome_length"] <= a["max_length"]),
        ("max_contigs", meta["contigs"] <= a["max_contigs"]),
        ("min_n50", meta["contig_n50"].isna() | (meta["contig_n50"] >= a["min_n50"])),
        ("checkm_completeness", comp.isna() | (comp >= a["min_checkm_completeness"])),
        ("checkm_contamination", cont.isna() | (cont <= a["max_checkm_contamination"])),
    ]
    out = meta.copy()
    out["qc_fail_reason"] = None
    remaining = pd.Series(True, index=out.index)
    funnel = [("labelled genomes with metadata", int(remaining.sum()), 0)]
    for step, ok in steps:
        fails = remaining & ~ok.fillna(False)
        out.loc[fails, "qc_fail_reason"] = step
        remaining &= ~fails
        funnel.append((step, int(remaining.sum()), int(fails.sum())))
    out["qc_pass"] = remaining
    out["checkm_missing"] = comp.isna()
    return out, pd.DataFrame(funnel, columns=["step", "remaining", "removed"])


def select_sample(
    qc: pd.DataFrame, labels: pd.DataFrame, drugs: list[str], per_stratum: int, seed: int
) -> pd.DataFrame:
    """Pick QC-passing genomes covering every drug x {R,S} stratum, spreading across STs.

    Within each stratum genomes are taken round-robin over STs (rarest ST first after a seeded
    shuffle), so a single dominant clone cannot fill the stratum. A genome chosen for one
    stratum counts towards every stratum it belongs to.
    """
    rng = random.Random(seed)
    passing = qc.loc[qc["qc_pass"], ["genome_id", "st"]].set_index("genome_id")["st"]
    lab = labels[labels["label"].isin(["R", "S"]) & labels["genome_id"].isin(passing.index)]
    chosen: dict[str, str] = {}

    for drug in drugs:
        for label in ("R", "S"):
            pool = sorted(
                lab.loc[(lab["antibiotic"] == drug) & (lab["label"] == label), "genome_id"]
            )
            have = sum(1 for g in pool if g in chosen)
            if have >= per_stratum:
                continue
            by_st: dict[str, list[str]] = {}
            for g in pool:
                if g not in chosen:
                    st = passing[g]
                    by_st.setdefault(st if isinstance(st, str) else "unknown", []).append(g)
            for genomes in by_st.values():
                rng.shuffle(genomes)
            order = sorted(by_st, key=lambda st: (len(by_st[st]), st))
            while have < per_stratum and any(by_st.values()):
                for st in order:
                    if by_st[st] and have < per_stratum:
                        chosen[by_st[st].pop()] = f"{drug}:{label}"
                        have += 1

    sample = pd.DataFrame(sorted(chosen.items()), columns=["genome_id", "selected_for"])
    sample["st"] = sample["genome_id"].map(passing)
    wide = lab.pivot_table(index="genome_id", columns="antibiotic", values="label", aggfunc="first")
    return sample.join(wide, on="genome_id")


def top_counts(series: pd.Series, n: int = 10) -> list[str]:
    counts = series.fillna("missing").value_counts()
    total = int(counts.sum())
    return [f"| {k} | {v:,} | {v / total:.1%} |" for k, v in counts.head(n).items()]


def write_report(
    qc: pd.DataFrame,
    funnel: pd.DataFrame,
    labels: pd.DataFrame,
    records: pd.DataFrame,
    sample: pd.DataFrame,
    drugs: list[str],
    meta_path: Path,
    out: Path,
) -> None:
    passing = set(qc.loc[qc["qc_pass"], "genome_id"])
    lab = labels[labels["label"].isin(["R", "S", "I"])]
    lines = [
        "# Cohort QC and data-quality report",
        "",
        f"- Generated (UTC): {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        f"- Genome metadata: `{meta_path}`",
        "- Thresholds: `config/qc.toml`",
        "",
        "## QC funnel",
        "",
        "| Step | Remaining | Removed |",
        "|---|---|---|",
        *[f"| {r.step} | {r.remaining:,} | {r.removed:,} |" for r in funnel.itertuples()],
        "",
        f"CheckM values missing (not failed): {int(qc['checkm_missing'].sum()):,} genomes.",
        f"Genomes without an assembly accession: {int(qc['assembly_accession'].isna().sum()):,}.",
        "",
        "## Modelling labels after QC",
        "",
        "| Drug | R | S | I | R share of R+S |",
        "|---|---|---|---|---|",
    ]
    for drug in drugs:
        d = lab[(lab["antibiotic"] == drug) & lab["genome_id"].isin(passing)][
            "label"
        ].value_counts()
        r, s, i = (int(d.get(k, 0)) for k in ("R", "S", "I"))
        share = f"{r / (r + s):.1%}" if r + s else "n/a"
        lines.append(f"| {drug} | {r:,} | {s:,} | {i:,} | {share} |")

    col = records[records["antibiotic"] == "colistin"].copy()
    col["method"] = col["laboratory_typing_method"].fillna("not reported")
    col = col[col["genome_id"].isin(passing)].merge(
        labels[labels["antibiotic"] == "colistin"][["genome_id", "label"]], on="genome_id"
    )
    method = col.drop_duplicates(["genome_id", "method"]).pivot_table(
        index="method", columns="label", values="genome_id", aggfunc="count", fill_value=0
    )
    lines += [
        "",
        "## Colistin labels by testing method (QC-passing genomes)",
        "",
        "EUCAST requires broth microdilution for colistin MICs.",
        "",
        "```",
        method.to_string(),
        "```",
    ]

    kept = qc[qc["qc_pass"]]
    years = kept["collection_year"].dropna().astype(int)
    bins = pd.cut(
        years,
        bins=[1900, 2004, 2009, 2014, 2019, 2100],
        labels=["<=2004", "2005-2009", "2010-2014", "2015-2019", ">=2020"],
    )
    lines += [
        "",
        "## Collection year (QC-passing genomes)",
        "",
        f"Missing year: {int(kept['collection_year'].isna().sum()):,} of {len(kept):,}.",
        "",
        "| Period | Genomes |",
        "|---|---|",
        *[f"| {k} | {v:,} |" for k, v in bins.value_counts().sort_index().items()],
    ]
    for title, column in (
        ("Sequence types", "st"),
        ("Isolation country", "isolation_country"),
        ("BioProject", "bioproject_accession"),
    ):
        lines += [
            "",
            f"## Top 10 {title.lower()}",
            "",
            f"| {title} | Genomes | Share |",
            "|---|---|---|",
            *top_counts(kept[column]),
        ]

    lines += [
        "",
        "## Pilot sample",
        "",
        f"{len(sample):,} genomes in `reports/cohort/pilot_sample.tsv`, "
        f"{sample['st'].nunique():,} distinct STs; "
        f"{int((~sample['genome_id'].isin(passing)).sum())} no longer pass QC.",
        "",
        "| Drug | R in sample | S in sample |",
        "|---|---|---|",
    ]
    for drug in drugs:
        vc = sample[drug].value_counts() if drug in sample else pd.Series(dtype=int)
        lines.append(f"| {drug} | {int(vc.get('R', 0))} | {int(vc.get('S', 0))} |")
    out.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", type=Path, help="data/raw/bvbrc_genome/<date>/genome.jsonl")
    parser.add_argument("--config", type=Path, default=Path("config/qc.toml"))
    parser.add_argument("--curated", type=Path, default=Path("data/curated"))
    parser.add_argument("--reports", type=Path, default=Path("reports/cohort"))
    parser.add_argument(
        "--resample", action="store_true", help="re-draw the pilot sample even if one exists"
    )
    args = parser.parse_args()

    cfg = tomllib.loads(args.config.read_text())
    drugs = cfg["sample"]["drugs"]
    labels = pd.read_parquet(args.curated / "cohort_labels.parquet")
    records = pd.read_parquet(args.curated / "ast_records.parquet")
    meta = load_metadata(args.metadata)

    qc, funnel = apply_qc(meta, cfg)
    sample_path = args.reports / "pilot_sample.tsv"
    frozen = sample_path.exists() and not args.resample
    if frozen:
        # The pilot sample is frozen once tools have run on it; re-draw only with --resample.
        sample = pd.read_csv(sample_path, sep="\t", dtype=str)
        print(f"using frozen pilot sample {sample_path} ({len(sample)} genomes)")
    else:
        per_stratum, seed = cfg["sample"]["per_stratum"], cfg["sample"]["seed"]
        sample = select_sample(qc, labels, drugs, per_stratum, seed)
        sample = sample.merge(qc[["genome_id", "assembly_accession"]], on="genome_id", how="left")

    qc.to_parquet(args.curated / "genome_qc.parquet", index=False)
    funnel.to_csv(args.reports / "qc_funnel.csv", index=False)
    if not frozen:
        sample.to_csv(sample_path, sep="\t", index=False)
    write_report(
        qc, funnel, labels, records, sample, drugs, args.metadata, args.reports / "qc_report.md"
    )
    print(funnel.to_string(index=False))
    print(f"\npilot sample: {len(sample):,} genomes, {sample['st'].nunique():,} STs")
    print(f"wrote {args.curated}/genome_qc.parquet and reports in {args.reports}/")


if __name__ == "__main__":
    main()
