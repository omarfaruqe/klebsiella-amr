"""Build the benchmark inputs: a genome x feature matrix and a label table with split metadata.

Features (binary presence per genome):
  - AMRFinderPlus AMR elements, intrinsic oqxA/oqxB removed. Acquired genes are collapsed to gene
    family (blaKPC-2 -> blaKPC, rmtB1 -> rmtB); point mutations keep the full call (gyrA_S83I).
  - Efflux-regulator variants are kept: a learned model may use them where the core rule does not.
  - mgrB status from the Sprint 2 check: mgrB_disrupted, mgrB_missense.
  Features present in fewer than --min-count genomes are dropped.

Labels (one row per genome x drug, S/R only, ADR-002/003/004):
  - colistin uses broth-microdilution labels only;
  - highest MIC per genome and drug, and a borderline flag when that MIC is exactly the highest
    susceptible value (s_le) or the lowest resistant doubling dilution (2 x r_gt);
  - a flag from config/label_flags.csv (study-level label problems found in review);
  - split metadata: BioProject, ST, collection year.

Outputs: data/features/X.parquet, data/features/labels.parquet (not tracked),
         reports/benchmark/feature_summary.md (tracked)
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
from pathlib import Path

import pandas as pd

from kamr.analysis.concordance import evaluation_labels
from kamr.analysis.rules import INTRINSIC
from kamr.curation.breakpoints import load_breakpoints

DRUGS = ("meropenem", "ciprofloxacin", "gentamicin", "colistin", "tigecycline")
ALLELE_SUFFIX = re.compile(r"(-\d+[A-Za-z]?|\d+)$")


def gene_family(symbol: str) -> str:
    """blaKPC-2 -> blaKPC, rmtB1 -> rmtB, qnrS1 -> qnrS, aac(6')-Ib-cr5 -> aac(6')-Ib-cr."""
    if symbol.startswith("tet(") or symbol.endswith(")"):
        return symbol
    return ALLELE_SUFFIX.sub("", symbol) or symbol


def amr_features(amr: pd.DataFrame) -> pd.DataFrame:
    """Long table genome_id, feature from AMRFinderPlus AMR elements."""
    sym = amr["Element symbol"].fillna("").astype(str)
    etype = amr["Type"].fillna("").astype(str).str.upper()
    subtype = amr["Subtype"].fillna("").astype(str).str.upper()
    keep = (etype == "AMR") & ~sym.str.match(INTRINSIC) & (sym != "")
    point = subtype.str.startswith("POINT")
    feat = sym.where(point, sym.map(gene_family))
    return pd.DataFrame({"genome_id": amr.loc[keep, "genome_id"], "feature": feat[keep]})


def mgrb_features(status: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gid, st, disrupted in status[["genome_id", "mgrb_status", "mgrb_disrupted"]].itertuples(
        index=False
    ):
        if str(disrupted).lower() == "true":
            rows.append((gid, "mgrB_disrupted"))
        elif st == "missense":
            rows.append((gid, "mgrB_missense"))
    return pd.DataFrame(rows, columns=["genome_id", "feature"])


def feature_matrix(long: pd.DataFrame, genomes: list[str], min_count: int) -> pd.DataFrame:
    long = long.drop_duplicates()
    counts = long.groupby("feature")["genome_id"].nunique()
    keep = counts[counts >= min_count].index
    long = long[long["feature"].isin(keep)]
    x = pd.crosstab(long["genome_id"], long["feature"]).clip(upper=1).astype("uint8")
    return x.reindex(index=genomes, fill_value=0).sort_index(axis=1)


def highest_mics(records: pd.DataFrame) -> pd.DataFrame:
    rec = records[records["eucast_call"].notna()].copy()
    rec["mic"] = pd.to_numeric(rec["measurement_value"], errors="coerce")
    sign = rec["measurement_sign"].fillna("").astype(str).str.strip()
    rec["exact"] = sign.isin(["", "=", "=="])
    rec = rec.dropna(subset=["mic"]).sort_values("mic")
    top = rec.groupby(["genome_id", "antibiotic"]).tail(1)
    return top.rename(columns={"antibiotic": "drug"})[["genome_id", "drug", "mic", "exact"]]


def borderline(mic: float, exact: bool, s_le: float, r_gt: float) -> bool:
    if pd.isna(mic) or not exact:
        return False
    return mic == s_le or mic == 2 * r_gt


def label_table(labels, records, qc, flags, bps, genomes) -> pd.DataFrame:
    ev = evaluation_labels(labels, records, set(genomes))
    ev = ev[ev["drug"].isin(DRUGS)]
    lab = ev.merge(highest_mics(records), on=["genome_id", "drug"], how="left")
    lab["exact"] = lab["exact"].fillna(False).astype(bool)

    def is_borderline(r) -> bool:
        bp = bps.get(r.drug)
        return bool(bp and bp.applies and borderline(r.mic, r.exact, bp.s_le, bp.r_gt))

    lab["borderline"] = [is_borderline(r) for r in lab.itertuples()]
    meta = qc.set_index("genome_id")[["bioproject_accession", "st", "collection_year"]]
    lab = lab.join(meta, on="genome_id")
    excl = flags[flags["action"] == "exclude"][["bioproject_accession", "antibiotic"]]
    keys = set(zip(excl["bioproject_accession"], excl["antibiotic"], strict=True))
    lab["flagged"] = [
        (b, d) in keys for b, d in zip(lab["bioproject_accession"], lab["drug"], strict=True)
    ]
    lab["y"] = (lab["label"] == "R").astype("uint8")
    return lab.sort_values(["drug", "genome_id"]).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path("results/cohort"))
    parser.add_argument("--curated", type=Path, default=Path("data/curated"))
    parser.add_argument("--genomes", type=Path, default=Path("reports/cohort/qc_pass_genomes.tsv"))
    parser.add_argument("--mgrb", type=Path, default=Path("reports/review/mgrb_status.tsv"))
    parser.add_argument("--flags", type=Path, default=Path("config/label_flags.csv"))
    parser.add_argument(
        "--breakpoints", type=Path, default=Path("config/breakpoints_eucast_v16.1.csv")
    )
    parser.add_argument("--out", type=Path, default=Path("data/features"))
    parser.add_argument("--report", type=Path, default=Path("reports/benchmark/feature_summary.md"))
    parser.add_argument("--min-count", type=int, default=5)
    args = parser.parse_args()

    genomes = sorted(pd.read_csv(args.genomes, sep="\t", dtype=str)["genome_id"])
    amr = pd.read_parquet(args.results / "amrfinder.parquet")
    amr = amr[amr["genome_id"].isin(genomes)]
    status = pd.read_csv(args.mgrb, sep="\t", dtype={"genome_id": str, "mgrb_disrupted": str})
    long = pd.concat([amr_features(amr), mgrb_features(status)], ignore_index=True)
    x = feature_matrix(long, genomes, args.min_count)

    bps, _ = load_breakpoints(args.breakpoints)
    lab = label_table(
        pd.read_parquet(args.curated / "cohort_labels.parquet"),
        pd.read_parquet(args.curated / "ast_records.parquet"),
        pd.read_parquet(args.curated / "genome_qc.parquet"),
        pd.read_csv(args.flags, dtype=str),
        bps,
        genomes,
    )

    args.out.mkdir(parents=True, exist_ok=True)
    x.to_parquet(args.out / "X.parquet")
    lab.to_parquet(args.out / "labels.parquet", index=False)

    n_point = sum("_" in c and not c.startswith("mgrB_") for c in x.columns)
    prevalence = x.mean().sort_values(ascending=False)
    uniq = lab.drop_duplicates("genome_id")
    lines = [
        "# Benchmark inputs",
        "",
        f"- Generated (UTC): {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        f"- Genomes: {x.shape[0]:,}; features: {x.shape[1]:,} (present in >= {args.min_count} "
        f"genomes; {n_point} point-mutation features, 2 mgrB features, rest gene families)",
        "- Colistin labels: broth microdilution only (ADR-003)",
        "- Borderline: highest MIC exactly at s_le or at 2 x r_gt (ADR-004)",
        "- Flagged: study-level exclusions from config/label_flags.csv (ADR-004)",
        "",
        "## Labels per drug",
        "",
        "| Drug | Labelled | R | S | Borderline | Flagged | Primary set (not flagged) |",
        "|---|---|---|---|---|---|---|",
    ]
    for drug in DRUGS:
        d = lab[lab["drug"] == drug]
        lines.append(
            f"| {drug} | {len(d):,} | {int(d.y.sum()):,} | {int((1 - d.y).sum()):,} | "
            f"{int(d.borderline.sum()):,} | {int(d.flagged.sum()):,} | "
            f"{int((~d.flagged).sum()):,} |"
        )
    lines += [
        "",
        "## Split metadata coverage (labelled genomes)",
        "",
        f"- Labelled genomes: {len(uniq):,}",
        f"- BioProject missing: {int(uniq['bioproject_accession'].isna().sum()):,}",
        f"- ST missing: {int(uniq['st'].isna().sum()):,}",
        f"- Collection year missing: {int(uniq['collection_year'].isna().sum()):,}",
        "",
        "## 20 most common features",
        "",
        "| Feature | Genomes | Share |",
        "|---|---|---|",
        *[f"| {f} | {int(x[f].sum()):,} | {p:.1%} |" for f, p in prevalence.head(20).items()],
    ]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
