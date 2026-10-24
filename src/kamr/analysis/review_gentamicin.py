"""Sprint 2 review: why gentamicin rules miss resistant isolates.

The EUCAST v16.1 gentamicin breakpoint is a single value (S <= 2, R > 2 mg/L), so a one-dilution
shift moves an isolate across categories. Three hypotheses for rule false negatives:
  1. borderline MICs (R at 4 mg/L, one dilution above the breakpoint);
  2. aminoglycoside genes the rule does not count (enrichment in false negatives vs susceptible
     isolates without a gentamicin marker);
  3. study effects (false negatives concentrated in a few BioProjects or testing methods).
False positives (marker present, phenotype S) are summarised by marker and MIC.

Outputs: reports/review/gentamicin.md, reports/review/gentamicin_genomes.tsv
"""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import pandas as pd

from kamr.analysis.review_meropenem import md_table
from kamr.analysis.rules import INTRINSIC, amrfinder_markers, predictions

MIC_BINS = [0, 0.5, 1, 2, 4, 8, 16, 1e9]
MIC_LABELS = ["<=0.5", "1", "2", "4", "8", "16", ">16"]


def highest_mic(records: pd.DataFrame, drug: str) -> pd.DataFrame:
    """Highest re-interpretable MIC per genome for one drug, with its testing method."""
    rec = records[(records["antibiotic"] == drug) & records["eucast_call"].notna()].copy()
    rec["mic"] = pd.to_numeric(rec["measurement_value"], errors="coerce")
    rec = rec.dropna(subset=["mic"]).sort_values("mic")
    top = rec.groupby("genome_id").tail(1).set_index("genome_id")
    return pd.DataFrame(
        {
            "mic": top["mic"],
            "mic_bin": pd.cut(top["mic"], MIC_BINS, labels=MIC_LABELS).astype(str),
            "method": top["laboratory_typing_method"].fillna("not reported"),
        }
    )


def aminoglycoside_genes(amr: pd.DataFrame) -> pd.DataFrame:
    """Long table genome_id, gene for every acquired aminoglycoside element (Type AMR)."""
    sym = amr["Element symbol"].fillna("").astype(str)
    cls = amr["Class"].fillna("").astype(str).str.upper()
    etype = amr["Type"].fillna("").astype(str).str.upper()
    keep = (etype == "AMR") & cls.str.contains("AMINOGLYCOSIDE") & ~sym.str.match(INTRINSIC)
    return pd.DataFrame(
        {"genome_id": amr.loc[keep, "genome_id"], "gene": sym[keep]}
    ).drop_duplicates()


def enrichment(genes: pd.DataFrame, fn: set[str], tn: set[str], top: int = 15) -> pd.DataFrame:
    """Gene carriage in rule false negatives vs susceptible genomes without a gentamicin marker."""
    g = genes[genes["genome_id"].isin(fn | tn)]
    counts = g.assign(group=g["genome_id"].map(lambda x: "FN" if x in fn else "TN")).pivot_table(
        index="gene", columns="group", values="genome_id", aggfunc="nunique", fill_value=0
    )
    for col in ("FN", "TN"):
        if col not in counts:
            counts[col] = 0
    out = pd.DataFrame(
        {
            "in_FN": counts["FN"].astype(int),
            "in_TN": counts["TN"].astype(int),
            "FN_share": (counts["FN"] / max(len(fn), 1)).round(3),
            "TN_share": (counts["TN"] / max(len(tn), 1)).round(3),
        }
    )
    return out.sort_values(["in_FN", "in_TN"], ascending=[False, True]).head(top)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amrfinder", type=Path, default=Path("results/cohort/amrfinder.parquet"))
    parser.add_argument("--curated", type=Path, default=Path("data/curated"))
    parser.add_argument("--genomes", type=Path, default=Path("reports/cohort/qc_pass_genomes.tsv"))
    parser.add_argument("--out", type=Path, default=Path("reports/review"))
    args = parser.parse_args()

    genomes = set(pd.read_csv(args.genomes, sep="\t", dtype=str)["genome_id"])
    amr = pd.read_parquet(args.amrfinder)
    amr = amr[amr["genome_id"].isin(genomes)]
    labels = pd.read_parquet(args.curated / "cohort_labels.parquet")
    records = pd.read_parquet(args.curated / "ast_records.parquet")
    qc = pd.read_parquet(args.curated / "genome_qc.parquet").set_index("genome_id")

    lab = labels[(labels["antibiotic"] == "gentamicin") & labels["label"].isin(["R", "S"])]
    lab = lab[lab["genome_id"].isin(genomes)].set_index("genome_id")[["label"]]
    pred = predictions(amrfinder_markers(amr), sorted(lab.index), "AMRFinderPlus")
    pred = pred[pred["drug"] == "gentamicin"].set_index("genome_id")[["predicted", "evidence"]]
    t = (
        lab.join(pred)
        .join(highest_mic(records, "gentamicin"))
        .join(qc[["st", "bioproject_accession", "collection_year", "isolation_country"]])
    )
    t["outcome"] = t["label"] + "/" + t["predicted"]
    names = {"R/R": "TP", "R/S": "FN", "S/R": "FP", "S/S": "TN"}
    t["outcome"] = t["outcome"].map(names)
    fn, tn = set(t.index[t.outcome == "FN"]), set(t.index[t.outcome == "TN"])

    mic = pd.crosstab(t["outcome"], t["mic_bin"].fillna("no MIC"))
    mic = mic[[c for c in [*MIC_LABELS, "no MIC"] if c in mic.columns]]
    fn_mic = t.loc[list(fn)]
    fn_with = fn_mic[fn_mic["mic"].notna()]
    at4 = int((fn_with["mic_bin"] == "4").sum())

    genes = aminoglycoside_genes(amr)
    enr = enrichment(genes, fn, tn)
    fn_any_ag = genes[genes["genome_id"].isin(fn)]["genome_id"].nunique()

    r_all = t[t.label == "R"]
    bp = (
        pd.DataFrame(
            {
                "FN": t.loc[list(fn), "bioproject_accession"].fillna("missing").value_counts(),
                "all_R": r_all["bioproject_accession"].fillna("missing").value_counts(),
            }
        )
        .fillna(0)
        .astype(int)
    )
    bp["FN_share_of_project_R"] = (bp["FN"] / bp["all_R"].where(bp["all_R"] > 0)).round(3)
    bp = bp.sort_values("FN", ascending=False).head(8)
    method = pd.crosstab(t["method"].fillna("no MIC"), t["outcome"])
    method = method[[c for c in ("TP", "FN", "FP", "TN") if c in method.columns]]

    fp = t[t.outcome == "FP"]
    fp_markers = (
        fp["evidence"].str.split(";").explode().value_counts().head(10).to_frame("FP genomes")
    )

    lines = [
        "# Gentamicin discordance review",
        "",
        f"- Generated (UTC): {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        f"- Genomes with a gentamicin S/R label: {len(t):,} "
        f"(R {int((t.label == 'R').sum()):,}, S {int((t.label == 'S').sum()):,})",
        "- Rule: AMRFinderPlus core (Subclass GENTAMICIN or a 16S rRNA methyltransferase)",
        "- EUCAST v16.1 gentamicin (systemic): S <= 2, R > 2 mg/L",
        "",
        "| Outcome | Genomes |",
        "|---|---|",
        *[f"| {k} | {v:,} |" for k, v in t["outcome"].value_counts().items()],
        "",
        "## 1. MIC distribution by outcome (mg/L)",
        "",
        f"{len(fn_with):,} of {len(fn):,} false negatives have an MIC; {at4:,} of those "
        f"({at4 / max(len(fn_with), 1):.0%}) are at 4 mg/L, one dilution above the breakpoint.",
        "",
        *md_table(mic, "Outcome"),
        "",
        "## 2. Aminoglycoside genes: false negatives vs susceptible non-carriers",
        "",
        f"{fn_any_ag:,} of {len(fn):,} false negatives carry at least one acquired aminoglycoside "
        "gene that the gentamicin rule does not count. FN = rule false negatives (R, no "
        "gentamicin marker); TN = susceptible genomes without a gentamicin marker.",
        "",
        *md_table(enr, "Gene"),
        "",
        "## 3. Study and method effects",
        "",
        "BioProjects with the most false negatives:",
        "",
        *md_table(bp, "BioProject"),
        "",
        "Testing method of the highest MIC record:",
        "",
        *md_table(method, "Method"),
        "",
        "## 4. False positives (gentamicin marker present, phenotype S)",
        "",
        f"{len(fp):,} genomes. Most frequent markers:",
        "",
        *md_table(fp_markers, "Marker"),
        "",
        "Per-genome detail of false negatives and false positives: "
        "`reports/review/gentamicin_genomes.tsv`.",
    ]
    args.out.mkdir(parents=True, exist_ok=True)
    ag = genes.groupby("genome_id")["gene"].apply(lambda s: ";".join(sorted(s)))
    detail = t[t.outcome.isin(["FN", "FP"])].assign(aminoglycoside_genes=lambda d: d.index.map(ag))
    keep = [
        "label",
        "outcome",
        "evidence",
        "aminoglycoside_genes",
        "mic",
        "method",
        "st",
        "bioproject_accession",
        "collection_year",
        "isolation_country",
    ]
    detail[keep].sort_values(["outcome", "mic"]).to_csv(
        args.out / "gentamicin_genomes.tsv", sep="\t", index_label="genome_id"
    )
    (args.out / "gentamicin.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
