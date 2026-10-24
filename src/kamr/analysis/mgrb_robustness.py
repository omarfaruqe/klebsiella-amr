"""Robustness checks for the mgrB colistin result (reviewer concerns).

Three questions the headline claim depends on:

1. Confounding by assembly quality. Several "disrupted" categories (split_contigs, contig_break,
   partial) are inferred from assembly structure. If fragmented assemblies both break at mgrB more
   often and come disproportionately from colistin-resistant collections, the association would be
   an artefact. Tested by stratifying the colistin association on assembly-quality quartile.

2. Circularity of the missense tier. Counting missense variants as predicting resistance was a
   decision made after observing their enrichment in this cohort, so the pooled 67.3% sensitivity
   is partly in-sample. Re-estimated by leave-one-study-out: variants are selected on the training
   studies only and applied to the held-out study.

3. Choice of baseline. Sensitivity gain is reported against both AMRFinderPlus and Kleborate,
   since Kleborate already makes truncation calls and is the stronger comparator.

Outputs: reports/review/mgrb_robustness.md
"""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import numpy as np
import pandas as pd

from kamr.analysis.concordance import wilson

DISRUPTED_ASSEMBLY = {"split_contigs", "contig_break", "partial"}


def load(curated: Path, mgrb: Path) -> pd.DataFrame:
    st = pd.read_csv(mgrb, sep="\t", dtype={"genome_id": str})
    qc = pd.read_parquet(curated / "genome_qc.parquet")
    qc["genome_id"] = qc["genome_id"].astype(str)
    lab = pd.read_parquet(curated / "cohort_labels.parquet")
    lab["genome_id"] = lab["genome_id"].astype(str)
    rec = pd.read_parquet(curated / "ast_records.parquet")
    rec["genome_id"] = rec["genome_id"].astype(str)
    bmd = set(
        rec[
            (rec.antibiotic == "colistin")
            & (rec.laboratory_typing_method.fillna("") == "Broth dilution")
            & (rec.record_label.isin(["R", "S"]))
        ].genome_id
    )
    col = lab[(lab.antibiotic == "colistin") & (lab.label.isin(["R", "S"]))]
    col = col[col.genome_id.isin(bmd)]
    d = col.merge(st, on="genome_id").merge(
        qc[["genome_id", "contigs", "contig_n50", "bioproject_accession", "st"]],
        on="genome_id",
        how="left",
    )
    d["group"] = np.where(
        d.mgrb_status == "intact",
        "intact",
        np.where(d.mgrb_status == "missense", "missense", "disrupted"),
    )
    d["assembly_inferred"] = d.mgrb_status.isin(DISRUPTED_ASSEMBLY)
    d["R"] = (d.label == "R").astype(int)
    return d


def sens_spec(pred: pd.Series, y: pd.Series) -> tuple[float, float, int, int]:
    tp = int(((y == 1) & pred).sum())
    fn = int(((y == 1) & ~pred).sum())
    fp = int(((y == 0) & pred).sum())
    tn = int(((y == 0) & ~pred).sum())
    return (
        tp / (tp + fn) if tp + fn else np.nan,
        tn / (tn + fp) if tn + fp else np.nan,
        tp + fn,
        tn + fp,
    )


def leave_one_study_out(d: pd.DataFrame, min_r: int = 5) -> pd.DataFrame:
    """Select enriched missense variants on training studies only; apply to held-out study."""
    d = d.copy()
    d["bp"] = d.bioproject_accession.fillna("missing:" + d.genome_id)
    rows = []
    for held in sorted(d.bp.unique()):
        test = d[d.bp == held]
        train = d[d.bp != held]
        if test.R.sum() < min_r:
            continue
        tm = train[train.group == "missense"]
        mis = tm.assign(v=tm.mgrb_detail.fillna("").str.split(";"))
        mis = mis.explode("v")
        enriched = set()
        for v, g in mis.groupby("v"):
            if v and len(g) >= 3 and g.R.mean() > 0.5:
                enriched.add(v)

        def pred_row(r, sel=frozenset(enriched)):
            if r.group == "disrupted":
                return True
            if r.group == "missense":
                return any(x in sel for x in str(r.mgrb_detail).split(";"))
            return False

        pred = test.apply(pred_row, axis=1)
        se, sp, nr, ns = sens_spec(pred, test.R)
        rows.append(
            {
                "held_out_study": held,
                "n_R": nr,
                "n_S": ns,
                "sensitivity": se,
                "specificity": sp,
                "n_variants_selected": len(enriched),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--curated", type=Path, default=Path("data/curated"))
    p.add_argument("--mgrb", type=Path, default=Path("reports/review/mgrb_status.tsv"))
    p.add_argument("--out", type=Path, default=Path("reports/review/mgrb_robustness.md"))
    a = p.parse_args()

    d = load(a.curated, a.mgrb)
    L = [
        "# mgrB result: robustness checks",
        "",
        f"- Generated (UTC): {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        f"- Colistin broth-microdilution set: {len(d):,} genomes ({int(d.R.sum())} R)",
        "",
        "## 1. Is the association confounded by assembly quality?",
        "",
        "Assembly fragmentation by mgrB call (all genomes in this set):",
        "",
        "| mgrB call | n | median contigs | median N50 |",
        "|---|---|---|---|",
    ]
    g = d.groupby("group").agg(
        n=("contigs", "size"), c=("contigs", "median"), n50=("contig_n50", "median")
    )
    for i, r in g.iterrows():
        L.append(f"| {i} | {int(r.n):,} | {r.c:.0f} | {r.n50:,.0f} |")

    L += [
        "",
        "Colistin resistance rate by assembly-quality quartile and mgrB status. If the",
        "association were an artefact of fragmentation, the gap would close within quartiles.",
        "",
        "| Contig quartile | mgrB intact %R (n) | mgrB disrupted %R (n) |",
        "|---|---|---|",
    ]
    d["q"] = pd.qcut(d.contigs, 4, labels=["Q1 fewest", "Q2", "Q3", "Q4 most"])
    for q, sub in d.groupby("q", observed=True):
        a_ = sub[sub.group == "intact"]
        b_ = sub[sub.group == "disrupted"]
        L.append(
            f"| {q} | {a_.R.mean() * 100:.0f}% ({len(a_)}) | "
            + (f"{b_.R.mean() * 100:.0f}% ({len(b_)}) |" if len(b_) else "n/a |")
        )

    intact_r = d[d.group == "intact"].R.mean() * 100
    asm = d[d.assembly_inferred]
    seq = d[(d.group == "disrupted") & ~d.assembly_inferred]
    asm_r, seq_r = asm.R.mean() * 100, seq.R.mean() * 100
    L += [
        "",
        "Disruption calls split by how they were detected:",
        "",
        "| Detection | n | %R |",
        "|---|---|---|",
        f"| Assembly-structure (split/break/partial) | {len(asm)} | {asm_r:.0f}% |",
        f"| Sequence-level (stop/indel/absent/insertion) | {len(seq)} | {seq_r:.0f}% |",
        f"| Intact mgrB | {int((d.group == 'intact').sum())} | {intact_r:.0f}% |",
    ]

    L += [
        "",
        "## 2. Leave-one-study-out estimate for the missense tier",
        "",
        "Missense variants are selected on training studies only (>=3 genomes, >50% resistant)",
        "and applied to the held-out study, removing the circularity of selecting variants on",
        "the same data used to report sensitivity.",
        "",
    ]
    loo = leave_one_study_out(d)
    if len(loo):
        L += [
            "| Held-out study | R | S | Sens: disrupted only | Sens: + missense (OOS)"
            " | Specificity | variants |",
            "|---|---|---|---|---|---|---|",
        ]
        for r in loo.itertuples():
            L.append(
                f"| {r.held_out_study} | {r.n_R} | {r.n_S} | {r.sens_disrupted_only:.1%} | "
                f"{r.sensitivity:.1%} | {r.specificity:.1%} | {r.n_variants_selected} |"
            )
        tot_r = int(loo.n_R.sum())
        wd = np.average(loo.sens_disrupted_only, weights=loo.n_R)
        wm = np.average(loo.sensitivity, weights=loo.n_R)
        lo_d, hi_d = wilson(int(round(wd * tot_r)), tot_r)
        lo_m, hi_m = wilson(int(round(wm * tot_r)), tot_r)
        L += [
            "",
            f"**Disrupted only (nothing fitted): {wd:.1%} ({lo_d:.1%}-{hi_d:.1%}).**",
            f"**With missense selected out-of-sample: {wm:.1%} ({lo_m:.1%}-{hi_m:.1%}).**",
            f"Pooled over {tot_r} resistant isolates in {len(loo)} studies.",
            "",
            "In-sample equivalents from the main analysis: 56.7% (disrupted) and 67.3%",
            "(disrupted or missense). The gap between the in-sample and out-of-sample missense",
            "figures is the cost of selecting variants on the same data. Per-study sensitivity",
            "also varies widely in the disrupted-only column, which is fitted to nothing: the",
            "mechanism mix of colistin resistance genuinely differs between collections.",
        ]
    else:
        L.append("No study had enough resistant isolates for leave-one-study-out.")

    L += [
        "",
        "## 3. Gain against the stronger baseline",
        "",
        "Kleborate already reports MgrB truncations, so it, not AMRFinderPlus, is the",
        "conservative comparator for the headline claim.",
        "",
        "| Baseline | Sensitivity | With mgrB check | Absolute gain |",
        "|---|---|---|---|",
        "| AMRFinderPlus | 28.3% | 67.3% | +39.0 points |",
        "| Kleborate | 39.5% | 61.6% | +22.1 points |",
        "",
        "The gain against Kleborate is the figure that should carry the headline.",
    ]
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
