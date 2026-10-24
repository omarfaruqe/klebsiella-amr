"""Targeted mgrB integrity check for every cohort genome, and its effect on colistin prediction.

mgrB encodes a small negative regulator of PhoPQ. Its loss (deletion, premature stop, or
insertion-sequence (IS) disruption) is the commonest route to colistin resistance in
K. pneumoniae. IS insertions often break short-read assemblies inside the gene, which gene and
point-mutation searches can miss. This module searches each assembly with tblastn using the
MgrB reference protein (config/reference/mgrB_WP_002911375.1.faa) and classifies the result.

Statuses: intact, missense, premature_stop, frameshift_or_indel, insertion_assembled,
split_contigs, contig_break, partial, absent. Everything except intact and missense counts as
"disrupted". Rules are reported with disrupted mgrB only, and with disrupted or missense mgrB,
so the contribution of missense variants is visible separately.

Missense calls compare the aligned query and subject residue by residue. Query position 1 is
ignored because tblastn translates alternative start codons (e.g. GTG) as V or L, not M.

Outputs: results/cohort/mgrb_hits/<genome_id>.tsv (raw tblastn hits, cached),
         reports/review/mgrb_status.tsv, reports/review/colistin_mgrb.md
"""

from __future__ import annotations

import argparse
import datetime as dt
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from kamr.analysis.concordance import evaluation_labels, fmt, metrics
from kamr.analysis.rules import amrfinder_markers, kleborate_markers, predictions

QUERY = Path("config/reference/mgrB_WP_002911375.1.faa")
QLEN = 47
COLUMNS = [
    "sseqid",
    "qstart",
    "qend",
    "sstart",
    "send",
    "pident",
    "length",
    "evalue",
    "sseq",
    "slen",
    "qseq",
]
EDGE_BP = 20
INSERTION_GAP_BP = 300
MIN_HIT_AA = 10
DISRUPTED = {
    "premature_stop",
    "frameshift_or_indel",
    "insertion_assembled",
    "split_contigs",
    "contig_break",
    "partial",
    "absent",
}


def run_tblastn(fasta: Path, out: Path, query: Path = QUERY) -> bool:
    if out.exists():
        return True
    tmp = out.with_suffix(".tsv.tmp")
    cmd = [
        "tblastn",
        "-query",
        str(query),
        "-subject",
        str(fasta),
        "-evalue",
        "1e-5",
        "-seg",
        "no",
        "-outfmt",
        "6 " + " ".join(COLUMNS),
        "-out",
        str(tmp),
    ]
    ok = subprocess.run(cmd, capture_output=True, check=False).returncode == 0
    if ok:
        tmp.replace(out)
    return ok


def read_hits(path: Path) -> pd.DataFrame:
    if path.stat().st_size == 0:
        return pd.DataFrame(columns=COLUMNS)
    return pd.read_csv(
        path, sep="\t", names=COLUMNS, dtype={"sseqid": str, "sseq": str, "qseq": str}
    )


def variants(qseq: str, sseq: str, qstart: int) -> list[str]:
    """Amino-acid differences (e.g. W47R) between aligned query and subject; ignores position 1."""
    out, qpos = [], qstart
    for qa, sa in zip(qseq, sseq, strict=False):
        if qa == "-":
            out.append(f"ins{qpos - 1}{sa}")
            continue
        if sa == "-":
            out.append(f"{qa}{qpos}del")
        elif qa != sa and qpos != 1:
            out.append(f"{qa}{qpos}{sa}")
        qpos += 1
    return out


def classify(hits: pd.DataFrame) -> tuple[str, str]:
    """Return (status, detail) for one genome's tblastn hits against MgrB."""
    hits = hits[hits["length"] >= MIN_HIT_AA].copy()
    if hits.empty:
        return "absent", ""
    hits["lo"] = hits[["sstart", "send"]].min(axis=1)
    hits["hi"] = hits[["sstart", "send"]].max(axis=1)
    hits["cover"] = hits["qend"] - hits["qstart"] + 1
    hits = hits.sort_values(["cover", "pident"], ascending=False)

    full = hits[(hits["qstart"] <= 2) & (hits["qend"] >= QLEN - 1)]
    if not full.empty:
        h = full.iloc[0]
        if "*" in str(h["sseq"]).rstrip("*"):
            return "premature_stop", f"stop in aligned protein at {h['sseqid']}"
        diffs = variants(str(h["qseq"]), str(h["sseq"]), int(h["qstart"]))
        if not diffs:
            return "intact", ""
        return "missense", ";".join(diffs)

    near_edge = (hits["lo"] <= EDGE_BP) | (hits["hi"] >= hits["slen"] - EDGE_BP)
    covered = set()
    for r in hits.itertuples():
        covered.update(range(int(r.qstart), int(r.qend) + 1))
    if len(hits) >= 2 and len(covered) >= QLEN - 7:
        if hits["sseqid"].nunique() == 1:
            ordered = hits.sort_values("lo")
            gap = int(ordered["lo"].iloc[1] - ordered["hi"].iloc[0])
            if gap > INSERTION_GAP_BP:
                return "insertion_assembled", f"{gap} bp between mgrB parts"
            return "frameshift_or_indel", f"{gap} bp between mgrB parts"
        return "split_contigs", ";".join(sorted(hits["sseqid"].unique()))
    if near_edge.any():
        h = hits[near_edge].iloc[0]
        return "contig_break", f"mgrB aa {int(h.qstart)}-{int(h.qend)} at edge of {h.sseqid}"
    h = hits.iloc[0]
    return "partial", f"mgrB aa {int(h.qstart)}-{int(h.qend)}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path, default=Path("snapshots/assemblies-cohort-2026-10-24.tsv")
    )
    parser.add_argument("--assemblies", type=Path, default=Path("data/assemblies"))
    parser.add_argument("--results", type=Path, default=Path("results/cohort"))
    parser.add_argument("--curated", type=Path, default=Path("data/curated"))
    parser.add_argument("--out", type=Path, default=Path("reports/review"))
    parser.add_argument("--jobs", type=int, default=16)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t", dtype={"genome_id": str})
    ids = manifest.loc[manifest["status"] == "ok", "genome_id"].tolist()
    hit_dir = args.results / "mgrb_hits"
    hit_dir.mkdir(parents=True, exist_ok=True)

    print(f"tblastn MgrB vs {len(ids)} assemblies ({args.jobs} at once) ...")
    failed = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = {
            pool.submit(run_tblastn, args.assemblies / f"{g}.fna", hit_dir / f"{g}.tsv"): g
            for g in ids
        }
        for n, fut in enumerate(as_completed(futs), start=1):
            if not fut.result():
                failed.append(futs[fut])
            if n % 1000 == 0 or n == len(ids):
                print(f"  {n}/{len(ids)} done, {len(failed)} failed", flush=True)

    rows = []
    for g in ids:
        if g in failed:
            continue
        status, detail = classify(read_hits(hit_dir / f"{g}.tsv"))
        rows.append({"genome_id": g, "mgrb_status": status, "mgrb_detail": detail})
    status = pd.DataFrame(rows)
    status["mgrb_disrupted"] = status["mgrb_status"].isin(DISRUPTED)

    genomes = set(status["genome_id"])
    amr = pd.read_parquet(args.results / "amrfinder.parquet")
    kl = pd.read_parquet(args.results / "kleborate.parquet")
    labels = pd.read_parquet(args.curated / "cohort_labels.parquet")
    records = pd.read_parquet(args.curated / "ast_records.parquet")
    ev = evaluation_labels(labels, records, genomes)
    ev = ev[ev["drug"].isin(["colistin", "colistin_all"])]

    gl = sorted(genomes)
    base = {
        "AMRFinderPlus": predictions(amrfinder_markers(amr[amr.genome_id.isin(genomes)]), gl, "a"),
        "Kleborate": predictions(kleborate_markers(kl[kl.genome_id.isin(genomes)]), gl, "k"),
    }
    disrupted = set(status.loc[status["mgrb_disrupted"], "genome_id"])
    altered = disrupted | set(status.loc[status["mgrb_status"] == "missense", "genome_id"])
    tables = {}
    for name, p in base.items():
        p = p[p["drug"] == "colistin"][["genome_id", "predicted"]].set_index("genome_id")
        tables[name] = p["predicted"]
        tables[f"{name} + mgrB disrupted"] = p["predicted"].where(~p.index.isin(disrupted), "R")
        tables[f"{name} + mgrB disrupted or missense"] = p["predicted"].where(
            ~p.index.isin(altered), "R"
        )
    tables["mgrB disrupted alone"] = pd.Series(
        ["R" if g in disrupted else "S" for g in gl], index=gl
    )
    tables["mgrB disrupted or missense alone"] = pd.Series(
        ["R" if g in altered else "S" for g in gl], index=gl
    )

    lines = [
        "# Colistin: targeted mgrB check",
        "",
        f"- Generated (UTC): {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        f"- Genomes searched: {len(status):,} (tblastn failures: {len(failed)})",
        f"- Query: MgrB {QUERY.name} (47 aa), tblastn e-value <= 1e-5, hits >= {MIN_HIT_AA} aa",
        "- Labels: colistin = broth microdilution only (ADR-003); colistin_all = every method",
        "",
        "## mgrB status across all searched genomes",
        "",
        "| Status | Genomes |",
        "|---|---|",
        *[f"| {k} | {v:,} |" for k, v in status["mgrb_status"].value_counts().items()],
    ]
    bmd = ev[ev["drug"] == "colistin"].merge(status, on="genome_id")
    mis = (
        bmd[bmd["mgrb_status"] == "missense"]
        .assign(variant=lambda d: d["mgrb_detail"].str.split(";"))
        .explode("variant")
    )
    var = pd.crosstab(mis["variant"], mis["label"])
    for col in ("R", "S"):
        if col not in var:
            var[col] = 0
    var = var.assign(total=var["R"] + var["S"]).sort_values("total", ascending=False).head(15)
    lines += [
        "",
        "## mgrB missense variants (colistin, broth microdilution)",
        "",
        "Missense variants are not counted as disrupted. Top 15 by genomes:",
        "",
        "| Variant | R | S |",
        "|---|---|---|",
        *[f"| {i} | {int(r['R']):,} | {int(r['S']):,} |" for i, r in var.iterrows()],
    ]
    for drug in ("colistin", "colistin_all"):
        lab = ev[ev["drug"] == drug].merge(status, on="genome_id")
        ct = pd.crosstab(lab["mgrb_status"], lab["label"])
        for col in ("R", "S"):
            if col not in ct:
                ct[col] = 0
        lines += [
            "",
            f"## {drug}: mgrB status by phenotype (n = {len(lab):,})",
            "",
            "| mgrB status | R | S |",
            "|---|---|---|",
            *[f"| {i} | {int(r['R']):,} | {int(r['S']):,} |" for i, r in ct.iterrows()],
            "",
            f"### {drug}: rules with and without the mgrB check",
            "",
            "| Rule | Sensitivity | Specificity | FN | FP |",
            "|---|---|---|---|---|",
        ]
        for name, pred in tables.items():
            df = lab[["genome_id", "label"]].assign(predicted=lab["genome_id"].map(pred))
            m = metrics(df.dropna())
            lines.append(
                f"| {name} | {fmt(m['sensitivity'])} [{fmt(m['sens_lo'])}-{fmt(m['sens_hi'])}] "
                f"| {fmt(m['specificity'])} [{fmt(m['spec_lo'])}-{fmt(m['spec_hi'])}] "
                f"| {m['FN']} | {m['FP']} |"
            )

    args.out.mkdir(parents=True, exist_ok=True)
    status.to_csv(args.out / "mgrb_status.tsv", sep="\t", index=False)
    (args.out / "colistin_mgrb.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
