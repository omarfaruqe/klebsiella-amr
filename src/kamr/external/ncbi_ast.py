"""Assemble an external colistin validation set from NCBI Pathogen Detection.

Leave-one-study-out within BV-BRC is internal validation: every study was drawn from the same
aggregator, curated the same way, and was available when the rule was built. This module builds
a set that is external in the only sense that matters -- isolates whose phenotypes we had no
access to when the rule was specified -- from NCBI Pathogen Detection.

Two NCBI tables are needed, because neither alone is sufficient:

  ncbi-pathogen-detect.pdbrowser.ast        submitted AST results, with MIC and testing standard
  ncbi-pathogen-detect.pdbrowser.isolates   isolate records, with the assembly accession

The AST table carries the MIC but not the assembly accession; the isolates table carries the
assembly accession but reduces AST to a bare phenotype. We need the MIC, because the submitted
phenotype was interpreted under whatever standard the submitter used, and the whole point of
section 2.2 is that this matters. So we join them. See SQL below.

Overlap is the central risk. BV-BRC ingests from NCBI, so a large share of NCBI isolates are
already in our cohort, and scoring the rule on those would restate the internal result with
extra steps. We exclude by BioSample accession and by assembly accession, and report how much
was excluded -- the overlap fraction is itself worth knowing, because it bounds how independent
any "external" NCBI-derived validation set can be.

Labels are re-derived from MIC under the same pinned EUCAST table used for the internal cohort,
never taken from the submitted phenotype, except for isolates with no usable MIC. Those are
retained in a separate stratum so the main result does not depend on foreign interpretation.

Outputs:
  data/external/ncbi/colistin_external.csv   held-out isolates with harmonised labels
  data/external/ncbi/assembly_manifest.txt   assembly accessions to download
  reports/external/ncbi_overlap.md           overlap and label-provenance report
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path

import pandas as pd

from kamr.curation.breakpoints import MIC_UNITS, interpret_mic, load_breakpoints, normalize_unit

DRUG = "colistin"

# Run this in BigQuery (a free sandbox project suffices) and export the result as CSV.
# Kept here so the query that produced the input file is under version control with the code
# that consumes it.
SQL = """
SELECT
  a.biosample_acc,
  i.asm_acc,
  a.bioproject_acc,
  a.collection_date,
  a.geo_loc_name,
  a.taxgroup_name,
  a.antibiotic,
  a.phenotype,
  a.mic,
  a.measurement_sign,
  -- the ast table has no units column: mic is mg/L by definition, disk_diffusion is mm.
  -- Emitted explicitly so the consumer's unit guard stays meaningful.
  CASE WHEN a.mic IS NOT NULL THEN 'mg/L' END AS measurement_units,
  a.standard
FROM `ncbi-pathogen-detect.pdbrowser.ast` AS a
JOIN `ncbi-pathogen-detect.pdbrowser.isolates` AS i
  ON a.target_acc = i.target_acc
WHERE a.antibiotic = 'colistin'
  AND a.taxgroup_name = 'Klebsiella pneumoniae'
  AND i.asm_acc IS NOT NULL
"""

# Submitted phenotype vocabulary (NCBI BioSample antibiogram controlled values). Used only when
# no usable MIC is present. "nonsusceptible" is deliberately not mapped to R: it conflates I and
# R, and folding it either way would misstate the denominator.
PHENOTYPE_MAP = {
    "resistant": "R",
    "susceptible": "S",
    "intermediate": "I",
    "susceptible-dose dependent": "I",
    "susceptible dose dependent": "I",
}

# Canonical field -> accepted spellings. The BigQuery tables and the browser downloads name the
# same fields differently ("biosample_acc" vs "BioSample"), so both are accepted.
COLUMNS = {
    "biosample": ("biosample_acc", "biosample"),
    "assembly": ("asm_acc", "assembly", "assembly_accession"),
    "bioproject": ("bioproject_acc", "bioproject"),
    "collection_date": ("collection_date",),
    "location": ("geo_loc_name", "location"),
    "submitted_phenotype": ("phenotype", "resistance_phenotype"),
    "mic": ("mic", "mic_mg/l"),
    "sign": ("measurement_sign", "sign"),
    "unit": ("measurement_units", "measurement_unit", "unit"),
    "standard": ("testing_standard", "standard"),
}


def read_export(path: Path) -> pd.DataFrame:
    """Read a BigQuery CSV export or an Isolates/AST Browser TSV, however its columns are cased."""
    sep = "\t" if path.suffix.lower() in (".tsv", ".tab") else ","
    df = pd.read_csv(path, sep=sep, dtype=str).fillna("")
    lookup = {
        c.strip().lower().replace(" ", "_").replace("(mg/l)", "").strip("_"): c for c in df.columns
    }
    out = pd.DataFrame(index=df.index)
    missing = []
    for name, spellings in COLUMNS.items():
        src = next((lookup[k] for k in spellings if k in lookup), None)
        if src is None:
            missing.append(spellings[0])
            out[name] = ""
        else:
            out[name] = df[src].astype(str).str.strip()
    if "assembly" in out and not out["assembly"].str.strip().any():
        raise SystemExit(
            "no assembly accessions in the export; the AST Browser download does not carry them, "
            "so use the documented SQL join against pdbrowser.isolates instead"
        )
    if missing:
        print(f"note: columns absent from export, treated as empty: {', '.join(sorted(missing))}")
    # NCBI's ast table carries no units column: its `mic` field is mg/L by definition, and disk
    # diffusion lives in a separate `disk_diffusion` field. Without this, every MIC would fail the
    # unit guard below and silently fall back to the submitter's own interpretation, defeating the
    # breakpoint harmonisation this module exists to perform. Stated loudly because it is an
    # assumption about an external schema, not something the file asserts.
    has_mic = out["mic"].str.strip() != ""
    if has_mic.any() and not (out["unit"].str.strip() != "").any():
        out.loc[has_mic, "unit"] = "mg/L"
        print(
            f"note: export carries no units column; treating {int(has_mic.sum()):,} MIC values as "
            "mg/L per the NCBI ast schema, where disk diffusion is a separate field"
        )
    return out


def cohort_identifiers(genome_dir: Path) -> tuple[set[str], set[str]]:
    """BioSample and assembly accessions of our BV-BRC cohort, from the newest snapshot."""
    days = sorted(p for p in genome_dir.iterdir() if p.is_dir()) if genome_dir.exists() else []
    if not days:
        raise SystemExit(f"no BV-BRC genome snapshot under {genome_dir}")
    path = days[-1] / "genome.jsonl"
    biosamples, assemblies = set(), set()
    with path.open() as fh:
        for line in fh:
            rec = json.loads(line)
            for key, sink in (
                ("biosample_accession", biosamples),
                ("assembly_accession", assemblies),
            ):
                val = str(rec.get(key) or "").strip()
                if val:
                    sink.add(val)
    print(f"cohort snapshot {path}: {len(biosamples):,} BioSamples, {len(assemblies):,} assemblies")
    return biosamples, assemblies


def strip_version(acc: str) -> str:
    """GCA_000240185.2 and GCA_000240185.1 are the same assembly for overlap purposes."""
    return acc.split(".", 1)[0].strip()


def harmonise(df: pd.DataFrame, breakpoints_path: Path) -> pd.DataFrame:
    """Re-derive S/I/R from MIC under the pinned EUCAST table, falling back to submitted calls."""
    bps, verified = load_breakpoints(breakpoints_path)
    bp = bps.get(DRUG)
    if bp is None or not bp.applies:
        raise SystemExit(f"no applicable {DRUG} breakpoint in {breakpoints_path}")
    if not verified:
        print("warning: breakpoint table contains unverified rows")

    labels, sources = [], []
    for row in df.itertuples(index=False):
        unit = normalize_unit(row.unit)
        call = interpret_mic(row.sign, row.mic, bp) if unit in MIC_UNITS else None
        if call is not None:
            labels.append(call)
            sources.append("mic_eucast")
        else:
            labels.append(PHENOTYPE_MAP.get(row.submitted_phenotype.strip().lower(), ""))
            sources.append("submitted" if labels[-1] else "unusable")
    out = df.copy()
    out["label"] = labels
    out["label_source"] = sources
    out["breakpoint"] = bp.version
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("export", type=Path, help="CSV/TSV from the documented BigQuery join")
    p.add_argument("--genomes", type=Path, default=Path("data/raw/bvbrc_genome"))
    bp_default = Path("config/breakpoints_eucast_v16.1.csv")
    p.add_argument("--breakpoints", type=Path, default=bp_default)
    p.add_argument("--out", type=Path, default=Path("data/external/ncbi"))
    p.add_argument("--reports", type=Path, default=Path("reports/external"))
    p.add_argument("--print-sql", action="store_true", help="print the BigQuery SQL and exit")
    a = p.parse_args()

    if a.print_sql:
        print(SQL.strip())
        return

    raw = read_export(a.export)
    n_rows = len(raw)
    raw = raw[raw["assembly"].str.strip() != ""]
    # one AST result per isolate; duplicates occur when a lab submits repeat measurements
    raw = raw.drop_duplicates(subset=["biosample", "assembly", "mic", "sign", "standard"])

    biosamples, assemblies = cohort_identifiers(a.genomes)
    assemblies_nv = {strip_version(x) for x in assemblies}
    in_cohort = raw["biosample"].isin(biosamples) | raw["assembly"].map(strip_version).isin(
        assemblies_nv
    )
    overlap_n = int(in_cohort.sum())
    held = raw[~in_cohort].copy()

    held = harmonise(held, a.breakpoints)
    usable = held[held["label"].isin(["R", "S", "I"])].copy()
    usable["y"] = (usable["label"] == "R").astype("uint8")

    a.out.mkdir(parents=True, exist_ok=True)
    a.reports.mkdir(parents=True, exist_ok=True)
    usable.to_csv(a.out / "colistin_external.csv", index=False)
    manifest = sorted(set(usable["assembly"]))
    (a.out / "assembly_manifest.txt").write_text("\n".join(manifest) + "\n")

    by_source = usable["label_source"].value_counts().to_dict()
    by_label = usable["label"].value_counts().to_dict()
    by_standard = usable["standard"].replace("", "(not stated)").value_counts().head(8)
    frac = overlap_n / len(raw) if len(raw) else float("nan")

    L = [
        "# External colistin set from NCBI Pathogen Detection",
        "",
        f"- Generated (UTC): {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        f"- Export: `{a.export}` ({n_rows:,} rows, {len(raw):,} after deduplication)",
        f"- Breakpoint: {usable['breakpoint'].iloc[0] if len(usable) else 'n/a'}",
        "",
        "## Overlap with the BV-BRC cohort",
        "",
        f"- Already in our cohort: {overlap_n:,} of {len(raw):,} ({frac:.1%})",
        f"- Held out as external: {len(held):,}",
        f"- With a usable label: {len(usable):,}",
        f"- Distinct assemblies to download: {len(manifest):,}",
        "",
        "Matched on BioSample accession or on assembly accession ignoring the version suffix.",
        "The overlap fraction bounds how independent any NCBI-derived validation set can be,",
        "given that BV-BRC ingests from NCBI.",
        "",
        "## Label provenance",
        "",
        "| Source | Isolates |",
        "|---|---|",
    ]
    for src in ("mic_eucast", "submitted"):
        L.append(f"| {src} | {by_source.get(src, 0):,} |")
    L += [
        "",
        "Labels from `mic_eucast` were re-derived from the submitted MIC under our pinned EUCAST",
        "table. Labels from `submitted` had no usable MIC and retain the submitter's own",
        "interpretation; they are reported as a separate stratum so the headline result does not",
        "rest on foreign breakpoint choices.",
        "",
        "## Composition",
        "",
        "| Label | Isolates |",
        "|---|---|",
    ]
    for lab in ("R", "I", "S"):
        L.append(f"| {lab} | {by_label.get(lab, 0):,} |")
    L += ["", "Testing standards as submitted:", "", "| Standard | Isolates |", "|---|---|"]
    L += [f"| {k} | {v:,} |" for k, v in by_standard.items()]
    text = "\n".join(L) + "\n"
    (a.reports / "ncbi_overlap.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
