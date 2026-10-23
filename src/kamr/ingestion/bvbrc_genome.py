"""Download BV-BRC genome metadata for every labelled genome in the cohort.

Writes a raw JSONL snapshot under data/raw/bvbrc_genome/<date>/ and a provenance manifest under
snapshots/, in the same style as bvbrc_amr.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from collections.abc import Iterable
from pathlib import Path

import pandas as pd
import requests

API = "https://www.bv-brc.org/api/genome/"
BATCH = 200
FIELDS = (
    "genome_id",
    "genome_name",
    "species",
    "taxon_id",
    "genome_status",
    "genome_quality",
    "assembly_accession",
    "bioproject_accession",
    "biosample_accession",
    "sra_accession",
    "contigs",
    "contig_n50",
    "genome_length",
    "gc_content",
    "checkm_completeness",
    "checkm_contamination",
    "mlst",
    "collection_year",
    "collection_date",
    "isolation_country",
    "isolation_source",
    "host_name",
)


def build_query(ids: Iterable[str]) -> str:
    ids = list(ids)
    return f"in(genome_id,({','.join(ids)}))&select({','.join(FIELDS)})&limit({len(ids) + 10})"


def batches(items: list[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(items), size):
        yield items[start : start + size]


def fetch(ids: list[str], size: int = BATCH, timeout: int = 120) -> list[dict]:
    records: list[dict] = []
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/rqlquery+x-www-form-urlencoded",
    }
    with requests.Session() as session:
        for n, chunk in enumerate(batches(ids, size), start=1):
            resp = session.post(API, data=build_query(chunk), headers=headers, timeout=timeout)
            resp.raise_for_status()
            page = resp.json()
            records.extend(page)
            print(f"  batch {n:>3}: requested {len(chunk):>4}, received {len(page):>4}")
    return records


def labelled_genome_ids(labels_path: Path) -> list[str]:
    labels = pd.read_parquet(labels_path)
    keep = labels[labels["label"].isin(["R", "S", "I"])]
    return sorted(keep["genome_id"].astype(str).unique())


def write_snapshot(records: list[dict], requested: int, out_root: Path, manifests: Path) -> Path:
    now = dt.datetime.now(dt.UTC)
    day = now.strftime("%Y-%m-%d")
    out_dir = out_root / day
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "genome.jsonl"
    with out_file.open("w") as fh:
        for rec in sorted(records, key=lambda r: str(r.get("genome_id"))):
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
    manifest = {
        "source": "BV-BRC",
        "endpoint": API,
        "query_template": build_query(["<genome_id>", "..."]),
        "batch_size": BATCH,
        "retrieved_utc": now.isoformat(timespec="seconds"),
        "requested_genomes": requested,
        "records": len(records),
        "file": str(out_file),
        "sha256": hashlib.sha256(out_file.read_bytes()).hexdigest(),
    }
    manifests.mkdir(parents=True, exist_ok=True)
    manifest_file = manifests / f"bvbrc-genome-{day}.json"
    manifest_file.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {out_file}\nwrote {manifest_file}")
    return out_file


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, default=Path("data/curated/cohort_labels.parquet"))
    parser.add_argument("--out", type=Path, default=Path("data/raw/bvbrc_genome"))
    parser.add_argument("--manifests", type=Path, default=Path("snapshots"))
    args = parser.parse_args()

    ids = labelled_genome_ids(args.labels)
    print(f"Fetching BV-BRC genome metadata for {len(ids):,} labelled genomes ...")
    records = fetch(ids)
    missing = len(set(ids) - {str(r.get("genome_id")) for r in records})
    print(f"received {len(records):,} records; {missing:,} requested genomes not returned")
    write_snapshot(records, len(ids), args.out, args.manifests)


if __name__ == "__main__":
    main()
