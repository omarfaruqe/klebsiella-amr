"""Download a raw snapshot of BV-BRC AMR phenotype records for one taxon.

Raw records are written untouched (JSON Lines) under data/raw/, and a provenance
manifest is written under snapshots/ so it can be committed.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path

import requests

API = "https://www.bv-brc.org/api/genome_amr/"
TAXON_KPN = 573
PAGE_SIZE = 10_000
FIELDS = (
    "id",
    "genome_id",
    "genome_name",
    "taxon_id",
    "antibiotic",
    "resistant_phenotype",
    "measurement",
    "measurement_sign",
    "measurement_value",
    "measurement_unit",
    "laboratory_typing_method",
    "testing_standard",
    "testing_standard_year",
    "evidence",
    "source",
)


def build_query(taxon_id: int, limit: int, offset: int) -> str:
    """Return an RQL query string for one page of genome_amr records."""
    return (
        f"eq(taxon_id,{taxon_id})"
        f"&select({','.join(FIELDS)})"
        f"&sort(+id)"
        f"&limit({limit},{offset})"
    )


def fetch_all(taxon_id: int, page_size: int = PAGE_SIZE, timeout: int = 120) -> list[dict]:
    """Page through the API until a short page is returned."""
    records: list[dict] = []
    offset = 0
    with requests.Session() as session:
        session.headers.update({"Accept": "application/json"})
        while True:
            url = f"{API}?{build_query(taxon_id, page_size, offset)}"
            resp = session.get(url, timeout=timeout)
            resp.raise_for_status()
            page = resp.json()
            records.extend(page)
            print(f"  fetched {len(page):>6} records (offset {offset})")
            if len(page) < page_size:
                break
            offset += page_size
    return records


def write_snapshot(records: list[dict], out_root: Path, manifest_dir: Path, taxon_id: int) -> Path:
    """Write raw records plus a provenance manifest; return the data file path."""
    now = dt.datetime.now(dt.UTC)
    day = now.strftime("%Y-%m-%d")
    out_dir = out_root / day
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "genome_amr.jsonl"
    with out_file.open("w") as fh:
        for rec in records:
            fh.write(json.dumps(rec, sort_keys=True) + "\n")

    manifest = {
        "source": "BV-BRC",
        "endpoint": API,
        "query_template": build_query(taxon_id, PAGE_SIZE, 0),
        "taxon_id": taxon_id,
        "retrieved_utc": now.isoformat(timespec="seconds"),
        "records": len(records),
        "unique_genomes": len({r.get("genome_id") for r in records}),
        "file": str(out_file),
        "sha256": hashlib.sha256(out_file.read_bytes()).hexdigest(),
    }
    manifest_dir.mkdir(parents=True, exist_ok=True)
    manifest_file = manifest_dir / f"bvbrc-genome_amr-{day}.json"
    manifest_file.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {out_file}\nwrote {manifest_file}")
    return out_file


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--taxon", type=int, default=TAXON_KPN)
    parser.add_argument("--out", type=Path, default=Path("data/raw/bvbrc"))
    parser.add_argument("--manifests", type=Path, default=Path("snapshots"))
    args = parser.parse_args()

    print(f"Querying BV-BRC genome_amr for taxon {args.taxon} ...")
    records = fetch_all(args.taxon)
    write_snapshot(records, args.out, args.manifests, args.taxon)


if __name__ == "__main__":
    main()
