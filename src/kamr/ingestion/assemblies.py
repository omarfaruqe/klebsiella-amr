"""Download pilot assemblies from BV-BRC by genome_id, with validation and a provenance manifest.

Source order: the BV-BRC genome_sequence API, then BV-BRC FTP (https) as a fallback. In the
pilot the FTP route timed out from the server while the API answered, so the API goes first.
Each file is checked (FASTA format, total length within tolerance of BV-BRC genome_length)
and written atomically, so an interrupted run can simply be restarted.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

FTP_URL = "https://ftp.bvbrc.org/genomes/{gid}/{gid}.fna"
API_URL = "https://www.bv-brc.org/api/genome_sequence/?eq(genome_id,{gid})&limit(25000)"
LENGTH_TOLERANCE = 0.05
_local = threading.local()


def session() -> requests.Session:
    if not hasattr(_local, "session"):
        _local.session = requests.Session()
    return _local.session


def fasta_stats(text: str) -> tuple[int, int]:
    """Return (number of sequences, total sequence length) for FASTA text."""
    n_seqs, total = 0, 0
    for line in text.splitlines():
        if line.startswith(">"):
            n_seqs += 1
        else:
            total += len(line.strip())
    return n_seqs, total


def check_fasta(text: str, expected_length: float | None) -> str:
    """Return 'ok', or the reason the FASTA is not acceptable."""
    if not text.startswith(">"):
        return "not_fasta"
    n_seqs, total = fasta_stats(text)
    if n_seqs == 0 or total == 0:
        return "empty"
    if expected_length and abs(total - expected_length) / expected_length > LENGTH_TOLERANCE:
        return "length_mismatch"
    return "ok"


def fetch_fasta(gid: str) -> tuple[str, str]:
    """Return (fasta_text, source). Raises if both sources fail."""
    errors = []
    for source, url, headers, timeout in (
        ("bvbrc_api", API_URL.format(gid=gid), {"Accept": "application/dna+fasta"}, 120),
        ("bvbrc_ftp", FTP_URL.format(gid=gid), {}, 30),
    ):
        try:
            resp = session().get(url, headers=headers, timeout=timeout)
            resp.raise_for_status()
            if resp.text.startswith(">"):
                return resp.text, source
            errors.append(f"{source}: not FASTA")
        except requests.RequestException as exc:
            errors.append(f"{source}: {exc.__class__.__name__}")
    raise RuntimeError("; ".join(errors))


def download_one(gid: str, out_dir: Path, expected_length: float | None) -> dict:
    out = out_dir / f"{gid}.fna"
    row = {
        "genome_id": gid,
        "file": str(out),
        "source": "",
        "status": "",
        "n_seqs": 0,
        "length": 0,
        "sha256": "",
        "error": "",
    }
    if out.exists():
        text = out.read_text()
        row["source"] = "existing"
    else:
        try:
            text, row["source"] = fetch_fasta(gid)
        except RuntimeError as exc:
            row["status"], row["error"] = "download_failed", str(exc)
            return row
    row["status"] = check_fasta(text, expected_length)
    row["n_seqs"], row["length"] = fasta_stats(text)
    row["sha256"] = hashlib.sha256(text.encode()).hexdigest()
    if row["source"] != "existing" and row["status"] in ("ok", "length_mismatch"):
        tmp = out.with_suffix(".fna.tmp")
        tmp.write_text(text)
        tmp.replace(out)
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, default=Path("reports/cohort/pilot_sample.tsv"))
    parser.add_argument("--qc", type=Path, default=Path("data/curated/genome_qc.parquet"))
    parser.add_argument("--out", type=Path, default=Path("data/assemblies"))
    parser.add_argument("--manifests", type=Path, default=Path("snapshots"))
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--limit", type=int, default=0, help="only the first N genomes (testing)")
    args = parser.parse_args()

    ids = pd.read_csv(args.sample, sep="\t", dtype={"genome_id": str})["genome_id"].tolist()
    if args.limit:
        ids = ids[: args.limit]
    qc = pd.read_parquet(args.qc, columns=["genome_id", "genome_length"])
    expected = dict(zip(qc["genome_id"].astype(str), qc["genome_length"], strict=False))
    args.out.mkdir(parents=True, exist_ok=True)

    print(f"Downloading {len(ids)} assemblies with {args.workers} workers ...")
    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(download_one, g, args.out, expected.get(g)) for g in ids]
        for n, fut in enumerate(as_completed(futures), start=1):
            row = fut.result()
            rows.append(row)
            if row["status"] != "ok" or n % 50 == 0 or n == len(ids):
                print(f"  [{n:>4}/{len(ids)}] {row['genome_id']}: {row['status']} {row['error']}")

    table = pd.DataFrame(rows).sort_values("genome_id")
    day = dt.datetime.now(dt.UTC).strftime("%Y-%m-%d")
    args.manifests.mkdir(parents=True, exist_ok=True)
    manifest = args.manifests / f"assemblies-pilot-{day}.tsv"
    table.drop(columns=["file"]).to_csv(manifest, sep="\t", index=False)
    print("\nstatus counts:")
    print(table["status"].value_counts().to_string())
    print(f"\nwrote {manifest}")


if __name__ == "__main__":
    main()
