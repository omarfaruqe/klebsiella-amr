"""Download the held-out external assemblies from NCBI by accession.

Uses the NCBI Datasets v2 REST endpoint, one accession per request, because that keeps the run
resumable: an interrupted download leaves completed FASTA files in place and the next run skips
them. Each archive ships an md5sum.txt, which is checked before the FASTA is accepted, so a
truncated transfer fails loudly instead of producing a short genome that would later look like a
fragmented assembly and be silently excluded by QC.

NCBI permits three requests per second without an API key, so concurrency is deliberately low.

Outputs:
  data/external/ncbi/assemblies/<accession>.fna
  snapshots/ncbi-assemblies-<date>.json   provenance manifest with per-file SHA-256
  snapshots/ncbi-assemblies-<date>.tsv    genome_id/status manifest for kamr.pipeline.run_tools

The TSV deliberately matches the manifest shape run_tools already consumes, so the external set
goes through the identical tool invocations as the internal cohort rather than a parallel code
path that could drift from it.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import io
import json
import threading
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

from kamr.ingestion.assemblies import fasta_stats

API = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
    "?include_annotation_type=GENOME_FASTA"
)
# K. pneumoniae is ~5.2-5.8 Mb; this window only rejects truncation and gross mis-identification,
# not the fragmented-but-complete assemblies that QC handles downstream on its own terms.
MIN_LENGTH = 4_000_000
MAX_LENGTH = 8_000_000
WORKERS = 3
RETRIES = 4
_local = threading.local()


def session() -> requests.Session:
    if not hasattr(_local, "session"):
        _local.session = requests.Session()
    return _local.session


def expected_md5(archive: zipfile.ZipFile, member: str) -> str | None:
    """The md5 that md5sum.txt records for `member`, if it is listed."""
    try:
        text = archive.read("md5sum.txt").decode()
    except KeyError:
        return None
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[1].lstrip("./") == member.lstrip("./"):
            return parts[0]
    return None


def extract_fasta(blob: bytes, acc: str) -> tuple[str, str]:
    """Return (FASTA text, md5 status) from a Datasets archive, raising on anything unexpected."""
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        names = [
            n
            for n in z.namelist()
            if n.startswith(f"ncbi_dataset/data/{acc}/") and n.endswith("_genomic.fna")
        ]
        if len(names) != 1:
            raise ValueError(f"expected one genomic FASTA for {acc}, found {len(names)}")
        member = names[0]
        raw = z.read(member)
        want = expected_md5(z, member)
    got = hashlib.md5(raw).hexdigest()  # noqa: S324 - integrity check against NCBI, not security
    if want is None:
        status = "md5 not listed"
    elif want != got:
        raise ValueError(f"md5 mismatch for {acc}: archive says {want}, content is {got}")
    else:
        status = "md5 verified"
    return raw.decode(), status


def download_one(acc: str, out_dir: Path) -> dict:
    """Fetch one assembly, or report why it could not be used. Idempotent."""
    dest = out_dir / f"{acc}.fna"
    if dest.exists():
        text = dest.read_text()
        n_seqs, total = fasta_stats(text)
        return {
            "accession": acc,
            "status": "already present",
            "contigs": n_seqs,
            "length": total,
            "sha256": hashlib.sha256(dest.read_bytes()).hexdigest(),
        }
    last = ""
    for attempt in range(RETRIES):
        try:
            resp = session().get(API.format(acc=acc), timeout=180)
            if resp.status_code in (429, 500, 502, 503, 504):
                last = f"HTTP {resp.status_code}"
                time.sleep(2**attempt)
                continue
            resp.raise_for_status()
            text, md5_status = extract_fasta(resp.content, acc)
            n_seqs, total = fasta_stats(text)
            if not text.startswith(">"):
                return {"accession": acc, "status": "not FASTA"}
            if not MIN_LENGTH <= total <= MAX_LENGTH:
                return {
                    "accession": acc,
                    "status": f"length {total:,} outside {MIN_LENGTH:,}-{MAX_LENGTH:,}",
                    "contigs": n_seqs,
                    "length": total,
                }
            tmp = dest.with_suffix(".fna.part")
            tmp.write_text(text)
            tmp.replace(dest)
            return {
                "accession": acc,
                "status": md5_status,
                "contigs": n_seqs,
                "length": total,
                "sha256": hashlib.sha256(dest.read_bytes()).hexdigest(),
            }
        except (requests.RequestException, ValueError, zipfile.BadZipFile) as exc:
            last = f"{type(exc).__name__}: {exc}"
            time.sleep(2**attempt)
    return {"accession": acc, "status": f"failed after {RETRIES} attempts ({last})"}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--manifest", type=Path, default=Path("data/external/ncbi/assembly_manifest.txt")
    )
    p.add_argument("--out", type=Path, default=Path("data/external/ncbi/assemblies"))
    p.add_argument("--snapshots", type=Path, default=Path("snapshots"))
    p.add_argument("--workers", type=int, default=WORKERS)
    a = p.parse_args()

    accessions = [ln.strip() for ln in a.manifest.read_text().splitlines() if ln.strip()]
    a.out.mkdir(parents=True, exist_ok=True)
    print(f"{len(accessions):,} accessions from {a.manifest}")

    results = []
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        futures = {pool.submit(download_one, acc, a.out): acc for acc in accessions}
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result()
            results.append(res)
            if i % 25 == 0 or i == len(accessions):
                print(f"  {i:>4}/{len(accessions)}")

    ok = [r for r in results if r.get("sha256")]
    bad = [r for r in results if not r.get("sha256")]
    now = dt.datetime.now(dt.UTC)
    manifest = {
        "source": "NCBI Datasets v2",
        "endpoint": API.format(acc="<accession>"),
        "retrieved_utc": now.isoformat(timespec="seconds"),
        "requested": len(accessions),
        "downloaded": len(ok),
        "failed": len(bad),
        "files": sorted(results, key=lambda r: r["accession"]),
    }
    a.snapshots.mkdir(parents=True, exist_ok=True)
    stem = f"ncbi-assemblies-{now.strftime('%Y-%m-%d')}"
    out = a.snapshots / f"{stem}.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n")

    # genome_id/status manifest in the shape kamr.pipeline.run_tools expects, so the external
    # assemblies run through the same tool invocations as the internal cohort
    tsv = a.snapshots / f"{stem}.tsv"
    lines = ["genome_id\tstatus"]
    for r in sorted(results, key=lambda r: r["accession"]):
        lines.append(f"{r['accession']}\t{'ok' if r.get('sha256') else 'failed'}")
    tsv.write_text("\n".join(lines) + "\n")

    print(f"\nusable: {len(ok):,}   unusable: {len(bad):,}")
    for r in bad:
        print(f"  {r['accession']}: {r['status']}")
    print(f"wrote {out}\nwrote {tsv}")


if __name__ == "__main__":
    main()
