"""Run AMRFinderPlus and Kleborate on every downloaded pilot assembly, in parallel.

Each genome is one job: AMRFinderPlus (nucleotide mode, --plus, K. pneumoniae) then Kleborate
(-p kpsc). Outputs are written to a temporary path and moved into place on success, so a
restarted run skips finished genomes and redoes only incomplete ones.

Layout under --results (default results/pilot/):
  amrfinder/<genome_id>.tsv
  kleborate/<genome_id>/            Kleborate output tables
  logs/<genome_id>.<tool>.log       stdout + stderr of each tool
  run_log.tsv                       one line per job: status and seconds per tool
  versions.txt                      tool and database versions used
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

ORGANISM = "Klebsiella_pneumoniae"


def tool_versions() -> str:
    def first_line(cmd: list[str]) -> str:
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, check=False)
            return (out.stdout or out.stderr).strip().splitlines()[-1]
        except (OSError, IndexError):
            return "unavailable"

    db = subprocess.run(["amrfinder", "-l"], capture_output=True, text=True, check=False)
    db_line = next(
        (ln for ln in (db.stdout + db.stderr).splitlines() if ln.startswith("Database version")),
        "Database version: unknown",
    )
    return "\n".join(
        [
            f"amrfinder {first_line(['amrfinder', '--version'])}",
            f"amrfinder {db_line}",
            f"kleborate {first_line(['kleborate', '--version'])}",
        ]
    )


def run_logged(cmd: list[str], log: Path) -> tuple[bool, float]:
    start = time.monotonic()
    with log.open("w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, check=False)
    return proc.returncode == 0, round(time.monotonic() - start, 1)


def run_genome(gid: str, fasta: Path, results: Path, threads: int) -> dict:
    row = {
        "genome_id": gid,
        "amrfinder": "skipped",
        "amrfinder_s": 0.0,
        "kleborate": "skipped",
        "kleborate_s": 0.0,
    }
    logs = results / "logs"

    amr_out = results / "amrfinder" / f"{gid}.tsv"
    if not amr_out.exists():
        tmp = amr_out.with_suffix(".tsv.tmp")
        ok, secs = run_logged(
            [
                "amrfinder",
                "-n",
                str(fasta),
                "--organism",
                ORGANISM,
                "--plus",
                "--threads",
                str(threads),
                "--name",
                gid,
                "-o",
                str(tmp),
            ],
            logs / f"{gid}.amrfinder.log",
        )
        row["amrfinder_s"] = secs
        if ok and tmp.exists():
            tmp.replace(amr_out)
            row["amrfinder"] = "ok"
        else:
            row["amrfinder"] = "failed"

    kl_out = results / "kleborate" / gid
    if not (kl_out / ".done").exists():
        tmp_dir = results / "kleborate" / f".{gid}.tmp"
        shutil.rmtree(tmp_dir, ignore_errors=True)
        ok, secs = run_logged(
            ["kleborate", "-a", str(fasta), "-o", str(tmp_dir), "-p", "kpsc"],
            logs / f"{gid}.kleborate.log",
        )
        row["kleborate_s"] = secs
        if ok and any(tmp_dir.glob("*.txt")):
            (tmp_dir / ".done").touch()
            shutil.rmtree(kl_out, ignore_errors=True)
            tmp_dir.replace(kl_out)
            row["kleborate"] = "ok"
        else:
            row["kleborate"] = "failed"
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="snapshots/assemblies-pilot-<date>.tsv")
    parser.add_argument("--assemblies", type=Path, default=Path("data/assemblies"))
    parser.add_argument("--results", type=Path, default=Path("results/pilot"))
    parser.add_argument("--jobs", type=int, default=8, help="genomes processed at once")
    parser.add_argument("--threads", type=int, default=4, help="AMRFinderPlus threads per job")
    parser.add_argument("--limit", type=int, default=0, help="only the first N genomes (testing)")
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t", dtype={"genome_id": str})
    ids = manifest.loc[manifest["status"] == "ok", "genome_id"].tolist()
    if args.limit:
        ids = ids[: args.limit]
    for sub in ("amrfinder", "kleborate", "logs"):
        (args.results / sub).mkdir(parents=True, exist_ok=True)
    (args.results / "versions.txt").write_text(tool_versions() + "\n")

    run_log = args.results / "run_log.tsv"
    if not run_log.exists():
        run_log.write_text("genome_id\tamrfinder\tamrfinder_s\tkleborate\tkleborate_s\n")

    print(f"Running {len(ids)} genomes: {args.jobs} at once, {args.threads} threads each")
    start = time.monotonic()
    failed = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(run_genome, g, args.assemblies / f"{g}.fna", args.results, args.threads): g
            for g in ids
        }
        for n, fut in enumerate(as_completed(futures), start=1):
            row = fut.result()
            failed += "failed" in (row["amrfinder"], row["kleborate"])
            with run_log.open("a") as fh:
                fh.write("\t".join(str(row[k]) for k in row) + "\n")
            elapsed = (time.monotonic() - start) / 60
            print(
                f"  [{n:>4}/{len(ids)}] {row['genome_id']}: amrfinder {row['amrfinder']} "
                f"({row['amrfinder_s']}s), kleborate {row['kleborate']} ({row['kleborate_s']}s)"
                f"  | {elapsed:.1f} min, {failed} failed",
                flush=True,
            )
    print(
        f"\nfinished in {(time.monotonic() - start) / 60:.1f} min; {failed} genomes with failures"
    )


if __name__ == "__main__":
    main()
