"""Test whether mgrB contig breaks look like insertion-sequence (IS) insertions.

The mechanistic claim behind the split_contigs, contig_break and partial categories is that an IS
element has landed in mgrB and the short-read assembler could not span it, so the contig ends
inside the gene. That claim makes predictions which generic assembly fragmentation does not:

1. Position. The break should sit immediately at the mgrB fragment, not at a random distance.
   Compared against intact-mgrB genomes, where mgrB should lie far from any contig end.
2. Repetitiveness within the genome. IS elements are present in multiple copies, so sequence at
   the mgrB-adjacent contig end should have high-identity copies elsewhere in the same assembly.
3. Recurrence across genomes. The same few elements should appear at mgrB in independent isolates,
   so flanking sequences should cluster rather than be idiosyncratic.

Prediction 2 needs a control, because contig ends are repeat-enriched in general: we compare
mgrB-adjacent contig ends against randomly chosen other contig ends in the same assemblies.

Outputs: reports/review/mgrb_flanks.md, reports/review/mgrb_flanks.tsv
"""

from __future__ import annotations

import argparse
import datetime as dt
import random
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from kamr.analysis.mgrb import MIN_HIT_AA, read_hits

ASSEMBLY_INFERRED = {"split_contigs", "contig_break", "partial"}
FLANK_BP = 600
MIN_FLANK_BP = 50  # most breaks sit <100 bp from mgrB, so short flanks must be usable
MIN_COPY_ID = 95.0
MIN_COPY_FRACTION = 0.8  # alignment must cover this much of the flank, which is often short
AWAY_FROM_MGRB_BP = 2000


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def best_fragments(hits: pd.DataFrame) -> pd.DataFrame:
    """Alignment fragments of >= MIN_HIT_AA, with contig span normalised to lo/hi."""
    h = hits[hits["length"] >= MIN_HIT_AA].copy()
    if h.empty:
        return h
    h["lo"] = h[["sstart", "send"]].min(axis=1)
    h["hi"] = h[["sstart", "send"]].max(axis=1)
    return h


def distance_to_contig_end(frag: pd.Series) -> int:
    """Base pairs from the nearer edge of this mgrB fragment to the nearer contig end."""
    return int(min(frag["lo"] - 1, frag["slen"] - frag["hi"]))


def flank_sequence(seq: str, frag: pd.Series, flank_bp: int = FLANK_BP) -> tuple[str, str]:
    """Sequence at the contig end nearest the mgrB fragment; returns (side, sequence)."""
    left, right = int(frag["lo"] - 1), int(frag["slen"] - frag["hi"])
    if left <= right:
        return "left", seq[: int(frag["lo"]) - 1][-flank_bp:]
    return "right", seq[int(frag["hi"]) :][:flank_bp]


def blastn_counts(
    query: Path, subject: Path, exclude: tuple[str, int, int] | None, min_bp: int
) -> int:
    """High-identity copies of the query in the subject assembly, outside the mgrB region."""
    cmd = [
        "blastn",
        "-query",
        str(query),
        "-subject",
        str(subject),
        "-outfmt",
        "6 sseqid sstart send pident length",
        "-evalue",
        "1e-10",
        "-dust",
        "no",
    ]
    out = subprocess.run(cmd, capture_output=True, text=True, check=False)
    n = 0
    for line in out.stdout.splitlines():
        sseqid, sstart, send, pident, length = line.split("\t")
        if float(pident) < MIN_COPY_ID or int(length) < min_bp:
            continue
        if exclude is not None:
            contig, lo, hi = exclude
            mid = (int(sstart) + int(send)) / 2
            if sseqid == contig and lo - AWAY_FROM_MGRB_BP <= mid <= hi + AWAY_FROM_MGRB_BP:
                continue
        n += 1
    return n


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n{seq}\n")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mgrb", type=Path, default=Path("reports/review/mgrb_status.tsv"))
    p.add_argument("--hits", type=Path, default=Path("results/cohort/mgrb_hits"))
    p.add_argument("--assemblies", type=Path, default=Path("data/assemblies"))
    p.add_argument("--out", type=Path, default=Path("reports/review"))
    p.add_argument("--controls-per-genome", type=int, default=1)
    p.add_argument("--seed", type=int, default=20261024)
    a = p.parse_args()

    rng = random.Random(a.seed)
    status = pd.read_csv(a.mgrb, sep="\t", dtype={"genome_id": str})
    disrupted = status[status.mgrb_status.isin(ASSEMBLY_INFERRED)].genome_id.tolist()
    intact = status[status.mgrb_status == "intact"].genome_id.tolist()
    rng.shuffle(intact)
    intact = intact[: len(disrupted)]

    rows, flanks, controls = [], {}, {}
    tmp = Path(tempfile.mkdtemp())
    for group, ids in (("disrupted", disrupted), ("intact", intact)):
        for gid in ids:
            hp, ap = a.hits / f"{gid}.tsv", a.assemblies / f"{gid}.fna"
            if not hp.exists() or not ap.exists():
                continue
            frags = best_fragments(read_hits(hp))
            if frags.empty:
                continue
            frag = frags.sort_values("length", ascending=False).iloc[0]
            seqs = read_fasta(ap)
            seq = seqs.get(str(frag["sseqid"]))
            if seq is None:
                continue
            dist = distance_to_contig_end(frag)
            side, flank = flank_sequence(seq, frag)
            row = {
                "genome_id": gid,
                "group": group,
                "mgrb_status": status.loc[status.genome_id == gid, "mgrb_status"].iloc[0],
                "contig": str(frag["sseqid"]),
                "contig_len": int(frag["slen"]),
                "dist_to_contig_end": dist,
                "side": side,
                "flank_bp": len(flank),
            }
            if group == "disrupted" and len(flank) >= MIN_FLANK_BP:
                flanks[gid] = flank
                min_bp = max(MIN_FLANK_BP, int(len(flank) * MIN_COPY_FRACTION))
                qf = tmp / "q.fna"
                write_fasta(qf, {gid: flank})
                row["copies_in_genome"] = blastn_counts(
                    qf, ap, (str(frag["sseqid"]), int(frag["lo"]), int(frag["hi"])), min_bp
                )
                # Control: a random other contig end, truncated to the same length as this
                # genome's mgrB flank, so the copy-count comparison is not length-biased.
                others = [c for c in seqs if c != str(frag["sseqid"]) and len(seqs[c]) > len(flank)]
                if others:
                    c = rng.choice(others)
                    ctrl = seqs[c][-len(flank) :]
                    controls[gid] = ctrl
                    write_fasta(qf, {f"{gid}_ctrl": ctrl})
                    row["control_copies_in_genome"] = blastn_counts(
                        qf, ap, (c, len(seqs[c]) - len(flank), len(seqs[c])), min_bp
                    )

            rows.append(row)
    t = pd.DataFrame(rows)
    a.out.mkdir(parents=True, exist_ok=True)
    t.to_csv(a.out / "mgrb_flanks.tsv", sep="\t", index=False)

    # cross-genome recurrence of mgrB flanks
    shared = 0
    if len(flanks) >= 2:
        fa = tmp / "flanks.fna"
        write_fasta(fa, flanks)
        out = subprocess.run(
            [
                "blastn",
                "-query",
                str(fa),
                "-subject",
                str(fa),
                "-outfmt",
                "6 qseqid sseqid pident length",
                "-evalue",
                "1e-10",
                "-dust",
                "no",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        partners: dict[str, set[str]] = {}
        for line in out.stdout.splitlines():
            q, s, pid, ln = line.split("\t")
            shorter = min(len(flanks[q]), len(flanks[s]))
            if q != s and float(pid) >= 90 and int(ln) >= max(MIN_FLANK_BP, shorter * 0.8):
                partners.setdefault(q, set()).add(s)
        shared = sum(1 for g in flanks if partners.get(g))

    d = t[t.group == "disrupted"]
    i = t[t.group == "intact"]
    analysed = d[d.flank_bp >= MIN_FLANK_BP]
    L = [
        "# Do mgrB contig breaks look like IS insertions?",
        "",
        f"- Generated (UTC): {dt.datetime.now(dt.UTC).isoformat(timespec='seconds')}",
        f"- Assembly-inferred disruptions examined: {len(d):,}; "
        f"matched intact controls: {len(i):,}",
        f"- Flank: up to {FLANK_BP} bp of contig end adjacent to mgrB, minimum {MIN_FLANK_BP} bp; "
        f"copies counted at >={MIN_COPY_ID}% identity over >={MIN_COPY_FRACTION:.0%} of the flank, "
        f"excluding +/-{AWAY_FROM_MGRB_BP} bp around mgrB",
        "- Controls are random other contig ends truncated to the same length as each genome's "
        "mgrB flank",
        "",
        "## 1. Position of mgrB relative to the contig end",
        "",
        "| Genomes | median distance (bp) | within 100 bp | within 1 kb |",
        "|---|---|---|---|",
    ]
    for name, sub in (("disrupted (assembly-inferred)", d), ("intact mgrB", i)):
        if len(sub):
            L.append(
                f"| {name} (n={len(sub)}) | {sub.dist_to_contig_end.median():,.0f} | "
                f"{(sub.dist_to_contig_end <= 100).mean():.0%} | "
                f"{(sub.dist_to_contig_end <= 1000).mean():.0%} |"
            )
    L += [
        "",
        "If breaks were incidental to fragmentation, mgrB would sit at a similar distance",
        "from contig ends in both groups.",
        "",
        "## 2. Are mgrB-adjacent contig ends repetitive within their own genome?",
        "",
        f"Flank lengths analysed: median {analysed.flank_bp.median():.0f} bp "
        f"(range {analysed.flank_bp.min():.0f}-{analysed.flank_bp.max():.0f}); "
        f"{len(d) - len(analysed)} of {len(d)} genomes had under {MIN_FLANK_BP} bp "
        "between mgrB and the contig end and were not analysed here. Controls are "
        "length-matched.",
        "",
    ]
    if "copies_in_genome" in d and d.copies_in_genome.notna().any():
        cc = d.copies_in_genome.dropna()
        ctrl = (
            d.control_copies_in_genome.dropna()
            if "control_copies_in_genome" in d
            else pd.Series(dtype=float)
        )
        L += [
            "| Sequence | n | median copies elsewhere | >=1 copy | >=3 copies |",
            "|---|---|---|---|---|",
            f"| mgrB-adjacent contig end | {len(cc)} | {cc.median():.0f} | "
            f"{(cc >= 1).mean():.0%} | {(cc >= 3).mean():.0%} |",
        ]
        if len(ctrl):
            L.append(
                f"| random other contig end (control) | {len(ctrl)} | {ctrl.median():.0f} | "
                f"{(ctrl >= 1).mean():.0%} | {(ctrl >= 3).mean():.0%} |"
            )
        L += [
            "",
            "Contig ends are repeat-enriched in general, so the control row is the",
            "comparison of interest, not the absolute copy number.",
        ]
    L += [
        "",
        "## 3. Do the same sequences recur at mgrB across independent genomes?",
        "",
        f"{shared:,} of {len(flanks):,} mgrB-adjacent flanks ({shared / max(len(flanks), 1):.0%}) "
        f"match a flank from a different genome at >=90% identity over >=80% of the shorter "
        f"sequence (minimum {MIN_FLANK_BP} bp).",
        "",
        "Idiosyncratic assembly artefacts would not be shared between isolates; mobile elements",
        "would. This does not identify which element is involved, which requires either an IS",
        "reference database or long-read assemblies of these isolates.",
        "",
        "Per-genome detail: `reports/review/mgrb_flanks.tsv`.",
    ]
    (a.out / "mgrb_flanks.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
