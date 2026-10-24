# Do mgrB contig breaks look like IS insertions?

- Generated (UTC): 2026-10-24T18:18:53+00:00
- Assembly-inferred disruptions examined: 328; matched intact controls: 328
- Flank: up to 600 bp of contig end adjacent to mgrB, minimum 50 bp; copies counted at >=95.0% identity over >=80% of the flank, excluding +/-2000 bp around mgrB
- Controls are random other contig ends truncated to the same length as each genome's mgrB flank

## 1. Position of mgrB relative to the contig end

| Genomes | median distance (bp) | within 100 bp | within 1 kb |
|---|---|---|---|
| disrupted (assembly-inferred) (n=328) | 79 | 63% | 87% |
| intact mgrB (n=328) | 80,470 | 0% | 2% |

If breaks were incidental to fragmentation, mgrB would sit at a similar distance
from contig ends in both groups.

## 2. Are mgrB-adjacent contig ends repetitive within their own genome?

Flank lengths analysed: median 95 bp (range 50-600); 73 of 328 genomes had under 50 bp between mgrB and the contig end and were not analysed here. Controls are length-matched.

| Sequence | n | median copies elsewhere | >=1 copy | >=3 copies |
|---|---|---|---|---|
| mgrB-adjacent contig end | 255 | 5 | 82% | 72% |
| random other contig end (control) | 255 | 2 | 66% | 35% |

Contig ends are repeat-enriched in general, so the control row is the
comparison of interest, not the absolute copy number.

## 3. Do the same sequences recur at mgrB across independent genomes?

249 of 255 mgrB-adjacent flanks (98%) match a flank from a different genome at >=90% identity over >=80% of the shorter sequence (minimum 50 bp).

Idiosyncratic assembly artefacts would not be shared between isolates; mobile elements
would. This does not identify which element is involved, which requires either an IS
reference database or long-read assemblies of these isolates.

Per-genome detail: `reports/review/mgrb_flanks.tsv`.
