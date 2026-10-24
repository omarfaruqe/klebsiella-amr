# mgrB result: robustness checks

- Generated (UTC): 2026-10-24T16:03:56+00:00
- Colistin broth-microdilution set: 1,528 genomes (367 R)

## 1. Is the association confounded by assembly quality?

Assembly fragmentation by mgrB call (all genomes in this set):

| mgrB call | n | median contigs | median N50 |
|---|---|---|---|
| disrupted | 194 | 116 | 219,734 |
| intact | 1,251 | 104 | 256,611 |
| missense | 83 | 108 | 269,763 |

Colistin resistance rate by assembly-quality quartile and mgrB status. If the
association were an artefact of fragmentation, the gap would close within quartiles.

| Contig quartile | mgrB intact %R (n) | mgrB disrupted %R (n) |
|---|---|---|
| Q1 fewest | 4% (337) | 81% (32) |
| Q2 | 8% (315) | 69% (48) |
| Q3 | 17% (292) | 83% (65) |
| Q4 most | 18% (307) | 82% (49) |

Disruption calls split by how they were detected:

| Detection | n | %R |
|---|---|---|
| Assembly-structure (split/break/partial) | 143 | 76% |
| Sequence-level (stop/indel/absent/insertion) | 51 | 86% |
| Intact mgrB | 1251 | 12% |

## 2. Leave-one-study-out estimate for the missense tier

Missense variants are selected on training studies only (>=3 genomes, >50% resistant)
and applied to the held-out study, removing the circularity of selecting variants on
the same data used to report sensitivity.

| Held-out study | R | S | Sens: disrupted only | Sens: + missense (OOS) | Specificity | variants |
|---|---|---|---|---|---|---|
| PRJEB10018 | 162 | 747 | 39.5% | 39.5% | 96.0% | 3 |
| PRJEB10561 | 17 | 14 | 70.6% | 70.6% | 100.0% | 7 |
| PRJEB22890 | 63 | 93 | 50.8% | 52.4% | 89.2% | 6 |
| PRJNA288601 | 20 | 110 | 55.0% | 65.0% | 99.1% | 8 |
| PRJNA292904 | 7 | 15 | 57.1% | 57.1% | 100.0% | 8 |
| PRJNA307517 | 19 | 4 | 63.2% | 63.2% | 100.0% | 8 |
| PRJNA504600 | 27 | 0 | 48.1% | 48.1% | nan% | 8 |

**Disrupted only (nothing fitted): 47.0% (41.5%-52.5%).**
**With missense selected out-of-sample: 47.9% (42.5%-53.4%).**
Pooled over 315 resistant isolates in 7 studies.

In-sample equivalents from the main analysis: 56.7% (disrupted) and 67.3%
(disrupted or missense). The gap between the in-sample and out-of-sample missense
figures is the cost of selecting variants on the same data. Per-study sensitivity
also varies widely in the disrupted-only column, which is fitted to nothing: the
mechanism mix of colistin resistance genuinely differs between collections.

## 3. Gain against the stronger baseline

Kleborate already reports MgrB truncations, so it, not AMRFinderPlus, is the
conservative comparator for the headline claim.

| Baseline | Sensitivity | With mgrB check | Absolute gain |
|---|---|---|---|
| AMRFinderPlus | 28.3% | 67.3% | +39.0 points |
| Kleborate | 39.5% | 61.6% | +22.1 points |

The gain against Kleborate is the figure that should carry the headline.
