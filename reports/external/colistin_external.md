# External validation of the colistin rules (NCBI Pathogen Detection)

- Generated (UTC): 2026-10-25T04:59:55+00:00
- Genomes: 263 (81 R, 182 S)
- Query: MgrB mgrB_WP_002911375.1.faa, tblastn e-value <= 1e-5, hits >= 10 aa
- Rules imported unchanged from `kamr.analysis.mgrb.colistin_rules`; nothing was
  refitted, reselected or thresholded for this cohort.

These isolates are not in the BV-BRC cohort, matched on BioSample accession and on
assembly accession ignoring the version suffix. Their phenotypes were submitted to NCBI
as BioSample antibiograms and were not available to us when the rule was specified.

## mgrB status

| Status | Genomes |
|---|---|
| intact | 223 |
| premature_stop | 11 |
| partial | 10 |
| split_contigs | 9 |
| contig_break | 3 |
| absent | 3 |
| missense | 3 |
| insertion_assembled | 1 |

## All external isolates

| Rule | n | Sensitivity | Specificity | FN | FP |
|---|---|---|---|---|---|
| AMRFinderPlus | 263 | 50.6% [40.0%-61.2%] | 99.5% [97.0%-99.9%] | 40 | 1 |
| AMRFinderPlus + mgrB disrupted | 263 | 77.8% [67.6%-85.5%] | 99.5% [97.0%-99.9%] | 18 | 1 |
| AMRFinderPlus + mgrB disrupted or missense | 263 | 79.0% [68.9%-86.5%] | 99.5% [97.0%-99.9%] | 17 | 1 |
| Kleborate | 263 | 48.1% [37.6%-58.9%] | 100.0% [97.9%-100.0%] | 42 | 0 |
| Kleborate + mgrB disrupted | 263 | 59.3% [48.4%-69.3%] | 100.0% [97.9%-100.0%] | 33 | 0 |
| Kleborate + mgrB disrupted or missense | 263 | 60.5% [49.6%-70.4%] | 100.0% [97.9%-100.0%] | 32 | 0 |
| mgrB disrupted alone | 263 | 45.7% [35.3%-56.5%] | 100.0% [97.9%-100.0%] | 44 | 0 |
| mgrB disrupted or missense alone | 263 | 49.4% [38.8%-60.0%] | 100.0% [97.9%-100.0%] | 41 | 0 |

## Mechanism composition among the 81 resistant isolates

| Explanation | Isolates | Share |
|---|---|---|
| acquired determinant | 26 | 32.1% |
| mgrB disrupted | 22 | 27.2% |
| both | 15 | 18.5% |
| unexplained | 18 | 22.2% |

Acquired determinants detected by AMRFinderPlus, by gene:

| Determinant | Resistant isolates |
|---|---|
| pmrB_T157P | 10 |
| mcr-3.22 | 4 |
| mgrB_L4Ter | 4 |
| mcr-3.21 | 4 |
| mgrB_Q30Ter | 3 |
| mgrB_L18Ter | 2 |
| lapB_V43G | 1 |
| crrB_P151L | 1 |
| crrB_L94M | 1 |
| crrB_G183V | 1 |
| crrB_Q10L | 1 |
| mgrB_C28Y | 1 |
| mcr-1.1 | 1 |
| mgrB_L24insTer9 | 1 |
| mgrB_K3Ter | 1 |
| mgrB_V32G | 1 |
| mgrB_N42IfsTer30 | 1 |
| phoQ_D438N | 1 |
| pmrA_G53C | 1 |
| pmrA_G53S | 1 |
| pmrB_S85R | 1 |

## Restricted to MIC-derived labels (n = 251)

Labels re-derived from the submitted MIC under our pinned EUCAST table, excluding
isolates that retain the submitter's own interpretation.

| Rule | n | Sensitivity | Specificity | FN | FP |
|---|---|---|---|---|---|
| AMRFinderPlus | 251 | 46.4% [35.1%-58.0%] | 99.5% [97.0%-99.9%] | 37 | 1 |
| AMRFinderPlus + mgrB disrupted | 251 | 78.3% [67.2%-86.4%] | 99.5% [97.0%-99.9%] | 15 | 1 |
| AMRFinderPlus + mgrB disrupted or missense | 251 | 79.7% [68.8%-87.5%] | 99.5% [97.0%-99.9%] | 14 | 1 |
| Kleborate | 251 | 43.5% [32.4%-55.2%] | 100.0% [97.9%-100.0%] | 39 | 0 |
| Kleborate + mgrB disrupted | 251 | 56.5% [44.8%-67.6%] | 100.0% [97.9%-100.0%] | 30 | 0 |
| Kleborate + mgrB disrupted or missense | 251 | 58.0% [46.2%-68.9%] | 100.0% [97.9%-100.0%] | 29 | 0 |
| mgrB disrupted alone | 251 | 52.2% [40.6%-63.5%] | 100.0% [97.9%-100.0%] | 33 | 0 |
| mgrB disrupted or missense alone | 251 | 55.1% [43.4%-66.2%] | 100.0% [97.9%-100.0%] | 31 | 0 |
