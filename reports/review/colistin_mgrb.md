# Colistin: targeted mgrB check

- Generated (UTC): 2026-10-24T14:14:00+00:00
- Genomes searched: 6,652 (tblastn failures: 0)
- Query: MgrB mgrB_WP_002911375.1.faa (47 aa), tblastn e-value <= 1e-5, hits >= 10 aa
- Labels: colistin = broth microdilution only (ADR-003); colistin_all = every method

## mgrB status across all searched genomes

| Status | Genomes |
|---|---|
| intact | 6,056 |
| split_contigs | 162 |
| partial | 153 |
| missense | 129 |
| premature_stop | 73 |
| absent | 34 |
| frameshift_or_indel | 27 |
| contig_break | 13 |
| insertion_assembled | 5 |

## mgrB missense variants (colistin, broth microdilution)

Missense variants are not counted as disrupted. Top 15 by genomes:

| Variant | R | S |
|---|---|---|
| M27K | 25 | 3 |
| T21N | 9 | 0 |
| K2E | 5 | 3 |
| N42Y | 7 | 0 |
| K43I | 7 | 0 |
| G37V | 5 | 0 |
| C28W | 4 | 1 |
| N25D | 0 | 3 |
| C39G | 3 | 0 |
| W6V | 1 | 2 |
| L4I | 1 | 2 |
| R5T | 1 | 2 |
| W6G | 2 | 0 |
| L4Y | 2 | 0 |
| R5G | 2 | 0 |

## colistin: mgrB status by phenotype (n = 1,528)

| mgrB status | R | S |
|---|---|---|
| absent | 14 | 2 |
| contig_break | 1 | 0 |
| frameshift_or_indel | 10 | 4 |
| insertion_assembled | 5 | 0 |
| intact | 144 | 1,107 |
| missense | 70 | 13 |
| partial | 58 | 16 |
| premature_stop | 15 | 1 |
| split_contigs | 50 | 18 |

### colistin: rules with and without the mgrB check

| Rule | Sensitivity | Specificity | FN | FP |
|---|---|---|---|---|
| AMRFinderPlus | 28.3% [24.0%-33.2%] | 97.8% [96.7%-98.5%] | 263 | 26 |
| AMRFinderPlus + mgrB disrupted | 56.7% [51.6%-61.7%] | 95.2% [93.8%-96.3%] | 159 | 56 |
| AMRFinderPlus + mgrB disrupted or missense | 67.3% [62.3%-71.9%] | 94.3% [92.8%-95.5%] | 120 | 66 |
| Kleborate | 39.5% [34.6%-44.6%] | 96.3% [95.0%-97.2%] | 222 | 43 |
| Kleborate + mgrB disrupted | 45.5% [40.5%-50.6%] | 95.3% [94.0%-96.4%] | 200 | 54 |
| Kleborate + mgrB disrupted or missense | 61.6% [56.5%-66.4%] | 94.5% [93.0%-95.7%] | 141 | 64 |
| mgrB disrupted alone | 41.7% [36.8%-46.8%] | 96.5% [95.2%-97.4%] | 214 | 41 |
| mgrB disrupted or missense alone | 60.8% [55.7%-65.6%] | 95.3% [94.0%-96.4%] | 144 | 54 |

## colistin_all: mgrB status by phenotype (n = 1,901)

| mgrB status | R | S |
|---|---|---|
| absent | 14 | 5 |
| contig_break | 2 | 3 |
| frameshift_or_indel | 11 | 4 |
| insertion_assembled | 5 | 0 |
| intact | 179 | 1,417 |
| missense | 72 | 16 |
| partial | 58 | 16 |
| premature_stop | 16 | 3 |
| split_contigs | 61 | 19 |

### colistin_all: rules with and without the mgrB check

| Rule | Sensitivity | Specificity | FN | FP |
|---|---|---|---|---|
| AMRFinderPlus | 26.6% [22.5%-31.0%] | 97.9% [97.0%-98.5%] | 307 | 31 |
| AMRFinderPlus + mgrB disrupted | 54.3% [49.5%-59.0%] | 95.4% [94.2%-96.4%] | 191 | 68 |
| AMRFinderPlus + mgrB disrupted or missense | 63.6% [58.9%-68.1%] | 94.5% [93.3%-95.6%] | 152 | 81 |
| Kleborate | 38.3% [33.7%-43.0%] | 96.4% [95.3%-97.2%] | 258 | 54 |
| Kleborate + mgrB disrupted | 43.8% [39.1%-48.6%] | 95.4% [94.2%-96.4%] | 235 | 68 |
| Kleborate + mgrB disrupted or missense | 58.1% [53.4%-62.8%] | 94.5% [93.3%-95.6%] | 175 | 81 |
| mgrB disrupted alone | 40.0% [35.4%-44.7%] | 96.6% [95.6%-97.4%] | 251 | 50 |
| mgrB disrupted or missense alone | 57.2% [52.4%-61.8%] | 95.5% [94.4%-96.5%] | 179 | 66 |
