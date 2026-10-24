# Cohort genotype-phenotype concordance

- Generated (UTC): 2026-10-24T12:42:44+00:00
- Genomes evaluated: 6652 of 6652 in `reports/cohort/qc_pass_genomes.tsv`
- Labels: EUCAST v16.1 re-interpreted laboratory S/R (I and conflicts excluded)
- Colistin: broth microdilution only (ADR-003); colistin_all uses every method
- Rules: see src/kamr/analysis/rules.py

Sensitivity = resistant isolates the rule calls R; specificity = susceptible isolates the
rule calls S. 95% Wilson intervals in brackets.
Genome set: all QC-passing genomes with laboratory labels. Related isolates are not de-duplicated, so large outbreak studies weigh heavily; see the benchmark for grouped estimates.

| Drug | Tool | n (R/S) | Sensitivity | Specificity | FN | FP |
|---|---|---|---|---|---|---|
| meropenem | AMRFinderPlus | 5541 (1854/3687) | 96.7% [95.7%-97.4%] | 91.2% [90.2%-92.0%] | 62 | 326 |
| meropenem | AMRFinderPlus+regulators | 5541 (1854/3687) | 96.7% [95.7%-97.4%] | 91.2% [90.2%-92.0%] | 62 | 326 |
| meropenem | Kleborate | 5541 (1854/3687) | 96.5% [95.6%-97.2%] | 91.3% [90.3%-92.1%] | 65 | 322 |
| ciprofloxacin | AMRFinderPlus | 4640 (3441/1199) | 98.8% [98.4%-99.1%] | 88.2% [86.2%-89.9%] | 42 | 142 |
| ciprofloxacin | AMRFinderPlus+regulators | 4640 (3441/1199) | 99.0% [98.6%-99.2%] | 77.6% [75.2%-79.9%] | 36 | 268 |
| ciprofloxacin | Kleborate | 4640 (3441/1199) | 98.7% [98.3%-99.0%] | 88.4% [86.5%-90.1%] | 44 | 139 |
| gentamicin | AMRFinderPlus | 4760 (2113/2647) | 85.1% [83.6%-86.6%] | 94.9% [94.0%-95.7%] | 314 | 135 |
| gentamicin | AMRFinderPlus+regulators | 4760 (2113/2647) | 85.1% [83.6%-86.6%] | 94.9% [94.0%-95.7%] | 314 | 135 |
| gentamicin | Kleborate | 4760 (2113/2647) | 81.4% [79.6%-83.0%] | 95.8% [95.0%-96.5%] | 394 | 111 |
| colistin | AMRFinderPlus | 1528 (367/1161) | 28.3% [24.0%-33.2%] | 97.8% [96.7%-98.5%] | 263 | 26 |
| colistin | AMRFinderPlus+regulators | 1528 (367/1161) | 28.3% [24.0%-33.2%] | 97.8% [96.7%-98.5%] | 263 | 26 |
| colistin | Kleborate | 1528 (367/1161) | 39.5% [34.6%-44.6%] | 96.3% [95.0%-97.2%] | 222 | 43 |
| colistin_all | AMRFinderPlus | 1901 (418/1483) | 26.6% [22.5%-31.0%] | 97.9% [97.0%-98.5%] | 307 | 31 |
| colistin_all | AMRFinderPlus+regulators | 1901 (418/1483) | 26.6% [22.5%-31.0%] | 97.9% [97.0%-98.5%] | 307 | 31 |
| colistin_all | Kleborate | 1901 (418/1483) | 38.3% [33.7%-43.0%] | 96.4% [95.3%-97.2%] | 258 | 54 |
| tigecycline | AMRFinderPlus | 948 (261/687) | 2.7% [1.3%-5.4%] | 100.0% [99.4%-100.0%] | 254 | 0 |
| tigecycline | AMRFinderPlus+regulators | 948 (261/687) | 12.6% [9.1%-17.2%] | 84.4% [81.5%-86.9%] | 228 | 107 |
| tigecycline | Kleborate | 948 (261/687) | 0.0% [0.0%-1.5%] | 100.0% [99.4%-100.0%] | 261 | 0 |

## AMRFinderPlus (core) vs Kleborate prediction agreement

| Drug | Agreement |
|---|---|
| meropenem | 99.8% |
| ciprofloxacin | 99.9% |
| gentamicin | 97.8% |
| colistin | 87.8% |
| colistin_all | 88.7% |
| tigecycline | 99.3% |

Discordant genome-drug pairs (core AMRFinderPlus or Kleborate wrong): 1740, in `reports/cohort/discordant.tsv` with each tool's evidence for review.
