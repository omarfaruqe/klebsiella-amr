# Pilot genotype-phenotype concordance

- Generated (UTC): 2026-10-24T01:44:46+00:00
- Genomes evaluated: 523 of 524 in the pilot sample (excluded, incomplete tool output: 573.46141)
- Labels: EUCAST v16.1 re-interpreted laboratory S/R (I and conflicts excluded)
- Colistin: broth microdilution only (ADR-003); colistin_all uses every method
- Rules: see src/kamr/analysis/rules.py

Sensitivity = resistant isolates the rule calls R; specificity = susceptible isolates the
rule calls S. 95% Wilson intervals in brackets. The pilot over-samples rare STs, so these
are tool-validation figures, not population estimates.

| Drug | Tool | n (R/S) | Sensitivity | Specificity | FN | FP |
|---|---|---|---|---|---|---|
| meropenem | AMRFinderPlus | 428 (169/259) | 92.3% [87.3%-95.4%] | 85.7% [80.9%-89.5%] | 13 | 37 |
| meropenem | AMRFinderPlus+regulators | 428 (169/259) | 92.3% [87.3%-95.4%] | 85.7% [80.9%-89.5%] | 13 | 37 |
| meropenem | Kleborate | 428 (169/259) | 92.3% [87.3%-95.4%] | 86.1% [81.4%-89.8%] | 13 | 36 |
| ciprofloxacin | AMRFinderPlus | 314 (190/124) | 94.7% [90.6%-97.1%] | 92.7% [86.8%-96.1%] | 10 | 9 |
| ciprofloxacin | AMRFinderPlus+regulators | 314 (190/124) | 95.8% [91.9%-97.9%] | 79.0% [71.0%-85.3%] | 8 | 26 |
| ciprofloxacin | Kleborate | 314 (190/124) | 94.7% [90.6%-97.1%] | 93.5% [87.8%-96.7%] | 10 | 8 |
| gentamicin | AMRFinderPlus | 397 (154/243) | 91.6% [86.1%-95.0%] | 94.2% [90.6%-96.5%] | 13 | 14 |
| gentamicin | AMRFinderPlus+regulators | 397 (154/243) | 91.6% [86.1%-95.0%] | 94.2% [90.6%-96.5%] | 13 | 14 |
| gentamicin | Kleborate | 397 (154/243) | 83.8% [77.1%-88.8%] | 97.9% [95.3%-99.1%] | 25 | 5 |
| colistin | AMRFinderPlus | 163 (77/86) | 24.7% [16.4%-35.4%] | 95.3% [88.6%-98.2%] | 58 | 4 |
| colistin | AMRFinderPlus+regulators | 163 (77/86) | 24.7% [16.4%-35.4%] | 95.3% [88.6%-98.2%] | 58 | 4 |
| colistin | Kleborate | 163 (77/86) | 44.2% [33.6%-55.3%] | 95.3% [88.6%-98.2%] | 43 | 4 |
| colistin_all | AMRFinderPlus | 214 (104/110) | 24.0% [16.8%-33.1%] | 96.4% [91.0%-98.6%] | 79 | 4 |
| colistin_all | AMRFinderPlus+regulators | 214 (104/110) | 24.0% [16.8%-33.1%] | 96.4% [91.0%-98.6%] | 79 | 4 |
| colistin_all | Kleborate | 214 (104/110) | 34.6% [26.2%-44.2%] | 96.4% [91.0%-98.6%] | 68 | 4 |
| tigecycline | AMRFinderPlus | 200 (100/100) | 2.0% [0.6%-7.0%] | 100.0% [96.3%-100.0%] | 98 | 0 |
| tigecycline | AMRFinderPlus+regulators | 200 (100/100) | 19.0% [12.5%-27.8%] | 85.0% [76.7%-90.7%] | 81 | 15 |
| tigecycline | Kleborate | 200 (100/100) | 0.0% [0.0%-3.7%] | 100.0% [96.3%-100.0%] | 100 | 0 |

## AMRFinderPlus (core) vs Kleborate prediction agreement

| Drug | Agreement |
|---|---|
| meropenem | 99.8% |
| ciprofloxacin | 99.7% |
| gentamicin | 94.7% |
| colistin | 84.7% |
| colistin_all | 86.4% |
| tigecycline | 99.0% |

Discordant genome-drug pairs (core AMRFinderPlus or Kleborate wrong): 275, in `reports/pilot/discordant.tsv` with each tool's evidence for review.
