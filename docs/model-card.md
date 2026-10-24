# Model card: K. pneumoniae AMR genotype-phenotype benchmark

Date: 2026-10-24. Code: `src/kamr/benchmark/`. Results: `reports/benchmark/results.md`.

## What was evaluated
Three predictors of laboratory resistance (R vs S) for five antibiotics, from genome assemblies:

| Predictor | Description |
|---|---|
| Rule | AMRFinderPlus 4.2.7 core rule per drug (`src/kamr/analysis/rules.py`); for colistin, plus disrupted or missense mgrB (Sprint 2) |
| Logistic regression | L2, C = 1, class-balanced, on 226 binary genomic features |
| Gradient boosting | scikit-learn HistGradientBoostingClassifier, 200 iterations, 15 leaves, class-balanced |

Features: presence of AMRFinderPlus gene families (OXA-48-like kept apart from other OXA), 100 point
mutations, 33 AMRFinderPlus subclasses and 2 mgrB features (`reports/benchmark/feature_summary.md`).

## Data
- 6,652 public K. pneumoniae assemblies from BV-BRC that pass assembly QC, with laboratory AST
  (snapshot 2026-10-22/23, manifests in `snapshots/`).
- Labels: laboratory MICs re-interpreted under EUCAST v16.1 (ADR-002); I and conflicting labels
  excluded; colistin restricted to broth microdilution (ADR-003); one study's gentamicin labels
  excluded as implausible (ADR-004, `config/label_flags.csv`).
- Labelled genomes per drug (R/S): meropenem 1,854/3,687; ciprofloxacin 3,441/1,199;
  gentamicin 1,945/2,483 after flags; colistin 367/1,161; tigecycline 261/687.

## Evaluation
- Splits: random 5-fold; study-aware 5-fold (grouped by BioProject); lineage-aware 5-fold (grouped
  by sequence type); temporal (train up to 2014, test after; tigecycline 2015).
- Metrics: sensitivity, specificity (95% Wilson intervals), balanced accuracy, AUROC, all on
  held-out genomes; reported on all labels and without borderline-MIC labels (ADR-004).
- A split is scored only if its test set has at least 10 R and 10 S genomes.

## Results in brief
- Meropenem, ciprofloxacin, gentamicin: rules reach 0.91-0.94 balanced accuracy; models add 1-2
  points under random splits and little or nothing under grouped splits. Without borderline
  labels the gentamicin rule (0.959) outperforms both models.
- Colistin: random-split model accuracy (0.84) falls to 0.71-0.73 under lineage-aware splits,
  below the rule with the mgrB check (0.81).
- Tigecycline: models appear to beat the rule (0.77 vs 0.51) under random splits, but reach only
  0.62 under study-aware splits.
- Random-split optimism (random minus the lower grouped split): 0.7-4 points for well-explained
  drugs, 11-15 points for colistin and tigecycline.

## Intended use
- A reproducible reference for how much validation design changes AMR prediction estimates on
  public data, and for where known mechanisms do and do not explain resistance.
- Teaching and method comparison.

## Not intended for
- Clinical decisions or susceptibility reporting. These are research predictions from public
  surveillance data and are not validated diagnostics.
- Estimating resistance prevalence: the cohort reflects what studies chose to sequence and test.

## Limitations
- Public data are biased towards outbreak, carbapenem-resistant and high-risk-clone isolates;
  two BioProjects hold 44% of genomes.
- Labels come from many laboratories and methods; borderline MICs and study-level errors remain
  (ADR-004). Colistin and tigecycline have few resistant isolates.
- 35% of genomes lack a collection year; the temporal split is not evaluable for colistin
  (1 R in test) or tigecycline (0 R).
- Features are presence/absence from assemblies: gene copy number, expression, promoter changes,
  IS insertions upstream of mgrB and plasmid context are not modelled.
- Hyperparameters were fixed, not tuned; results show the effect of validation design, not the
  best achievable model.

## How to reproduce
See "Reproduce the results" in the README. All results are deterministic for a fixed seed
(20261024) and were identical across two runs.
