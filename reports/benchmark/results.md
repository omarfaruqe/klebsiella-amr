# Leakage-aware benchmark

- Generated (UTC): 2026-10-24T15:56:33+00:00
- Features: 226 (reports/benchmark/feature_summary.md); flagged labels excluded
- Predictors: rule (AMRFinderPlus core; colistin + mgrB), logreg, gboost
- Splits: random 5-fold; study = 5-fold grouped by BioProject; lineage = 5-fold grouped by ST; temporal = train on earlier, test on the latest ~20% of dated genomes
- Temporal splits: meropenem train <= 2014, test > 2014; ciprofloxacin train <= 2014, test > 2014; gentamicin train <= 2014, test > 2014; colistin train <= 2014, test > 2014; not evaluable (1 R / 33 S in test); tigecycline train <= 2015, test > 2015; not evaluable (0 R / 28 S in test)
- Splits whose test set has fewer than 10 R or 10 S genomes are not evaluable and shown as n/a.

## Balanced accuracy, all test labels

| Drug | Model | random | study | lineage | temporal |
|---|---|---|---|---|---|
| ciprofloxacin | gboost | 0.954 | 0.948 | 0.950 | 0.941 |
| ciprofloxacin | logreg | 0.952 | 0.943 | 0.950 | 0.949 |
| ciprofloxacin | rule | 0.935 | 0.935 | 0.935 | 0.911 |
| colistin | gboost | 0.842 | 0.810 | 0.709 | n/a |
| colistin | logreg | 0.844 | 0.822 | 0.731 | n/a |
| colistin | rule | 0.759 | 0.759 | 0.759 | n/a |
| gentamicin | gboost | 0.925 | 0.886 | 0.911 | 0.912 |
| gentamicin | logreg | 0.924 | 0.882 | 0.916 | 0.911 |
| gentamicin | rule | 0.906 | 0.906 | 0.906 | 0.916 |
| meropenem | gboost | 0.950 | 0.935 | 0.943 | 0.941 |
| meropenem | logreg | 0.951 | 0.942 | 0.945 | 0.929 |
| meropenem | rule | 0.939 | 0.939 | 0.939 | 0.946 |
| tigecycline | gboost | 0.757 | 0.626 | 0.621 | n/a |
| tigecycline | logreg | 0.769 | 0.620 | 0.716 | n/a |
| tigecycline | rule | 0.513 | 0.513 | 0.513 | n/a |

## AUROC, all test labels (models only)

| Drug | Model | random | study | lineage | temporal |
|---|---|---|---|---|---|
| ciprofloxacin | gboost | 0.988 | 0.983 | 0.987 | 0.975 |
| ciprofloxacin | logreg | 0.987 | 0.983 | 0.986 | 0.982 |
| colistin | gboost | 0.892 | 0.867 | 0.811 | n/a |
| colistin | logreg | 0.906 | 0.879 | 0.845 | n/a |
| gentamicin | gboost | 0.975 | 0.949 | 0.965 | 0.977 |
| gentamicin | logreg | 0.973 | 0.955 | 0.967 | 0.973 |
| meropenem | gboost | 0.981 | 0.973 | 0.976 | 0.957 |
| meropenem | logreg | 0.981 | 0.975 | 0.977 | 0.965 |
| tigecycline | gboost | 0.822 | 0.636 | 0.697 | n/a |
| tigecycline | logreg | 0.823 | 0.639 | 0.736 | n/a |

## Balanced accuracy without borderline-MIC labels (ADR-004)

| Drug | Model | random | study | lineage | temporal |
|---|---|---|---|---|---|
| ciprofloxacin | gboost | 0.958 | 0.954 | 0.955 | 0.946 |
| ciprofloxacin | logreg | 0.959 | 0.951 | 0.958 | 0.956 |
| ciprofloxacin | rule | 0.936 | 0.936 | 0.936 | 0.910 |
| colistin | gboost | 0.844 | 0.823 | 0.713 | n/a |
| colistin | logreg | 0.846 | 0.834 | 0.738 | n/a |
| colistin | rule | 0.763 | 0.763 | 0.763 | n/a |
| gentamicin | gboost | 0.948 | 0.939 | 0.938 | 0.943 |
| gentamicin | logreg | 0.947 | 0.935 | 0.943 | 0.940 |
| gentamicin | rule | 0.959 | 0.959 | 0.959 | 0.960 |
| meropenem | gboost | 0.957 | 0.944 | 0.952 | 0.947 |
| meropenem | logreg | 0.958 | 0.951 | 0.953 | 0.934 |
| meropenem | rule | 0.951 | 0.951 | 0.951 | 0.950 |
| tigecycline | gboost | 0.757 | 0.626 | 0.621 | n/a |
| tigecycline | logreg | 0.769 | 0.620 | 0.716 | n/a |
| tigecycline | rule | 0.513 | 0.513 | 0.513 | n/a |

## Optimism of random splits

Gap = balanced accuracy under random split minus the lower of study and lineage splits.

| Drug | Model | Gap |
|---|---|---|
| ciprofloxacin | gboost | 0.007 |
| ciprofloxacin | logreg | 0.009 |
| ciprofloxacin | rule | 0.000 |
| colistin | gboost | 0.133 |
| colistin | logreg | 0.113 |
| colistin | rule | 0.000 |
| gentamicin | gboost | 0.039 |
| gentamicin | logreg | 0.042 |
| gentamicin | rule | 0.000 |
| meropenem | gboost | 0.016 |
| meropenem | logreg | 0.009 |
| meropenem | rule | 0.000 |
| tigecycline | gboost | 0.136 |
| tigecycline | logreg | 0.149 |
| tigecycline | rule | 0.000 |

Full metrics with sensitivity, specificity and 95% intervals: `reports/benchmark/results.csv`.
