# Paired comparisons between predictors

- Generated (UTC): 2026-10-24T16:04:09+00:00
- 2,000 paired bootstrap resamples of held-out genomes, seed 20261024
- Difference is balanced accuracy of A minus B; interval is the 2.5–97.5 percentile
- 'Yes' means the 95% interval excludes zero
- Splits with fewer than 10 genomes of either class are omitted

| Drug | Split | Comparison | Difference | 95% interval | Distinguishable |
|---|---|---|---|---|---|
| ciprofloxacin | lineage | rule vs logreg | -0.015 | -0.023 to -0.008 | Yes |
| ciprofloxacin | lineage | rule vs gboost | -0.015 | -0.023 to -0.008 | Yes |
| ciprofloxacin | lineage | logreg vs gboost | +0.000 | -0.005 to +0.005 | No |
| ciprofloxacin | random | rule vs logreg | -0.017 | -0.025 to -0.010 | Yes |
| ciprofloxacin | random | rule vs gboost | -0.020 | -0.027 to -0.013 | Yes |
| ciprofloxacin | random | logreg vs gboost | -0.002 | -0.007 to +0.002 | No |
| ciprofloxacin | study | rule vs logreg | -0.008 | -0.016 to -0.001 | Yes |
| ciprofloxacin | study | rule vs gboost | -0.013 | -0.020 to -0.006 | Yes |
| ciprofloxacin | study | logreg vs gboost | -0.005 | -0.010 to -0.000 | Yes |
| ciprofloxacin | temporal | rule vs logreg | -0.038 | -0.064 to -0.013 | Yes |
| ciprofloxacin | temporal | rule vs gboost | -0.031 | -0.051 to -0.010 | Yes |
| ciprofloxacin | temporal | logreg vs gboost | +0.007 | -0.010 to +0.026 | No |
| colistin | lineage | rule vs logreg | +0.028 | +0.003 to +0.052 | Yes |
| colistin | lineage | rule vs gboost | +0.050 | +0.026 to +0.075 | Yes |
| colistin | lineage | logreg vs gboost | +0.022 | +0.005 to +0.041 | Yes |
| colistin | random | rule vs logreg | -0.085 | -0.113 to -0.059 | Yes |
| colistin | random | rule vs gboost | -0.083 | -0.111 to -0.054 | Yes |
| colistin | random | logreg vs gboost | +0.002 | -0.014 to +0.018 | No |
| colistin | study | rule vs logreg | -0.063 | -0.088 to -0.037 | Yes |
| colistin | study | rule vs gboost | -0.051 | -0.078 to -0.026 | Yes |
| colistin | study | logreg vs gboost | +0.012 | -0.007 to +0.029 | No |
| gentamicin | lineage | rule vs logreg | -0.010 | -0.018 to -0.002 | Yes |
| gentamicin | lineage | rule vs gboost | -0.004 | -0.012 to +0.003 | No |
| gentamicin | lineage | logreg vs gboost | +0.005 | -0.001 to +0.012 | No |
| gentamicin | random | rule vs logreg | -0.018 | -0.025 to -0.010 | Yes |
| gentamicin | random | rule vs gboost | -0.019 | -0.026 to -0.011 | Yes |
| gentamicin | random | logreg vs gboost | -0.001 | -0.007 to +0.005 | No |
| gentamicin | study | rule vs logreg | +0.024 | +0.019 to +0.029 | Yes |
| gentamicin | study | rule vs gboost | +0.021 | +0.016 to +0.026 | Yes |
| gentamicin | study | logreg vs gboost | -0.003 | -0.007 to +0.000 | No |
| gentamicin | temporal | rule vs logreg | +0.005 | -0.017 to +0.026 | No |
| gentamicin | temporal | rule vs gboost | +0.004 | -0.019 to +0.027 | No |
| gentamicin | temporal | logreg vs gboost | -0.000 | -0.017 to +0.016 | No |
| meropenem | lineage | rule vs logreg | -0.006 | -0.011 to -0.001 | Yes |
| meropenem | lineage | rule vs gboost | -0.004 | -0.009 to +0.001 | No |
| meropenem | lineage | logreg vs gboost | +0.002 | -0.002 to +0.006 | No |
| meropenem | random | rule vs logreg | -0.012 | -0.016 to -0.007 | Yes |
| meropenem | random | rule vs gboost | -0.011 | -0.017 to -0.006 | Yes |
| meropenem | random | logreg vs gboost | +0.000 | -0.003 to +0.003 | No |
| meropenem | study | rule vs logreg | -0.003 | -0.007 to +0.002 | No |
| meropenem | study | rule vs gboost | +0.004 | -0.001 to +0.009 | No |
| meropenem | study | logreg vs gboost | +0.007 | +0.003 to +0.011 | Yes |
| meropenem | temporal | rule vs logreg | +0.016 | +0.000 to +0.038 | No |
| meropenem | temporal | rule vs gboost | +0.004 | +0.001 to +0.009 | Yes |
| meropenem | temporal | logreg vs gboost | -0.012 | -0.033 to +0.004 | No |
| tigecycline | lineage | rule vs logreg | -0.203 | -0.236 to -0.171 | Yes |
| tigecycline | lineage | rule vs gboost | -0.107 | -0.144 to -0.073 | Yes |
| tigecycline | lineage | logreg vs gboost | +0.096 | +0.062 to +0.129 | Yes |
| tigecycline | random | rule vs logreg | -0.255 | -0.284 to -0.226 | Yes |
| tigecycline | random | rule vs gboost | -0.244 | -0.278 to -0.209 | Yes |
| tigecycline | random | logreg vs gboost | +0.012 | -0.011 to +0.034 | No |
| tigecycline | study | rule vs logreg | -0.107 | -0.141 to -0.071 | Yes |
| tigecycline | study | rule vs gboost | -0.112 | -0.147 to -0.075 | Yes |
| tigecycline | study | logreg vs gboost | -0.006 | -0.037 to +0.026 | No |
