# Paired comparisons between predictors

- Generated (UTC): 2026-10-25T05:47:12+00:00
- 2,000 paired bootstrap resamples of held-out genomes, seed 20261024
- Difference is balanced accuracy of A minus B; interval is the 2.5–97.5 percentile
- Splits with fewer than 10 genomes of either class are omitted
- Family of 54 comparisons; q is the Benjamini-Hochberg adjusted p-value
- Bonferroni threshold for this family is p < 0.00093
- p is floored at 1/2,000 = 0.00050, the resolution of the bootstrap

Uncorrected, 33 of 54 intervals exclude zero; 31 survive Benjamini-Hochberg at 5% and 24 survive Bonferroni.

| Drug | Split | Comparison | Difference | 95% interval | p | BH q | Survives BH |
|---|---|---|---|---|---|---|---|
| ciprofloxacin | lineage | rule vs logreg | -0.015 | -0.023 to -0.008 | <0.0005 | 0.0011 | Yes |
| ciprofloxacin | lineage | rule vs gboost | -0.015 | -0.023 to -0.008 | <0.0005 | 0.0011 | Yes |
| ciprofloxacin | lineage | logreg vs gboost | +0.000 | -0.005 to +0.005 | 0.9810 | 0.9890 | No |
| ciprofloxacin | random | rule vs logreg | -0.017 | -0.025 to -0.010 | <0.0005 | 0.0011 | Yes |
| ciprofloxacin | random | rule vs gboost | -0.020 | -0.027 to -0.013 | <0.0005 | 0.0011 | Yes |
| ciprofloxacin | random | logreg vs gboost | -0.002 | -0.007 to +0.002 | 0.3540 | 0.4248 | No |
| ciprofloxacin | study | rule vs logreg | -0.008 | -0.016 to -0.001 | 0.0300 | 0.0506 | No |
| ciprofloxacin | study | rule vs gboost | -0.013 | -0.020 to -0.006 | <0.0005 | 0.0011 | Yes |
| ciprofloxacin | study | logreg vs gboost | -0.005 | -0.010 to -0.000 | 0.0380 | 0.0622 | No |
| ciprofloxacin | temporal | rule vs logreg | -0.038 | -0.064 to -0.013 | 0.0020 | 0.0042 | Yes |
| ciprofloxacin | temporal | rule vs gboost | -0.031 | -0.051 to -0.010 | <0.0005 | 0.0011 | Yes |
| ciprofloxacin | temporal | logreg vs gboost | +0.007 | -0.010 to +0.026 | 0.4230 | 0.4966 | No |
| colistin | lineage | rule vs logreg | +0.028 | +0.003 to +0.052 | 0.0260 | 0.0453 | Yes |
| colistin | lineage | rule vs gboost | +0.050 | +0.026 to +0.075 | <0.0005 | 0.0011 | Yes |
| colistin | lineage | logreg vs gboost | +0.022 | +0.005 to +0.041 | 0.0120 | 0.0231 | Yes |
| colistin | random | rule vs logreg | -0.085 | -0.113 to -0.059 | <0.0005 | 0.0011 | Yes |
| colistin | random | rule vs gboost | -0.083 | -0.111 to -0.054 | <0.0005 | 0.0011 | Yes |
| colistin | random | logreg vs gboost | +0.002 | -0.014 to +0.018 | 0.7940 | 0.8407 | No |
| colistin | study | rule vs logreg | -0.063 | -0.088 to -0.037 | <0.0005 | 0.0011 | Yes |
| colistin | study | rule vs gboost | -0.051 | -0.078 to -0.026 | <0.0005 | 0.0011 | Yes |
| colistin | study | logreg vs gboost | +0.012 | -0.007 to +0.029 | 0.2420 | 0.3267 | No |
| gentamicin | lineage | rule vs logreg | -0.010 | -0.018 to -0.002 | 0.0070 | 0.0140 | Yes |
| gentamicin | lineage | rule vs gboost | -0.004 | -0.012 to +0.003 | 0.2740 | 0.3536 | No |
| gentamicin | lineage | logreg vs gboost | +0.005 | -0.001 to +0.012 | 0.0960 | 0.1440 | No |
| gentamicin | random | rule vs logreg | -0.018 | -0.025 to -0.010 | <0.0005 | 0.0011 | Yes |
| gentamicin | random | rule vs gboost | -0.019 | -0.026 to -0.011 | <0.0005 | 0.0011 | Yes |
| gentamicin | random | logreg vs gboost | -0.001 | -0.007 to +0.005 | 0.7630 | 0.8240 | No |
| gentamicin | study | rule vs logreg | +0.024 | +0.019 to +0.029 | <0.0005 | 0.0011 | Yes |
| gentamicin | study | rule vs gboost | +0.021 | +0.016 to +0.026 | <0.0005 | 0.0011 | Yes |
| gentamicin | study | logreg vs gboost | -0.003 | -0.007 to +0.000 | 0.0700 | 0.1112 | No |
| gentamicin | temporal | rule vs logreg | +0.005 | -0.017 to +0.026 | 0.6580 | 0.7560 | No |
| gentamicin | temporal | rule vs gboost | +0.004 | -0.019 to +0.027 | 0.7040 | 0.7920 | No |
| gentamicin | temporal | logreg vs gboost | -0.000 | -0.017 to +0.016 | 0.9890 | 0.9890 | No |
| meropenem | lineage | rule vs logreg | -0.006 | -0.011 to -0.001 | 0.0180 | 0.0335 | Yes |
| meropenem | lineage | rule vs gboost | -0.004 | -0.009 to +0.001 | 0.1270 | 0.1805 | No |
| meropenem | lineage | logreg vs gboost | +0.002 | -0.002 to +0.006 | 0.2860 | 0.3592 | No |
| meropenem | random | rule vs logreg | -0.012 | -0.016 to -0.007 | <0.0005 | 0.0011 | Yes |
| meropenem | random | rule vs gboost | -0.011 | -0.017 to -0.006 | <0.0005 | 0.0011 | Yes |
| meropenem | random | logreg vs gboost | +0.000 | -0.003 to +0.003 | 0.9040 | 0.9388 | No |
| meropenem | study | rule vs logreg | -0.003 | -0.007 to +0.002 | 0.2750 | 0.3536 | No |
| meropenem | study | rule vs gboost | +0.004 | -0.001 to +0.009 | 0.1040 | 0.1518 | No |
| meropenem | study | logreg vs gboost | +0.007 | +0.003 to +0.011 | 0.0010 | 0.0022 | Yes |
| meropenem | temporal | rule vs logreg | +0.016 | +0.000 to +0.038 | 0.0900 | 0.1389 | No |
| meropenem | temporal | rule vs gboost | +0.004 | +0.001 to +0.009 | 0.0250 | 0.0450 | Yes |
| meropenem | temporal | logreg vs gboost | -0.012 | -0.033 to +0.004 | 0.1880 | 0.2603 | No |
| tigecycline | lineage | rule vs logreg | -0.203 | -0.236 to -0.171 | <0.0005 | 0.0011 | Yes |
| tigecycline | lineage | rule vs gboost | -0.107 | -0.144 to -0.073 | <0.0005 | 0.0011 | Yes |
| tigecycline | lineage | logreg vs gboost | +0.096 | +0.062 to +0.129 | <0.0005 | 0.0011 | Yes |
| tigecycline | random | rule vs logreg | -0.255 | -0.284 to -0.226 | <0.0005 | 0.0011 | Yes |
| tigecycline | random | rule vs gboost | -0.244 | -0.278 to -0.209 | <0.0005 | 0.0011 | Yes |
| tigecycline | random | logreg vs gboost | +0.012 | -0.011 to +0.034 | 0.3070 | 0.3768 | No |
| tigecycline | study | rule vs logreg | -0.107 | -0.141 to -0.071 | <0.0005 | 0.0011 | Yes |
| tigecycline | study | rule vs gboost | -0.112 | -0.147 to -0.075 | <0.0005 | 0.0011 | Yes |
| tigecycline | study | logreg vs gboost | -0.006 | -0.037 to +0.026 | 0.7190 | 0.7924 | No |
