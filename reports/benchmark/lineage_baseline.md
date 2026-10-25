# How much of model performance is lineage?

- Generated (UTC): 2026-10-25T05:53:10+00:00
- ST overlap: share of test genomes whose sequence type also appears in training
- Lineage-only: predicts the majority label of the genome's sequence type in training,
  falling back to the overall training majority for unseen types. Uses no resistance
  determinant of any kind.
- Best model: the better of logistic regression and gradient boosting from the benchmark

| Drug | Split | ST overlap | BioProject overlap | Lineage-only | Best model |
|---|---|---|---|---|---|
| meropenem | random | 91% | 97% | 0.737 | 0.951 |
| meropenem | study | 88% | 0% | 0.664 | 0.942 |
| meropenem | lineage | 0% | 95% | 0.500 | 0.945 |
| meropenem | temporal | 78% | 77% | 0.623 | 0.941 |
| ciprofloxacin | random | 89% | 95% | 0.728 | 0.954 |
| ciprofloxacin | study | 84% | 0% | 0.620 | 0.948 |
| ciprofloxacin | lineage | 0% | 92% | 0.500 | 0.950 |
| ciprofloxacin | temporal | 74% | 82% | 0.670 | 0.949 |
| gentamicin | random | 88% | 97% | 0.713 | 0.925 |
| gentamicin | study | 84% | 0% | 0.650 | 0.886 |
| gentamicin | lineage | 0% | 95% | 0.500 | 0.916 |
| gentamicin | temporal | 73% | 76% | 0.633 | 0.912 |
| colistin | random | 91% | 93% | 0.554 | 0.844 |
| colistin | study | 83% | 0% | 0.559 | 0.822 |
| colistin | lineage | 0% | 87% | 0.500 | 0.731 |
| tigecycline | random | 84% | 100% | 0.654 | 0.769 |
| tigecycline | study | 74% | 0% | 0.453 | 0.626 |
| tigecycline | lineage | 0% | 93% | 0.500 | 0.716 |

A lineage-only predictor that scores well under random splitting and near 0.5 under
lineage-aware splitting demonstrates that lineage alone carries predictive signal, and
that conventional validation lets models use it.
