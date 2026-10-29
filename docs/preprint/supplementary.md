# Supplementary material

**Why genomic prediction of colistin resistance does not transfer between *Klebsiella pneumoniae* collections**

Md Omar Faruqe and Md Humayun Kabir

---

## Supplementary tables

**Supplementary Table S1.** Colistin sensitivity and specificity for each rule, with and without
the targeted *mgrB* disruption check, on broth-microdilution labels (n = 1,528; 367 resistant).
Plotted in Figure 1.

| Rule | Sensitivity (95% CI) | Specificity (95% CI) |
|---|---|---|
| AMRFinderPlus | 28.3% (24.0–33.2) | 97.8% (96.7–98.5) |
| **AMRFinderPlus + *mgrB* disrupted** | **56.7% (51.6–61.7)** | **95.2% (93.8–96.3)** |
| Kleborate | 39.5% (34.6–44.6) | 96.3% (95.0–97.2) |
| **Kleborate + *mgrB* disrupted** | **45.5% (40.5–50.6)** | **95.3% (94.0–96.4)** |
| *mgrB* disrupted alone | 41.7% (36.8–46.8) | 96.5% (95.2–97.4) |

**Supplementary Table S2.** Sensitivity of the *mgrB* disruption rule in each held-out study
(leave-one-study-out). Plotted in Figure 2.

| Held-out study | Resistant | Sensitivity, disruption only |
|---|---|---|
| PRJEB10018 | 162 | 39.5% |
| PRJNA504600 | 27 | 48.1% |
| PRJEB22890 | 63 | 50.8% |
| PRJNA288601 | 20 | 55.0% |
| PRJNA292904 | 7 | 57.1% |
| PRJNA307517 | 19 | 63.2% |
| PRJEB10561 | 17 | 70.6% |
| **Pooled** | **315** | **47.0% (41.5–52.5)** |

**Supplementary Table S3.** Balanced accuracy of the per-drug rule, logistic regression and
gradient boosting under each split scheme. Plotted in Figure 3.

| Split | Rule + *mgrB* | Logistic regression | Gradient boosting |
|---|---|---|---|
| Random | 0.759 | **0.844** | **0.842** |
| Study-aware | 0.759 | **0.822** | **0.810** |
| Lineage-aware | **0.759** | 0.731 | 0.709 |
