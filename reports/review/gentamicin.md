# Gentamicin discordance review

- Generated (UTC): 2026-10-24T14:24:13+00:00
- Genomes with a gentamicin S/R label: 4,760 (R 2,113, S 2,647)
- Rule: AMRFinderPlus core (Subclass GENTAMICIN or a 16S rRNA methyltransferase)
- EUCAST v16.1 gentamicin (systemic): S <= 2, R > 2 mg/L

| Outcome | Genomes |
|---|---|
| TN | 2,512 |
| TP | 1,799 |
| FN | 314 |
| FP | 135 |

## 1. MIC distribution by outcome (mg/L)

299 of 314 false negatives have an MIC; 213 of those (71%) are at 4 mg/L, one dilution above the breakpoint.

| Outcome | <=0.5 | 1 | 2 | 4 | 8 | 16 | >16 | no MIC |
|---|---|---|---|---|---|---|---|---|
| FN | 0 | 0 | 0 | 213 | 63 | 15 | 8 | 15 |
| FP | 2 | 80 | 30 | 0 | 0 | 0 | 0 | 23 |
| TN | 171 | 494 | 688 | 0 | 0 | 0 | 0 | 1,159 |
| TP | 0 | 0 | 0 | 86 | 796 | 350 | 117 | 450 |

## 2. Aminoglycoside genes: false negatives vs susceptible non-carriers

283 of 314 false negatives carry at least one acquired aminoglycoside gene that the gentamicin rule does not count. FN = rule false negatives (R, no gentamicin marker); TN = susceptible genomes without a gentamicin marker.

| Gene | in_FN | in_TN | FN_share | TN_share |
|---|---|---|---|---|
| aadA2 | 227 | 568 | 72.3% | 22.6% |
| aac(6')-Ib | 199 | 312 | 63.4% | 12.4% |
| aph(3')-Ia | 191 | 357 | 60.8% | 14.2% |
| aph(3'')-Ib | 31 | 483 | 9.9% | 19.2% |
| aph(6)-Id | 30 | 479 | 9.6% | 19.1% |
| aac(6')-Ib-cr5 | 17 | 486 | 5.4% | 19.3% |
| aadA1 | 10 | 96 | 3.2% | 3.8% |
| aadA16 | 8 | 133 | 2.5% | 5.3% |
| aac(6')-Ib-cr | 5 | 5 | 1.6% | 0.2% |
| aph(3')-VI | 3 | 21 | 1.0% | 0.8% |
| aadA22 | 1 | 0 | 0.3% | 0.0% |
| aadA3 | 1 | 0 | 0.3% | 0.0% |
| aadA8 | 1 | 2 | 0.3% | 0.1% |
| aac(6')-Ia | 0 | 1 | 0.0% | 0.0% |
| aadA11 | 0 | 1 | 0.0% | 0.0% |

## 3. Study and method effects

BioProjects with the most false negatives:

| BioProject | FN | all_R | FN_share_of_project_R |
|---|---|---|---|
| PRJNA376414 | 207 | 841 | 24.6% |
| PRJNA504600 | 24 | 24 | 100.0% |
| PRJNA385863 | 23 | 24 | 95.8% |
| PRJEB5065 | 8 | 131 | 6.1% |
| PRJNA396774 | 5 | 48 | 10.4% |
| PRJEB28400 | 5 | 168 | 3.0% |
| PRJNA288601 | 5 | 53 | 9.4% |
| PRJEB22890 | 5 | 63 | 7.9% |

Testing method of the highest MIC record:

| Method | TP | FN | FP | TN |
|---|---|---|---|---|
| Agar dilution | 123 | 8 | 3 | 111 |
| Broth dilution | 981 | 282 | 33 | 1,001 |
| Disk diffusion | 0 | 1 | 0 | 15 |
| MIC | 245 | 8 | 76 | 226 |
| no MIC | 450 | 15 | 23 | 1,159 |

## 4. False positives (gentamicin marker present, phenotype S)

135 genomes. Most frequent markers:

| Marker | FP genomes |
|---|---|
| rmtB1 | 65 |
| aac(3)-IId | 14 |
| aac(3)-IIe | 14 |
| aac(6')-Ib | 13 |
| armA | 9 |
| ant(2'')-Ia | 4 |
| aac(6')-Ib4 | 4 |
| aac(3)-IVa | 3 |
| aac(6')-Ib' | 3 |
| aac(3)-Ia | 3 |

Per-genome detail of false negatives and false positives: `reports/review/gentamicin_genomes.tsv`.
