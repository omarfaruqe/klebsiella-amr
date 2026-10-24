# Meropenem discordance review

- Generated (UTC): 2026-10-24T13:09:03+00:00
- Genomes with a meropenem S/R label: 5,541 (R 1,854, S 3,687)
- Carbapenemase calls: AMRFinderPlus acquired genes with Subclass CARBAPENEM
- EUCAST v16.1 meropenem: S <= 2 mg/L, R > 8 mg/L

Carbapenemase present: 2,118 genomes, of which 326 test S (rule false positives).
No carbapenemase: 3,423 genomes, of which 62 test R (rule false negatives).

## 1. Carbapenemase family vs phenotype

| Family | genomes | R | S | S_share |
|---|---|---|---|---|
| KPC | 1,278 | 1,203 | 75 | 5.9% |
| OXA-48-like | 344 | 162 | 182 | 52.9% |
| NDM | 280 | 252 | 28 | 10.0% |
| multiple | 149 | 148 | 1 | 0.7% |
| VIM | 44 | 11 | 33 | 75.0% |
| IMP | 22 | 16 | 6 | 27.3% |
| other | 1 | 0 | 1 | 100.0% |

## 2. MIC of susceptible carbapenemase carriers (mg/L)

240 of 326 false positives have an MIC; the rest have only a reported S call.
- 191 of 240 (80%) sit at 1-2 mg/L, within two doubling dilutions of the S breakpoint.
- 224 of 240 (93%) are above the EUCAST carbapenemase screening cut-off (meropenem MIC > 0.125 mg/L): consistent with a carbapenemase even though clinically S.

| Family | <=0.125 | 0.25 | 0.5 | 1 | 2 | 4-8 |
|---|---|---|---|---|---|---|
| IMP | 0 | 0 | 1 | 3 | 2 | 0 |
| KPC | 7 | 1 | 2 | 33 | 28 | 1 |
| NDM | 1 | 1 | 0 | 1 | 3 | 0 |
| OXA-48-like | 6 | 3 | 15 | 48 | 49 | 0 |
| VIM | 1 | 4 | 5 | 12 | 11 | 0 |
| multiple | 1 | 0 | 0 | 0 | 0 | 0 |
| other | 0 | 0 | 0 | 1 | 0 | 0 |

## 3. Mechanism routes without a carbapenemase

Porin defects = AMRFinderPlus point/disruption calls in ompK35 or ompK36.

| Route | R, no carbapenemase | S, no carbapenemase |
|---|---|---|
| ESBL/AmpC + porin defect | 29 | 624 |
| ESBL/AmpC only | 22 | 1,382 |
| porin defect only | 8 | 162 |
| neither | 3 | 1,193 |

## 4. Where the false positives come from

Top BioProjects:

| BioProject | genomes |
|---|---|
| PRJEB10018 | 158 |
| PRJNA376414 | 56 |
| PRJEB24082 | 46 |
| PRJEB19229 | 19 |
| missing | 15 |
| PRJDB10842 | 5 |
| PRJEB28400 | 5 |
| PRJNA288601 | 3 |

Top STs:

| ST | genomes |
|---|---|
| ST15 | 60 |
| ST258 | 27 |
| ST11 | 26 |
| ST307 | 24 |
| ST133 | 12 |
| ST405 | 10 |
| ST16 | 9 |
| ST37 | 9 |

Per-genome detail: `reports/review/meropenem_genomes.tsv`.
