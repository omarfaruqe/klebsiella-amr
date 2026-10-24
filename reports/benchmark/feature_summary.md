# Benchmark inputs

- Generated (UTC): 2026-10-24T14:37:15+00:00
- Genomes: 6,652; features: 226 (present in >= 5 genomes; 100 point mutations, 33 AMRFinderPlus subclasses, 2 mgrB, rest gene families)
- Colistin labels: broth microdilution only (ADR-003)
- Borderline: highest MIC exactly at s_le or at 2 x r_gt (ADR-004)
- Flagged: study-level exclusions from config/label_flags.csv (ADR-004)

## Labels per drug

| Drug | Labelled | R | S | Borderline | Flagged | Primary set (not flagged) |
|---|---|---|---|---|---|---|
| meropenem | 5,541 | 1,854 | 3,687 | 269 | 0 | 5,541 |
| ciprofloxacin | 4,640 | 3,441 | 1,199 | 184 | 0 | 4,640 |
| gentamicin | 4,760 | 2,113 | 2,647 | 295 | 332 | 4,428 |
| colistin | 1,528 | 367 | 1,161 | 55 | 0 | 1,528 |
| tigecycline | 948 | 261 | 687 | 0 | 0 | 948 |

## Split metadata coverage (labelled genomes)

- Labelled genomes: 6,469
- BioProject missing: 322
- ST missing: 98
- Collection year missing: 2,245

## 20 most common features

| Feature | Genomes | Share |
|---|---|---|
| emrD | 6,643 | 99.9% |
| subclass:EFFLUX | 6,643 | 99.9% |
| subclass:FOSFOMYCIN | 6,641 | 99.8% |
| fosA | 6,641 | 99.8% |
| blaSHV | 6,459 | 97.1% |
| subclass:BETA-LACTAM | 6,369 | 95.7% |
| subclass:TRIMETHOPRIM | 4,473 | 67.2% |
| dfrA | 4,465 | 67.1% |
| subclass:SULFONAMIDE | 4,276 | 64.3% |
| sul | 4,276 | 64.3% |
| subclass:STREPTOMYCIN | 4,187 | 62.9% |
| subclass:CEPHALOSPORIN | 4,068 | 61.2% |
| subclass:KANAMYCIN | 3,684 | 55.4% |
| subclass:CHLORAMPHENICOL | 3,537 | 53.2% |
| parC_S80I | 3,488 | 52.4% |
| subclass:TOBRAMYCIN | 3,416 | 51.4% |
| blaTEM | 3,406 | 51.2% |
| subclass:AMIKACIN | 3,386 | 50.9% |
| blaCTX-M | 3,329 | 50.0% |
| blaOXA | 3,118 | 46.9% |
