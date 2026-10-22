# Feasibility spike log

Goal: decide which drugs have enough paired genome + laboratory AST data to model,
and how far a rules baseline gets.

## Retrievals
| Date (UTC) | Source | Query | Records | Manifest |
|---|---|---|---|---|
| 2026-10-22 | BV-BRC genome_amr | eq(taxon_id,573), all fields in FIELDS | 1,840,793 (39,598 genomes) | snapshots/bvbrc-genome_amr-2026-10-22.json |

## Per-drug inventory
Full table: reports/spike/bvbrc_drug_inventory.csv (laboratory evidence only, unique genomes).

Evidence types in the snapshot:
- Computational Method: 1,728,894 records (excluded as labels, ADR-001)
- Laboratory Method: 85,291 records, 7,276 unique genomes
- Missing evidence: 26,608 records (excluded until inspected)

Shortlisted drugs (laboratory genomes, Resistant / Susceptible / Intermediate / Unknown / conflicting):
| Drug | Genomes | R | S | I | Unknown | Conflicting |
|---|---|---|---|---|---|---|
| meropenem | 6,074 | 1,444 | 2,551 | 159 | 1,980 | 60 |
| ciprofloxacin | 4,809 | 2,548 | 1,180 | 103 | 1,006 | 28 |
| gentamicin | 4,890 | 1,628 | 2,586 | 118 | 588 | 29 |
| colistin | 1,959 | 338 | 1,023 | 0 | 646 | 48 |
| tigecycline (stretch) | 1,690 | 262 | 690 | 117 | 698 | 77 |

## Obstacles and decisions
- "Unknown" laboratory phenotypes all carry a measurement value (broth/agar dilution MIC,
  some disk diffusion). They are uninterpreted measurements, not missing data, and can be
  labelled by applying breakpoints. See ADR-002.
- Interpretations come from mixed standards (CLSI ~15.8k, EUCAST ~8.5k, none ~3.5k records
  for shortlisted drugs). Mixing them adds label noise, especially for colistin.
- The `source` field is empty for all laboratory records, so study-level provenance needs a
  join to genome metadata (BioProject) in Sprint 1.
- Excluded for class imbalance or intrinsic resistance: ampicillin (7 S), cefazolin, cefuroxime.
- Snapshot size is ~520 MB of JSONL; raw stays out of git, manifest records the SHA-256.
