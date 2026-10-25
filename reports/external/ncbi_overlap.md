# External colistin set from NCBI Pathogen Detection

- Generated (UTC): 2026-10-25T02:11:01+00:00
- Export: `data/external/ncbi/ast_export.csv` (355 rows, 355 after deduplication)
- Breakpoint: EUCAST 16.1

## Overlap with the BV-BRC cohort

- Already in our cohort: 84 of 355 (23.7%)
- Held out as external: 271
- With a usable label: 270
- Distinct assemblies to download: 270

Matched on BioSample accession or on assembly accession ignoring the version suffix.
The overlap fraction bounds how independent any NCBI-derived validation set can be,
given that BV-BRC ingests from NCBI.

## Label provenance

| Source | Isolates |
|---|---|
| mic_eucast | 257 |
| submitted | 13 |

Labels from `mic_eucast` were re-derived from the submitted MIC under our pinned EUCAST
table. Labels from `submitted` had no usable MIC and retain the submitter's own
interpretation; they are reported as a separate stratum so the headline result does not
rest on foreign breakpoint choices.

## Composition

| Label | Isolates |
|---|---|
| R | 82 |
| I | 1 |
| S | 187 |

Testing standards as submitted:

| Standard | Isolates |
|---|---|
| EUCAST | 140 |
| CLSI | 112 |
| (not stated) | 18 |
