# ADR-003: Assembly source, colistin method policy and validation splits

Status: Accepted
Date: 2026-10-23

## Context
Sprint 1 part 2 QC (reports/cohort/qc_report.md) on 6,817 labelled genomes found:
- 6,655 pass assembly QC (2.4% removed, mostly > 500 contigs).
- 4,190 QC-passing genomes have no NCBI assembly accession (BV-BRC assembled them from reads).
- Colistin labels come from several methods: broth microdilution 367 R / 1,164 S; unspecified
  "MIC" 27 R / 250 S; disk diffusion 18 R / 68 S; agar dilution 6 R / 14 S; not reported 17 R / 26 S.
  EUCAST v16.1 requires broth microdilution for colistin MICs.
- Two BioProjects (PRJEB10018, PRJNA376414) contribute 44% of QC-passing genomes.
- Collection year is missing for 34% of genomes; 3,102 of 4,386 dated genomes are from 2010-2014
  and only 25 from 2020 onwards.

## Decision
1. Assemblies are downloaded from BV-BRC by genome_id for all genomes, so the cohort does not
   depend on NCBI accessions. NCBI accessions are kept as metadata where present.
2. Colistin primary label set: broth microdilution only. Sensitivity analysis adds unspecified
   "MIC" records. Disk diffusion, agar dilution and method-not-reported colistin labels are
   excluded from modelling and listed in the data-quality report.
3. Validation splits for the benchmark:
   - random (reported only as the optimistic reference),
   - study-aware: grouped by BioProject so no project appears in both train and test,
   - lineage-aware: grouped by ST (or Mash cluster in Sprint 2),
   - temporal: restricted to dated genomes, cut-off chosen so the test set is the most recent
     ~20%; reported with its limited size and date range.
   The gap between random and grouped splits is a headline result.
4. The pilot sample (per_stratum = 100, ST round-robin) is for tool validation and discordance
   review, not for prevalence estimates: it deliberately over-represents rare STs.

## Consequences
- Part 3 downloads FASTA from BV-BRC rather than NCBI Datasets.
- Colistin modelling set shrinks to about 1,530 genomes with 367 R.
- Split definitions become code with tests; metrics are reported per split type.
