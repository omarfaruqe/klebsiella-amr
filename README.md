# klebsiella-amr

Genotype vs. phenotype in *Klebsiella pneumoniae*: a reproducible pipeline and a
leakage-aware benchmark.

## What this project does
- Curates public *K. pneumoniae* assemblies paired with laboratory AST results (BV-BRC, NCBI)
- Runs a pinned Nextflow pipeline: QC -> AMRFinderPlus + Kleborate -> MLST -> report
- Measures genotype-phenotype concordance per drug and investigates discordant isolates
- Benchmarks rules vs. logistic regression vs. gradient boosting under random,
  temporal and lineage-aware splits

## How it differs from existing tools
Pathogenwatch and similar services report genotypes. This project measures how well
those genotypes predict laboratory phenotype, and how much that estimate is inflated
by data leakage between related isolates.

## Status
Sprint 0: repository setup and data feasibility spike.

## Layout
    src/kamr/      Python package (ingestion, curation, analysis helpers)
    pipelines/     Nextflow workflow
    analysis/      notebooks
    snapshots/     provenance manifests for every raw data pull (tracked)
    reports/       generated tables and figures (tracked)
    data/          raw and curated data (not tracked)
    docs/          SRS, ADRs, spike log
