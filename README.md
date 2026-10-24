# klebsiella-amr

**How well do genes predict antibiotic resistance in *Klebsiella pneumoniae*, and how much do
standard validation choices overstate it?** A reproducible pipeline and benchmark on 6,652 public
genomes with laboratory susceptibility results.

## Key findings

1. **Random train/test splits overstate model performance most where rules fail.** For colistin
   and tigecycline, models look 11-15 points better under random splits than when whole studies
   or lineages are held out. For meropenem, ciprofloxacin and gentamicin, simple rules match or
   beat the models once leakage is removed.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="reports/figures/benchmark_splits_dark.png">
  <img alt="Balanced accuracy of rule, logistic regression and gradient boosting under random, study-aware and lineage-aware splits for five drugs" src="reports/figures/benchmark_splits_light.png">
</picture>

2. **Standard tools miss most colistin resistance; a targeted *mgrB* check recovers it.**
   AMRFinderPlus explains 28% of broth-microdilution colistin resistance. Adding a check for
   disrupted or mutated *mgrB* raises this to 67%, at a cost of 3.5 points of specificity. Many
   disruptions show up as assembly breaks inside the gene, consistent with insertion sequences.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="reports/figures/colistin_mgrb_dark.png">
  <img alt="Colistin sensitivity and specificity with 95% intervals for each rule, with and without the mgrB check" src="reports/figures/colistin_mgrb_light.png">
</picture>

3. **Many "errors" are label effects, not missing biology.**
   - Meropenem: 93% of carbapenemase carriers that test susceptible (with an MIC) are above the
     EUCAST carbapenemase screening cut-off; most carry OXA-48-like or VIM enzymes.
   - Gentamicin: 71% of unexplained resistant isolates sit at 4 mg/L, one dilution above the
     EUCAST breakpoint; the original CLSI-based labs had called 160 of 162 of them susceptible.
   - One study's gentamicin labels are biologically implausible (rmtB carriers at MIC 1 mg/L)
     and are flagged.

## How it differs from existing tools
Pathogenwatch, Kleborate and AMRFinderPlus report genotypes. This project measures how well those
genotypes predict laboratory phenotype on harmonised labels, explains the disagreements, and
shows how much the estimate depends on keeping related isolates out of the test set.

## Pipeline

| Step | Code | Output |
|---|---|---|
| Snapshot BV-BRC AST and genome metadata | `src/kamr/ingestion/` | `snapshots/` manifests with SHA-256 |
| Re-interpret MICs under EUCAST v16.1 | `src/kamr/curation/cohort.py` | `reports/cohort/label_report.md` |
| Assembly QC, cohort list | `src/kamr/curation/qc.py` | `reports/cohort/qc_report.md` |
| Download assemblies, run AMRFinderPlus and Kleborate | `src/kamr/ingestion/assemblies.py`, `src/kamr/pipeline/` | `reports/cohort/run_summary.md` |
| Genotype-phenotype concordance | `src/kamr/analysis/concordance.py` | `reports/cohort/concordance.md` |
| Discordance review (meropenem, colistin/mgrB, gentamicin) | `src/kamr/analysis/` | `reports/review/` |
| Leakage-aware benchmark | `src/kamr/benchmark/` | `reports/benchmark/results.md` |
| Figures | `src/kamr/reporting/figures.py` | `reports/figures/` |

Decisions are recorded as ADRs in `docs/adr/`, and every run is logged in `docs/spike-log.md`.
The model card is `docs/model-card.md`.

## Reproduce the results

Tested on Ubuntu x86_64 (32 cores) and macOS arm64; a cross-machine check
(`scripts/verify_toolchain.sh`) gives identical calls on both.

```bash
git clone https://github.com/omarfaruqe/klebsiella-amr.git && cd klebsiella-amr
mamba env create -f environment.yml          # macOS arm64: environment-macos-arm64.yml
conda activate kamr && python -m pip install -e . --no-deps && amrfinder -u
scripts/verify_toolchain.sh                  # expect FINGERPRINT 131c4f63ab10c2b3

python -m kamr.ingestion.bvbrc_amr                                   # AST snapshot (~0.5 GB)
python -m kamr.curation.cohort data/raw/bvbrc/<date>/genome_amr.jsonl
python -m kamr.ingestion.bvbrc_genome
python -m kamr.curation.qc data/raw/bvbrc_genome/<date>/genome.jsonl
python -m kamr.ingestion.assemblies --sample reports/cohort/qc_pass_genomes.tsv \
    --manifest-name assemblies-cohort --workers 8                    # ~38 GB
python -m kamr.pipeline.run_tools snapshots/assemblies-cohort-<date>.tsv \
    --results results/cohort --jobs 8 --threads 4                    # ~9 h on 32 cores
python -m kamr.pipeline.collate snapshots/assemblies-cohort-<date>.tsv \
    --results results/cohort --report reports/cohort/run_summary.md --title "Cohort run summary"
python -m kamr.analysis.concordance --results results/cohort \
    --sample reports/cohort/qc_pass_genomes.tsv --reports reports/cohort --population \
    --title "Cohort genotype-phenotype concordance"
python -m kamr.analysis.review_meropenem
python -m kamr.analysis.mgrb --jobs 16
python -m kamr.analysis.review_gentamicin
python -m kamr.benchmark.features
python -m kamr.benchmark.run
python -m kamr.reporting.figures
```

Tests: `pytest -q` (73 tests). Lint and format: `ruff check . && ruff format --check .`.

## Layout
    src/kamr/      Python package (ingestion, curation, pipeline, analysis, benchmark, reporting)
    config/        breakpoints (EUCAST v16.1), QC thresholds, label flags, MgrB reference
    snapshots/     provenance manifests for every data pull
    reports/       generated tables, reviews, benchmark results and figures
    docs/          ADRs, spike log, model card
    data/, results/  raw data, assemblies and tool outputs (not tracked)

## Data and tools
Genomes and AST: BV-BRC. Breakpoints: EUCAST Clinical Breakpoint Tables v16.1. Tools:
AMRFinderPlus 4.2.7 (database 2026-08-07.1), Kleborate 3.2.4 with Kaptive 3.2.2, BLAST 2.17,
scikit-learn 1.9.1. Research use only; not for clinical decisions.
