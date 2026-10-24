# klebsiella-amr

**How well do genes predict antibiotic resistance in *Klebsiella pneumoniae*, and how much do
standard validation choices overstate it?** A reproducible pipeline and benchmark on 6,652 public
genomes with laboratory susceptibility results.

## Key findings

1. **Standard tools miss most colistin resistance; a targeted *mgrB* check recovers part of it.**
   AMRFinderPlus explains 28% of broth-microdilution colistin resistance. Adding a check for
   disrupted *mgrB* raises this to 57%. Against Kleborate, which already calls MgrB truncations,
   the gain is smaller — 40% to 46% — and that is the honest increment over the state of the art.
   Many disruptions appear as assembly breaks inside the gene, consistent with insertion
   sequences; the association is not explained by assembly fragmentation (see the quartile
   analysis in `reports/review/mgrb_robustness.md`).

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="reports/figures/colistin_mgrb_dark.png">
  <img alt="Colistin sensitivity and specificity with 95% intervals for each rule, with and without the mgrB check" src="reports/figures/colistin_mgrb_light.png">
</picture>

2. **The same unfitted rule transfers unevenly between collections.** Leave-one-study-out,
   sensitivity ranges from 39.5% to 70.6%. Nothing is fitted, so this is not overfitting: the
   mechanistic composition of colistin resistance genuinely differs between isolate sets. In the
   largest collection, *mgrB* disruption accounts for under 40% of resistance. This bounds what
   any rule or model can achieve on new data.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="reports/figures/colistin_transfer_dark.png">
  <img alt="Sensitivity of the mgrB disruption rule in each held-out study, ranging from 39.5% to 70.6%" src="reports/figures/colistin_transfer_light.png">
</picture>

3. **Models beat the rule under random splitting and lose to it under lineage-aware splitting.**
   For colistin, logistic regression scores 0.844 balanced accuracy under random splits against
   the rule's 0.759, but 0.731 when sequence types are held out. Paired bootstrap intervals
   exclude zero in both directions, though the rule's lineage-aware margin over logistic
   regression is slim (+0.028, 95% CI +0.002 to +0.052) and is better read as parity. The models had the same *mgrB* features available; their
   advantage under laxer splitting came from lineage information, which will not be present on
   novel lineages.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="reports/figures/benchmark_splits_dark.png">
  <img alt="Balanced accuracy of rule, logistic regression and gradient boosting under random, study-aware and lineage-aware splits for five drugs" src="reports/figures/benchmark_splits_light.png">
</picture>

4. **A plausible determinant tier did not survive out-of-sample selection.** *mgrB* missense
   variants enriched among resistant isolates — including substitutions at both structural
   cysteines — appeared to add 19.1 points of sensitivity in-sample and 0.9 points when selected
   leave-one-study-out. The tier was discarded. It is reported as a cautionary result.

5. **Many "errors" are label effects, not missing biology.** 93% of carbapenemase carriers that
   test meropenem-susceptible exceed the EUCAST carbapenemase screening cut-off; 71% of
   unexplained gentamicin-resistant isolates sit one dilution above the breakpoint, where 160 of
   162 CLSI-tested records had been called susceptible by the originating laboratory.

## How it differs from existing tools
Pathogenwatch, Kleborate and AMRFinderPlus report genotypes. This project measures how well those
genotypes predict laboratory phenotype on harmonised labels, explains the disagreements, and
quantifies how far the estimate depends on which collection the test genomes came from.

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
| Paired comparisons between predictors | `src/kamr/benchmark/compare.py` | `reports/benchmark/comparisons.md` |
| mgrB robustness checks | `src/kamr/analysis/mgrb_robustness.py` | `reports/review/mgrb_robustness.md` |
| Figures | `src/kamr/reporting/figures.py` | `reports/figures/` |

Decisions are recorded as ADRs in `docs/adr/`, and every run is logged in `docs/spike-log.md`.
The model card is `docs/model-card.md`; a manuscript draft is `docs/preprint/manuscript.md`.

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
python -m kamr.analysis.mgrb_robustness
python -m kamr.analysis.review_gentamicin
python -m kamr.benchmark.features
python -m kamr.benchmark.run
python -m kamr.benchmark.compare
python -m kamr.reporting.figures
```

Tests: `pytest -q` (75 tests). Lint and format: `ruff check . && ruff format --check .`.

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
