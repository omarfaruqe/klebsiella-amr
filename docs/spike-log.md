# Feasibility spike log

Goal: decide which drugs have enough paired genome + laboratory AST data to model,
and how far a rules baseline gets.

## Retrievals
| Date (UTC) | Source | Query | Records | Manifest |
|---|---|---|---|---|
| 2026-10-22 | BV-BRC genome_amr | eq(taxon_id,573), all fields in FIELDS | 1,840,793 (39,598 genomes) | snapshots/bvbrc-genome_amr-2026-10-22.json |
| 2026-10-23 | BV-BRC genome_amr (server re-pull) | same query | 1,840,793 (39,598 genomes) | snapshots/bvbrc-genome_amr-2026-10-23.json |

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

## Toolchain validation: HS11286 (GCF_000240185.1), 2026-10-22
Versions: AMRFinderPlus 4.2.7 (database 2026-08-07.1, --organism Klebsiella_pneumoniae --plus),
Kleborate 3.1.3 (-p kpsc) with Kaptive 3.2.2.

Assembly: 7 contigs, 5,682,322 bp, N50 5,333,942; Kleborate QC flags 1 ambiguous base.

| Marker | AMRFinderPlus | Kleborate | Note |
|---|---|---|---|
| ST | n/a | ST11 | expected for HS11286 |
| K / O locus | n/a | KL103 / OL2a.1 | typeable |
| Carbapenemase | blaKPC-2 | KPC-2 | agree |
| ESBL | blaCTX-M-14 | CTX-M-14 | agree |
| QRDR | gyrA_S83I, parC_S80I | GyrA-83I, ParC-80I | agree |
| Porin | ompK35 frameshift | OmpK35-40% | agree |
| Colistin | none reported | PmrB-10% (truncation) | DISCORDANT; drives Kleborate resistance score 3 |
| Gentamicin | aac(3)-IId, rmtB1 | rmtB (no aac(3)-IId in AGly column) | DISCORDANT on aac(3)-IId |

Observations:
- Kleborate counts the PmrB truncation as a colistin marker. Colistin resistance via PmrB is
  usually gain-of-function missense; a truncation is more likely loss of function, so this may
  be a false-positive call. Colistin rules must be evaluated per tool against phenotype.
- Rule baselines will be built separately from each tool and from their union, and the
  tool-vs-tool disagreement rate is reported per drug.

## Environment notes
- Kleborate 3.1.3 imports kaptive.database; Kaptive 3.3 renamed it to kaptive.db and neither
  package caps the version. Pinned kaptive>=3.0,<3.3 in environment.yml (3.2.2 works).
- AMRFinderPlus took 273 s on one genome with default threads on the M4 Air.
- Server (Ubuntu x86_64, 32 cores, 31 GB RAM) rebuilt from environment.yml resolved Kleborate
  3.2.4; HS11286 results matched the Mac (ST11, resistance score 3, 9 classes, 13 genes;
  AMRFinderPlus blaKPC-2, gyrA_S83I, parC_S80I, aac(3)-IId). Pinned kleborate=3.2.4 and
  kaptive=3.2.2 exactly so both machines resolve the same versions.
- AMRFinderPlus nucleotide mode: 270 s/genome on the M4 Air (single core) vs 90 s on the
  server (multi-threaded). Batch runs go on the server.
- Installing Kleborate pulls ECTyper, whose post-link script downloads a ~990 MB Mash sketch
  from Zenodo; first environment builds take extra time for this.
- bioconda kleborate=3.2.4 requires stxtyper, which has no osx-arm64 build. The Mac
  installs Kleborate 3.2.4 and Kaptive 3.2.2 from PyPI (environment-macos-arm64.yml);
  the server uses environment.yml. Equivalence checked with scripts/verify_toolchain.sh.

## Cross-machine verification, 2026-10-23
`scripts/verify_toolchain.sh` on HS11286 (FASTA sha256 prefix 9528cd6940c58fae):

| Machine | OS / arch | Kleborate source | AMRFinderPlus | Kleborate | Fingerprint |
|---|---|---|---|---|---|
| Mac (faruqe) | Darwin arm64, 10 cores | PyPI (environment-macos-arm64.yml) | 269 s | 16 s | 131c4f63ab10c2b3 |
| Server (omf) | Linux x86_64, 32 cores | bioconda (environment.yml) | 63 s (8 threads) | 16 s | 131c4f63ab10c2b3 |

Both: AMRFinderPlus 4.2.7, database 2026-08-07.1, Kleborate 3.2.4, Kaptive 3.2.2, Mash 2.3,
seqkit 2.14.0, BLAST 2.17.0. Fingerprint is identical at 4 and 8 threads on the Mac, so thread
count does not change results. Decision: development on the Mac, batch runs on the server.

Kleborate 3.2 changed mutation notation from 3.1 (e.g. PmrB-10% -> pmrb:p.Gly20fs,
GyrA-83I -> GyrA:p.Ser83Ile, OmpK35-40% -> OmpK35:p.Trp133fs). AMRFinderPlus and Kleborate
also number the same OmpK35 frameshift differently (G134AfsTer11 vs Trp133fs). Feature
extraction must normalise variants to gene + variant class, not match raw strings.

## Snapshot reproducibility, 2026-10-23
The Mac's raw copy of the 2026-10-22 snapshot was lost when the project folder was replaced by
a GitHub archive (raw data is gitignored by design). The snapshot was re-downloaded on the server
with the same downloader. Result: 1,840,793 records, 39,598 genomes and an identical SHA-256
(d0e561a6...1505126). The sorted, key-ordered JSONL output makes pulls byte-reproducible while
BV-BRC content is unchanged. Raw data now lives on the server; the Mac keeps only manifests.

Sprint 0 complete.

# Sprint 1

## Part 1: labelled cohort (2026-10-23, server, 14 s)
`python -m kamr.curation.cohort data/raw/bvbrc/2026-10-23/genome_amr.jsonl` with
config/breakpoints_eucast_v16.1.csv (verified). 19,803 laboratory records for the 5 drugs.
Full table: reports/cohort/label_report.md.

| Drug | R before -> after | S before -> after | Conflicts before -> after | Recovered from Unknown | Modelling set (S+R) |
|---|---|---|---|---|---|
| meropenem | 1,444 -> 1,886 | 2,551 -> 3,798 | 60 -> 8 | 1,921 | 5,684 |
| gentamicin | 1,628 -> 2,147 | 2,586 -> 2,699 | 29 -> 7 | 560 | 4,846 |
| ciprofloxacin | 2,548 -> 3,501 | 1,180 -> 1,228 | 28 -> 1 | 959 | 4,729 |
| colistin | 338 -> 429 | 1,023 -> 1,528 | 48 -> 2 | 598 | 1,957 |
| tigecycline | 262 -> 262 | 690 -> 690 | 77 -> 0 | 0 | 952 |

"Before" = reported calls in the Sprint 0 inventory; "after" = EUCAST v16.1 re-interpretation.

Findings:
- Re-interpreting MIC-only "Unknown" records adds 30-50% labelled genomes per drug.
- Most apparent conflicts were mixed standards (CLSI vs EUCAST), not contradictory results;
  a single standard removes nearly all of them.
- Reported vs re-interpreted agreement is 92.9-97.2% where both exist. Disagreements are
  mostly standard differences (e.g. CLSI meropenem R from 4 mg/L vs EUCAST I for 4-8 mg/L,
  raising meropenem I from 159 to 382).
- Colistin agreement (100%) rests on only 269 records; most labels come from bare MICs.
  Testing method must be checked (EUCAST requires broth microdilution).
- Tigecycline: 621 MIC-only genomes stay unlabelled (no K. pneumoniae breakpoint).
- Ciprofloxacin: 19 genomes unlabelled because censored MICs are undecidable.

## Part 2: genome metadata and QC (2026-10-23, server)
Metadata for 6,817 labelled genomes (snapshots/bvbrc-genome-2026-10-23.json); 6,655 pass QC.
Full report: reports/cohort/qc_report.md. Decisions recorded in ADR-003:
- 61% of QC-passing genomes lack an NCBI assembly accession, so assemblies come from BV-BRC.
- Colistin primary labels restricted to broth microdilution (367 R / 1,164 S).
- Two BioProjects hold 44% of genomes: study-aware splits are required.
- Collection year missing for 34%; dated genomes concentrate in 2010-2014.
- First pilot (per_stratum = 50) gave 265 genomes over 223 STs; raised to per_stratum = 100.

## Part 3: pilot assemblies and tool run (2026-10-24 UTC, server)
- 524 pilot assemblies downloaded from the BV-BRC genome_sequence API in 3 min; all passed
  format and length checks (snapshots/assemblies-pilot-2026-10-24.tsv on the server). The
  BV-BRC FTP route timed out from the server, so the API is tried first.
- AMRFinderPlus + Kleborate, 8 jobs x 4 threads on 32 cores: 49 min wall time, 12 h CPU.
  523 of 524 genomes completed both tools.
- Failure: 573.46141. Kleborate rejected it as outside the K. pneumoniae species complex; its
  species module calls Citrobacter freundii (strong match). BV-BRC metadata contradicts itself:
  species "Klebsiella pneumoniae" but genome_name "Citrobacter freundii PBIO1966".
- QC fix: the species check now requires both the species field and the genome name to match
  (steps species_field and genome_name_species). This removes 3 of 6,655 QC-passing genomes
  (2 E. coli, 1 C. freundii); only 573.46141 is in the pilot sample.
- The pilot sample is not re-drawn, because the tool run used it; 573.46141 is excluded from
  pilot analysis. Kleborate's sequence-based species call is the final species QC.

## Part 4: pilot concordance (2026-10-24 UTC, server)
`python -m kamr.analysis.concordance` on 523 pilot genomes (573.46141 excluded).
Full tables: reports/pilot/concordance.md, reports/pilot/concordance.csv; discordant
genome-drug pairs with each tool's evidence: reports/pilot/discordant.tsv.

| Drug | Sensitivity AMRFinderPlus / Kleborate | Specificity AMRFinderPlus / Kleborate |
|---|---|---|
| meropenem | 92.3% / 92.3% | 85.7% / 86.1% |
| ciprofloxacin | 94.7% / 94.7% | 92.7% / 93.5% |
| gentamicin | 91.6% / 83.8% | 94.2% / 97.9% |
| colistin (broth microdilution) | 24.7% / 44.2% | 95.3% / 95.3% |
| tigecycline | 2.0% / 0.0% (19.0% with regulators) | 100% / 100% (85.0% with regulators) |

Findings:
- First run of the rules exposed two rule errors, not tool errors: AMRFinderPlus assigns the
  intrinsic chromosomal oqxAB pump to tigecycline as well as quinolones (88 of 91 tigecycline
  false positives), and it reports efflux-regulator variants (oqxR, ramR, acrR) that are mostly
  lineage polymorphisms (21 of 26 ciprofloxacin false positives). Fix: oqxAB ignored for all
  drugs; regulator variants moved to a separate "AMRFinderPlus+regulators" rule.
- Regulator variants add sensitivity only for tigecycline (2% -> 19%, 15 extra false positives);
  for ciprofloxacin they add 1.1 points of sensitivity and cost 13.7 points of specificity.
- Kleborate reports only acquired tigecycline genes (Tgc_acquired) and no regulator mutations.
- Meropenem, ciprofloxacin and gentamicin are largely explained by known mechanisms (controls).
- Colistin and tigecycline resistance is mostly unexplained by either tool: the targets where a
  learned model could add something in the benchmark. Expected colistin gap: mgrB disruption by
  insertion sequences, which assembly-based gene and point-mutation searches miss.
- Meropenem false positives (carbapenemase present, phenotype S) and false negatives (R without
  a carbapenemase) are listed for review in discordant.tsv.

## Part 5: full cohort (2026-10-24 UTC, server)
- QC re-run with the genome-name species check: 6,652 genomes pass (was 6,655); pilot sample
  kept frozen. List: reports/cohort/qc_pass_genomes.tsv.
- 6,652 assemblies downloaded from BV-BRC in 34 min, all valid
  (snapshots/assemblies-cohort-2026-10-24.tsv).
- AMRFinderPlus + Kleborate: 8 jobs x 4 threads, 530 min wall time, about 140 CPU-hours,
  0 failures; 523 pilot results reused. Summary: reports/cohort/run_summary.md.

Cohort concordance (reports/cohort/concordance.md; all QC-passing labelled genomes, not
de-duplicated for related isolates):

| Drug | n (R/S) | Sensitivity AMRFinderPlus / Kleborate | Specificity AMRFinderPlus / Kleborate |
|---|---|---|---|
| meropenem | 5,541 (1,854/3,687) | 96.7% / 96.5% | 91.2% / 91.3% |
| ciprofloxacin | 4,640 (3,441/1,199) | 98.8% / 98.7% | 88.2% / 88.4% |
| gentamicin | 4,760 (2,113/2,647) | 85.1% / 81.4% | 94.9% / 95.8% |
| colistin (BMD) | 1,528 (367/1,161) | 28.3% / 39.5% | 97.8% / 96.3% |
| tigecycline | 948 (261/687) | 2.7% / 0.0% (12.6% with regulators) | 100% / 100% (84.4% with regulators) |

Findings to review (reports/cohort/discordant.tsv, 1,740 genome-drug pairs):
- meropenem: 326 isolates with a carbapenemase but a susceptible phenotype.
- gentamicin: 314 resistant isolates with no detected gentamicin mechanism.
- colistin: 263 broth-microdilution resistant isolates unexplained by AMRFinderPlus.
- The pilot's conclusions hold on the full cohort with much narrower intervals.

# Sprint 2: discordance review

## Part 1: meropenem (2026-10-24, server; reports/review/meropenem.md)
- 2,118 carbapenemase carriers; 326 test S (rule false positives). 3,423 non-carriers; 62 test R.
- Share testing S by family: KPC 5.9% (75/1,278), NDM 10.0% (28/280), OXA-48-like 52.9%
  (182/344), VIM 75.0% (33/44), IMP 27.3% (6/22), multiple 0.7% (1/149). OXA-48-like and VIM
  account for 215 of 326 false positives.
- 240 of the 326 have an MIC: 191 (80%) at 1-2 mg/L and only 16 at or below 0.125 mg/L, the
  EUCAST carbapenemase screening cut-off. Most "false positives" are carbapenemase producers
  that are clinically susceptible, not genotype errors.
- Resistant non-carriers: 51 of 62 carry an ESBL/AmpC (29 with an OmpK35/36 defect, 22 without);
  but 624 susceptible non-carriers have ESBL/AmpC plus a porin defect, so the combination is a
  risk factor rather than a predictor.
- PRJEB10018 contributes 158 of 326 false positives (48%); ST15 is the commonest ST among them.
  Strong study clustering, supporting study-aware splits (ADR-003).

## Part 2: colistin and mgrB (2026-10-24, server; reports/review/colistin_mgrb.md)
- tblastn of MgrB (WP_002911375.1, from the AMRFinderPlus database) against all 6,652
  assemblies: 2 min on 16 jobs. Per-genome status: reports/review/mgrb_status.tsv.
- First run labelled 6,048 genomes "missense" at 97.9% identity: one artefactual mismatch at
  position 1, because tblastn translates the alternative start codon as V/L. Fixed by calling
  variants residue by residue and ignoring position 1. Now 6,056 intact, 129 missense, and
  467 disrupted (split_contigs 162, partial 153, premature_stop 73, absent 34,
  frameshift_or_indel 27, contig_break 13, insertion_assembled 5).
- Colistin, broth microdilution (n = 1,528; 367 R): disrupted mgrB 153 R / 41 S;
  missense 70 R / 13 S (e.g. M27K 25/3, T21N 9/0, N42Y 7/0, K43I 7/0, C28W 4/1, C39G 3/0);
  intact 144 R / 1,107 S.

| Rule | Sensitivity | Specificity |
|---|---|---|
| AMRFinderPlus | 28.3% | 97.8% |
| AMRFinderPlus + mgrB disrupted | 56.7% | 95.2% |
| AMRFinderPlus + mgrB disrupted or missense | 67.3% | 94.3% |
| Kleborate | 39.5% | 96.3% |
| Kleborate + mgrB disrupted or missense | 61.6% | 94.5% |
| mgrB disrupted or missense alone | 60.8% | 95.3% |

- A targeted mgrB check explains 143 more resistant isolates than AMRFinderPlus alone, for 40
  extra false positives. Many disruptions appear as contig breaks or split contigs, consistent
  with insertion sequences that short-read assemblies cannot span.
- 120 resistant isolates remain unexplained; most have an intact mgrB coding sequence. Next
  candidates: pmrB/phoQ/crrB variants and IS insertions in the mgrB promoter (not searched).

## Part 3: gentamicin (2026-10-24, server; reports/review/gentamicin.md)
- 4,760 labelled genomes: TP 1,799, FN 314, FP 135, TN 2,512 (AMRFinderPlus core rule).
- MIC: 213 of 299 false negatives with an MIC are at 4 mg/L, one dilution above the EUCAST
  breakpoint; true positives are mostly >= 8 (86 of 1,349 at 4). Among CLSI-tested records of
  these borderline false negatives, 160 of 162 were reported S by the original laboratory.
  The rule agrees with the original call; EUCAST re-interpretation made them R. Without the
  borderline isolates, sensitivity is about 95% rather than 85%.
- Genes: aac(6')-Ib (63% of FN vs 12% of susceptible non-carriers) and aph(3')-Ia (61% vs 14%)
  are enriched in false negatives; classical aac(6')-Ib does not confer gentamicin resistance, so
  this is more likely lineage or plasmid co-carriage than a missed mechanism. AMRFinderPlus assigns
  the GENTAMICIN subclass to some aac(6')-Ib alleles (Ib', Ib4) and not others.
- Study: PRJNA376414 contributes 207 of 314 false negatives; in PRJNA504600 (24/24) and
  PRJNA385863 (23/24) every resistant isolate is unexplained.
- False positives: 63 of 65 rmtB carriers reported S come from PRJEB28400, nearly all at MIC
  1 mg/L; implausible for rmtB, flagged as a likely metadata error.
- Decisions in ADR-004: report every benchmark result with and without borderline-MIC labels;
  keep a reviewed study-level flag table (config/label_flags.csv), first entry PRJEB28400.

# Sprint 3: leakage-aware benchmark

## Part 1: inputs (2026-10-24; reports/benchmark/feature_summary.md)
- 6,652 genomes x 226 binary features: gene families (OXA-48-like kept apart from other OXA),
  100 point mutations, 33 AMRFinderPlus subclasses, 2 mgrB features.
- First version lumped OXA-1 with OXA-48-like carbapenemases in one "blaOXA" feature; fixed
  before any modelling. Subclass features separate, e.g., ESBL SHV from chromosomal SHV-1.

## Part 2: benchmark (2026-10-24, server, 52 s; reports/benchmark/results.md)
Balanced accuracy, all test labels (flagged labels excluded):

| Drug | Rule | Best model, random | Best model, study | Best model, lineage |
|---|---|---|---|---|
| ciprofloxacin | 0.935 | 0.954 | 0.948 | 0.950 |
| meropenem | 0.939 | 0.951 | 0.942 | 0.945 |
| gentamicin | 0.906 | 0.925 | 0.886 | 0.916 |
| colistin | 0.808 | 0.844 | 0.822 | 0.731 |
| tigecycline | 0.513 | 0.769 | 0.626 | 0.716 |

- Well-explained drugs: models add 1-2 points on random splits, mostly lost under grouped
  splits; gentamicin rule beats both models under the study split, and without borderline
  labels (0.959 vs <= 0.948).
- Colistin: models fall below the rule (0.71-0.73 vs 0.81) under the lineage split; random-split
  optimism 11-13 points. Tigecycline: apparent model gain over the rule is mostly leakage;
  optimism 14-15 points.
- Temporal split (train <= 2014): meaningful only for meropenem, ciprofloxacin and gentamicin
  (546-569 test genomes). Colistin (34 test genomes, 1 R) and tigecycline (28, 0 R) are not
  evaluable; first report printed degenerate values (1.000 / 0.500), now reported as n/a with a
  minimum of 10 R and 10 S per test set.
- Headline: on public AMR data, random splits overstate model performance most where rules
  fail (colistin, tigecycline), because models learn lineage and study markers.
