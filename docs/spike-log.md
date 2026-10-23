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
