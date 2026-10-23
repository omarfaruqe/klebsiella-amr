# Cohort QC and data-quality report

- Generated (UTC): 2026-10-23T23:41:37+00:00
- Genome metadata: `data/raw/bvbrc_genome/2026-10-23/genome.jsonl`
- Thresholds: `config/qc.toml`

## QC funnel

| Step | Remaining | Removed |
|---|---|---|
| labelled genomes with metadata | 6,817 | 0 |
| species | 6,817 | 0 |
| assembly_stats_present | 6,817 | 0 |
| min_length | 6,801 | 16 |
| max_length | 6,774 | 27 |
| max_contigs | 6,674 | 100 |
| min_n50 | 6,674 | 0 |
| checkm_completeness | 6,670 | 4 |
| checkm_contamination | 6,655 | 15 |

CheckM values missing (not failed): 0 genomes.
Genomes without an assembly accession: 4,190.

## Modelling labels after QC

| Drug | R | S | I | R share of R+S |
|---|---|---|---|---|
| meropenem | 1,854 | 3,687 | 372 | 33.5% |
| ciprofloxacin | 3,441 | 1,199 | 60 | 74.2% |
| gentamicin | 2,113 | 2,647 | 36 | 44.4% |
| colistin | 418 | 1,486 | 0 | 22.0% |
| tigecycline | 261 | 687 | 116 | 27.5% |

## Colistin labels by testing method (QC-passing genomes)

EUCAST requires broth microdilution for colistin MICs.

```
label             R     S  conflict
method                             
Agar dilution     6    14         0
Broth dilution  367  1164         2
Disk diffusion   18    68         1
MIC              27   250         0
not reported     17    26         1
```

## Collection year (QC-passing genomes)

Missing year: 2,269 of 6,655.

| Period | Genomes |
|---|---|
| <=2004 | 202 |
| 2005-2009 | 452 |
| 2010-2014 | 3,102 |
| 2015-2019 | 605 |
| >=2020 | 25 |

## Top 10 sequence types

| Sequence types | Genomes | Share |
|---|---|---|
| ST307 | 762 | 11.5% |
| ST258 | 669 | 10.1% |
| ST15 | 484 | 7.3% |
| ST11 | 382 | 5.7% |
| ST512 | 298 | 4.5% |
| ST16 | 201 | 3.0% |
| ST147 | 194 | 2.9% |
| ST101 | 178 | 2.7% |
| ST37 | 130 | 2.0% |
| ST14 | 128 | 1.9% |

## Top 10 isolation country

| Isolation country | Genomes | Share |
|---|---|---|
| USA | 1,838 | 27.6% |
| Norway | 652 | 9.8% |
| Italy | 587 | 8.8% |
| Vietnam | 357 | 5.4% |
| Ireland | 252 | 3.8% |
| United Kingdom | 245 | 3.7% |
| Spain | 222 | 3.3% |
| Australia | 198 | 3.0% |
| Japan | 166 | 2.5% |
| Turkey | 154 | 2.3% |

## Top 10 bioproject

| BioProject | Genomes | Share |
|---|---|---|
| PRJEB10018 | 1,476 | 22.2% |
| PRJNA376414 | 1,467 | 22.0% |
| PRJEB27256 | 648 | 9.7% |
| PRJEB28400 | 332 | 5.0% |
| missing | 329 | 4.9% |
| PRJEB31361 | 248 | 3.7% |
| PRJEB5065 | 245 | 3.7% |
| PRJDB10842 | 166 | 2.5% |
| PRJEB24082 | 159 | 2.4% |
| PRJEB6891 | 156 | 2.3% |

## Pilot sample

524 genomes in `reports/cohort/pilot_sample.tsv`, 331 distinct STs.

| Drug | R in sample | S in sample |
|---|---|---|
| meropenem | 169 | 259 |
| ciprofloxacin | 190 | 124 |
| gentamicin | 154 | 243 |
| colistin | 104 | 111 |
| tigecycline | 100 | 100 |
