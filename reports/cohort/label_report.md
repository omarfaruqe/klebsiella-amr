# Cohort label report

- Generated (UTC): 2026-10-23T23:19:20+00:00
- Snapshot: `data/raw/bvbrc/2026-10-23/genome_amr.jsonl`
- Breakpoints: `config/breakpoints_eucast_v16.1.csv`
- Laboratory records for target drugs: 19,803

Label policy (ADR-002): laboratory records only; MICs re-interpreted with one breakpoint
table where it applies to K. pneumoniae, otherwise the reported call is used; genomes
with both S and R for a drug are conflicts and excluded from modelling; I is kept
separately.

| Drug | Genomes | R | S | I | Conflict | Unlabelled | Recovered from Unknown | Reported vs EUCAST agreement |
|---|---|---|---|---|---|---|---|---|
| meropenem | 6,074 | 1,886 | 3,798 | 382 | 8 | 0 | 1,921 | 94.5% (n=2,544) |
| gentamicin | 4,890 | 2,147 | 2,699 | 37 | 7 | 0 | 560 | 92.9% (n=2,525) |
| ciprofloxacin | 4,809 | 3,501 | 1,228 | 60 | 1 | 19 | 959 | 97.2% (n=2,026) |
| colistin | 1,959 | 429 | 1,528 | 0 | 2 | 0 | 598 | 100.0% (n=269) |
| tigecycline | 1,690 | 262 | 690 | 117 | 0 | 621 | 0 | n/a |
