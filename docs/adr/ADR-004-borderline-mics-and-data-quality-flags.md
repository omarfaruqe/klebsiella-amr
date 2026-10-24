# ADR-004: Borderline MICs and study-level data-quality flags

Status: Accepted
Date: 2026-10-24

## Context
The Sprint 2 discordance review (reports/review/) found two label problems that would distort the
benchmark if ignored:

1. Gentamicin: 213 of 299 rule false negatives with an MIC are at 4 mg/L, one dilution above the
   EUCAST v16.1 breakpoint (S <= 2, R > 2). Of the CLSI-tested records among them, 160 of 162 were
   reported susceptible by the original laboratory (older CLSI breakpoints call 4 mg/L S). The
   genotype rule agrees with the original call; re-interpretation under EUCAST made them R.
   Resistant isolates with a known mechanism mostly have MIC >= 8 (only 86 of 1,349 at 4).
   Without the borderline isolates, rule sensitivity is about 95% instead of 85%.
2. PRJEB28400: 63 of 65 gentamicin false positives carrying rmtB come from this one BioProject,
   almost all at MIC 1 mg/L. rmtB normally confers MICs > 256 mg/L, so these labels are
   implausible and likely a metadata error.

Meropenem shows the same pattern more mildly (ADR-002 already keeps I separate); its "false
positives" are mostly carbapenemase producers with MICs between the screening cut-off and the
clinical breakpoint (reports/review/meropenem.md).

## Decision
1. EUCAST v16.1 stays the primary standard (ADR-002).
2. Every benchmark result is reported twice: on all labels, and with "borderline" labels removed,
   where borderline means the highest MIC lies within one doubling dilution of the breakpoint on
   either side (gentamicin 2-4 mg/L; meropenem 2 and 16 mg/L are kept, as I is already excluded).
   The difference is reported as label uncertainty, not model error.
3. A study-level flag table (config/label_flags.csv) records label sets judged implausible by
   review. Flagged labels are excluded from the primary benchmark and listed in the report. First
   entry: PRJEB28400, gentamicin, reason "rmtB carriers at MIC 1 mg/L".
4. Flags are only added from a documented review with numbers, never to improve a metric.

## Consequences
- The benchmark needs an MIC-aware label table (highest MIC per genome and drug).
- Results become more honest about how much "error" is label noise near the breakpoint.
- Flag decisions are visible in git history and can be challenged.
