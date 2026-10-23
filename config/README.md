# Configuration

## breakpoints_eucast_v16.1.csv
Clinical MIC breakpoints (mg/L) used to re-interpret laboratory MICs (ADR-002).
S if MIC <= s_le, R if MIC > r_gt, otherwise I (susceptible, increased exposure).

Source: The European Committee on Antimicrobial Susceptibility Testing. Breakpoint tables for
interpretation of MICs and zone diameters. Version 16.1, valid from 2026-06-24.
https://www.eucast.org. Enterobacterales section.

Verified against the v16.1 PDF on 2026-10-23. Bracketed values (gentamicin, colistin) are
EUCAST's marker that the agent should be used in combination or with other measures; the
numeric breakpoints are used as given. Tigecycline has breakpoints only for E. coli and
C. koseri, so `applies_to_kpneumoniae=no` and reported calls are used instead.
