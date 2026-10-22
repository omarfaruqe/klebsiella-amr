# ADR-002: Target drugs and phenotype label policy

Status: Accepted
Date: 2026-10-22

## Context
The 2026-10-22 BV-BRC snapshot (docs/spike-log.md) has 7,276 K. pneumoniae genomes with
laboratory AST. Reported S/I/R calls mix CLSI and EUCAST, and every record without a call
still carries an MIC or zone measurement.

## Decision
Target drugs, chosen to cover distinct resistance mechanisms:
- meropenem: acquired carbapenemases; expected near-ceiling rules baseline (control)
- ciprofloxacin: QRDR point mutations (gyrA, parC) plus plasmid qnr / aac(6')-Ib-cr
- gentamicin: acquired aminoglycoside-modifying enzymes
- colistin: chromosomal mgrB / pmrAB / phoPQ changes plus rare mcr; expected rule failures
- tigecycline: stretch target only (small, noisy)

Label policy:
- Use laboratory-method records only.
- Re-interpret every available MIC with a single standard (EUCAST clinical breakpoints,
  version pinned in config) so labels are consistent across studies. Reported calls are kept
  for comparison, and agreement between reported and re-interpreted calls is reported.
- Binary target: Resistant vs Susceptible. Intermediate / "susceptible, increased exposure"
  is excluded from the primary benchmark and analysed separately.
- Genomes with conflicting labels for the same drug after re-interpretation are excluded
  and listed in the data-quality report.

## Consequences
- Labelled cohort grows by recovering "Unknown" records with MIC values.
- A breakpoint table becomes a versioned input that must be cited.
- Disk-diffusion-only records need zone breakpoints; if that is too costly they are dropped.
