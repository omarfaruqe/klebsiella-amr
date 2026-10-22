# ADR-001: Scope is K. pneumoniae assemblies only

Status: Accepted
Date: 2026-10-21

## Context
The original plan covered two organisms, raw reads, a multi-user web platform and a GPU
stack. Available time is about 15 hours per week, roughly 450 hours in total.

## Decision
- Organism: Klebsiella pneumoniae only.
- Input: genome assemblies only; no raw FASTQ processing.
- Labels: laboratory-method AST only; computational predictions are excluded as labels.
- Web layer: read-only, precomputed results; no uploads or accounts.
- Compute: CPU only. GPU use limited to an optional protein-embedding experiment.

## Consequences
Deferred: E. coli, read-level QC and assembly, user uploads, maps, monitoring stack.
These can be revisited after the benchmark is complete.
