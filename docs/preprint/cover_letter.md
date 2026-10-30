# Cover letter

To the Editors, *Microbial Genomics*

Dear Editors,

We submit for your consideration a manuscript entitled "Why genomic prediction of colistin
resistance does not transfer between *Klebsiella pneumoniae* collections".

Recent work has established that population structure confounds machine-learning prediction of
antimicrobial resistance (Yu et al., *PLOS Biology* 2025) and that published predictors lose
substantial accuracy under homology-aware validation (Hu et al., *Briefings in Bioinformatics*
2024). Our study differs from both in being mechanism-level rather than tool-level. We specify a
rule from a single known resistance mechanism rather than from a determinant catalogue, measure
directly what the leaked signal consists of rather than inferring it from a performance drop, and
validate the resulting rule, frozen, on an independent collection de-overlapped against our own.
That design lets us say not merely that apparent performance falls under stricter validation, but
what the models were using instead and why the loss occurs where it does.

Working from 6,652 public *K. pneumoniae* genomes with susceptibility results re-interpreted under
a single breakpoint standard, we find that a targeted *mgrB* integrity check — designed around the
assembly-break failure mode that gene screens miss — raises colistin sensitivity from 28.3% to
56.7% for AMRFinderPlus, and that the disrupted loci carry insertion-sequence signatures. The same
unfitted rule ranges from 39.5% to 70.6% sensitivity across held-out collections, so mechanistic
composition differs materially between public datasets. Machine-learning models exceed the rule
under random splitting but fall below it when sequence types are held out, a reversal that
survives correction across the family of 54 comparisons. A sequence-type-only predictor reaches
just 0.554, so that reversal is not simple lineage memorisation; and study-aware splits still
leave 74-88% sequence-type overlap while lineage-aware splits leave 87-95% BioProject overlap, so
neither scheme substitutes for the other. On 263 external isolates the frozen rule reaches 77.8% sensitivity at 99.5% specificity,
with the *mgrB* component performing as it did internally.

We believe the work suits *Microbial Genomics* because its contribution is simultaneously a
bacterial-genomics result about colistin resistance in *K. pneumoniae* and a methodological result
about how such predictions should be evaluated. It also reports two findings of practical use to
anyone assembling cohorts from public data: that six of 270 isolates in NCBI's *K. pneumoniae*
taxgroup are in fact *K. aerogenes*, and that temporal splitting constrains neither study nor
lineage overlap.

The manuscript is original, is not under consideration elsewhere, and all authors have approved
the submission. The authors declare no competing interests. The research received no external
funding. All code, pinned environments, provenance manifests with checksums, decision records and
derived tables are publicly available at https://github.com/omarfaruqe/klebsiella-amr and archived
at Zenodo; every reported figure regenerates deterministically from the documented pipeline.
Supplementary material is provided as a separate PDF.

We would be grateful for your consideration.

Yours sincerely,

Md Humayun Kabir (corresponding author)
Professor, Department of Computer Science and Engineering
University of Rajshahi, Rajshahi 6205, Bangladesh
hkj_cse@ru.ac.bd

on behalf of Md Omar Faruqe and Md Humayun Kabir
