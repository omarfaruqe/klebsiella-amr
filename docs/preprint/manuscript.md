# Why genomic prediction of colistin resistance does not transfer between *Klebsiella pneumoniae* collections

Md Omar Faruqe and Md Humayun Kabir

Department of Computer Science and Engineering, University of Rajshahi, Rajshahi 6205, Bangladesh

Correspondence: Md Humayun Kabir (hkj_cse@ru.ac.bd)

**Keywords:** antimicrobial resistance, *Klebsiella pneumoniae*, colistin, *mgrB*, genomic
prediction, data leakage

## Abstract

Genomic prediction of antimicrobial susceptibility is accurate for drugs whose
resistance is conferred by acquired genes, but poor for colistin, where resistance in
*Klebsiella pneumoniae* arises largely from loss of the chromosomal regulator *mgrB*. Whether
machine learning closes that gap has not been examined at the level of an individual mechanism.

We assembled 6,652 public *K. pneumoniae* genomes with laboratory susceptibility
results, re-interpreted every minimum inhibitory concentration under one standard (EUCAST v16.1),
and added a targeted *mgrB* check addressing the assembly-break failure mode gene screens miss.
Rules were benchmarked against logistic regression and gradient boosting under random,
study-aware, lineage-aware and temporal splits, then applied, frozen, to an independent cohort
from NCBI Pathogen Detection.

Mechanism-specified *mgrB* disruption raised colistin sensitivity from 28.3% to 56.7%
for AMRFinderPlus and from 39.5% to 45.5% for Kleborate, which already calls MgrB truncations; the
breaks carried insertion-sequence signatures. This unfitted rule ranged from 39.5% to 70.6%
sensitivity across held-out collections. Models exceeded the rule under random splitting (balanced
accuracy 0.844 versus 0.759) but fell below it when sequence types were held out (0.709–0.731), a
reversal surviving correction across 54 comparisons. A sequence-type-only predictor reached just
0.554, so the advantage was not simple lineage memorisation. On 263 external isolates the frozen
rule reached 77.8% sensitivity (95% CI 67.6–85.5) at 99.5% specificity, its *mgrB* component
performing as it had internally.

The apparent superiority of machine learning for colistin is an artefact of
validation design. More fundamentally, the mechanistic composition of colistin resistance differs
substantially between public collections, so neither rules nor models transfer reliably — a limit
random and study-aware validation conceal. Genomic AMR evaluations should report a
mechanism-specified rule baseline, harmonise breakpoints, and use lineage-aware splits.

## Data Summary

All data supporting this study are publicly available, and no new sequence data were generated.

Genome assemblies and antimicrobial susceptibility records for the internal cohort were obtained
from BV-BRC (https://www.bv-brc.org). The external validation cohort derives from NCBI Pathogen
Detection (https://www.ncbi.nlm.nih.gov/pathogens), retrieved by joining the `pdbrowser.ast` and
`pdbrowser.isolates` tables; the query used is version-controlled alongside the code that consumes
it. Assemblies for that cohort were downloaded from NCBI Datasets by assembly accession, and the
complete accession list with per-file SHA-256 checksums is provided in the repository. Phenotype
interpretation follows the EUCAST Clinical Breakpoint Tables v16.1 (https://www.eucast.org).

Analysis code, pinned software environments, provenance manifests recording queries, retrieval
timestamps and checksums, architecture decision records, derived tables and figure sources are
available at https://github.com/omarfaruqe/klebsiella-amr and archived at Zenodo
(DOI 10.5281/zenodo.23668062 for the version reported here; 10.5281/zenodo.23634239 resolves to
all versions). A model card is at `docs/model-card.md`. Every reported figure and table regenerates
deterministically from the documented pipeline for a fixed seed.

## Impact Statement

Genomic prediction of antimicrobial resistance is routinely evaluated by splitting genomes at
random, a design that places near-identical isolates on both sides of the split. This study shows,
for one well-defined mechanism, what that design conceals. A targeted check for disruption of
*mgrB* — the commonest route to colistin resistance in *Klebsiella pneumoniae*, and one that
short-read assemblies systematically hide — recovers substantial sensitivity that gene screens
miss, yet the same unfitted rule varies from 39.5% to 70.6% sensitivity between public
collections. Machine-learning models outperform the rule under random splitting and fall below it
when lineages are held out. For anyone building or evaluating genomic resistance predictors, three
consequences follow: reported performance depends on validation design as much as on method;
study-aware and lineage-aware splits do not substitute for one another; and external validation on
a de-overlapped collection is achievable from public data.

## 1. Introduction

*Klebsiella pneumoniae* is among the bacterial pathogens responsible for the greatest number of
deaths attributable to antimicrobial resistance worldwide [1]. For most drug–organism
combinations in this species, resistance is conferred by acquired genes or well-characterised
target mutations, and tools such as AMRFinderPlus [2] and Kleborate [3] detect them reliably.
Colistin is the conspicuous exception, and the exception matters: colistin is a last-line agent
against carbapenem-resistant isolates [4], and resistance arises predominantly not from an
acquired gene but from loss of function of *mgrB*, a 47-codon chromosomal gene encoding a negative
regulator of the PhoPQ two-component system. Loss of MgrB de-represses PhoPQ, modifying lipid A
and reducing colistin binding [4,5]. Plasmid-borne *mcr*
genes provide an acquired route and are globally disseminated [6], but account for a minority of
colistin-resistant *K. pneumoniae*: in the external cohort of Section 3.10, 9 of 81 resistant isolates
carried one.

The commonest route to *mgrB* loss is insertion-sequence (IS) transposition into the gene [5], which
creates a specific detection problem. An IS element inserted into a 144-base-pair coding sequence
is typically longer than the repeat a short-read assembler can resolve, so the assembly breaks
inside *mgrB*. A gene-presence screen finds no acquired colistin determinant; a point-mutation
screen finds no codon change, because the alignment simply ends. The genome appears wild type at
the one locus that matters.

Two further issues complicate measuring how well genotype predicts phenotype on public data.
Deposited labels mix interpretive standards, so identical MICs can be reported susceptible under
one CLSI revision and resistant under EUCAST; colistin is especially exposed, because its
susceptibility testing is itself difficult and contested [4], and the resulting noise is easily
mistaken for genotype error. And public collections are dominated by outbreak investigations and surveys of
high-risk clones, so closely related isolates are abundant and random train/test splitting places
near-identical genomes on both sides. Yu and colleagues recently demonstrated across more than
24,000 genomes from five species that models trained this way substantially conflate lineage
markers with resistance determinants, and recommended comparison against rule-based baselines as
future work [7].

Machine learning has been applied to this problem repeatedly, including MIC panel prediction for
*K. pneumoniae* [8], and has recently been benchmarked at scale: across 78 species–antibiotic
datasets, published predictors reached F1-macro ≥0.9 in 64% of randomly split experiments but only
25% when folds were made homology-aware, while a rule-based reference method performed best on the
most divergent genomes [9]. What that literature largely lacks is a mechanism-specified rule
baseline for a drug whose determinant is known but routinely missed, which is what makes colistin
the informative case.

We take that recommendation as the design principle and ask a mechanism-level question: when a
resistance determinant is detectable but routinely missed, does supplying it to a rule outperform
letting a model find it, and does either transfer to new collections? We assemble a 6,652-genome
cohort with harmonised labels, add a targeted *mgrB* check designed around the assembly-break
failure mode, and evaluate rules and models under progressively stricter validation.

## 2. Methods

### 2.1 Cohort assembly and provenance

All *K. pneumoniae* (taxon 573) AMR phenotype records were downloaded from the BV-BRC `genome_amr`
endpoint (snapshot 2026-10-22; 1,840,793 records, 39,598 genomes). Retrieval re-run the following
day reproduced byte-identical output (SHA-256 d0e561a6…1505126). Every data pull writes a manifest
recording query, timestamp, record count and checksum.

Only records with `evidence = "Laboratory Method"` were used as labels (85,291 records, 7,276
genomes); computational predictions deposited in the same database were excluded as circular.

### 2.2 Label harmonisation

Labels follow EUCAST Clinical Breakpoint Tables v16.1 [10]. Each record carrying an MIC was
re-interpreted (S if MIC ≤ s_le, R if MIC > r_gt, else I) rather than accepting the deposited
call. Censored measurements that cannot be resolved against the breakpoint yield no label.

This recovered substantial data: records deposited without an S/I/R call were uninterpreted
measurements rather than missing data, and re-interpretation added 1,921 labelled genomes for
meropenem, 959 for ciprofloxacin, 598 for colistin and 560 for gentamicin. It also resolved most
apparent contradictions — meropenem label conflicts fell from 60 genomes to 8 — because the
majority were the same genome interpreted under different standards by different studies. Where
both a deposited and a re-interpreted call existed they agreed in 92.9–97.2% of records.

Three policies were fixed in advance and recorded as decision records in the repository:
intermediate and conflicting labels are excluded; colistin labels are restricted to broth
microdilution, the only method EUCAST accepts for colistin MICs; and study-level label sets judged
implausible on review are flagged and excluded with evidence recorded (one flag applied, Section 3.9).
Tigecycline has no EUCAST breakpoint for *K. pneumoniae*, so deposited calls were used.

### 2.3 Assembly quality control and genotyping

Assemblies were required to be *K. pneumoniae* by both species field and genome name, 5.0–6.5 Mb,
≤500 contigs, N50 ≥20 kb, CheckM [11] completeness ≥95% and contamination ≤5%. Agreement between
the species field and the genome name was required, so that internally inconsistent BV-BRC records
are excluded.

Assemblies were downloaded from BV-BRC [12] by genome identifier rather than by NCBI assembly
accession, because many QC-passing genomes carry no such accession, and were validated for format
and length. AMRFinderPlus 4.2.7 [2] (database 2026-08-07.1) was run in nucleotide mode with
`--organism Klebsiella_pneumoniae --plus`, and Kleborate 3.2.4 [3] with Kaptive 3.2.2 [13] using
the `kpsc` preset. A cross-machine verification script reproduces identical calls on macOS arm64
and Ubuntu x86_64.

### 2.4 Rule-based prediction

Per-drug rules were specified from tool output before phenotype comparison: acquired carbapenemase
for meropenem; QRDR mutation or acquired quinolone gene for ciprofloxacin; aminoglycoside-modifying
enzyme or 16S rRNA methyltransferase for gentamicin; *mcr* or *mgrB*/*pmrAB*/*phoPQ* determinants
for colistin; acquired *tet(X)* or *tmexCD-toprJ* for tigecycline.

Two corrections were applied after a 524-genome pilot, both tool-annotation artefacts rather than
biology. AMRFinderPlus assigns the intrinsic chromosomal *oqxAB* pump to tigecycline as well as
quinolones; since essentially all *K. pneumoniae* carry *oqxAB* this produced 88 of 91 pilot
tigecycline false positives, and *oqxAB* is now excluded for all drugs. Efflux-regulator variants
(*oqxR*, *ramR*, *acrR*) produced 21 of 26 pilot ciprofloxacin false positives and are reported as
a separate tier rather than folded into the core rule.

### 2.5 Targeted *mgrB* integrity check

The MgrB reference protein (WP_002911375.1, 47 aa) was searched against each assembly with tblastn [14]
(e-value ≤1e-5, hits ≥10 aa) and the best hit classified as **intact** (full-length, no residue
differences); **missense** (full-length with substitutions); **premature_stop**;
**frameshift_or_indel** (two fragments on one contig <300 bp apart); **insertion_assembled** (two
fragments >300 bp apart); **split_contigs**; **contig_break** (alignment terminating at a contig
edge); **partial**; or **absent**. All categories except intact and missense count as disrupted.
The `split_contigs`, `contig_break` and `partial` categories are the assembly signature described
in Section 1.

Residue differences are enumerated by direct comparison of aligned sequences, ignoring query
position 1, because tblastn translates alternative start codons (GTG, TTG) as V or L rather than
M. Classifying variants by percentage identity instead produces a systematic artefact at the start
codon, quantified in Section 3.3.

### 2.6 Testing the insertion-sequence interpretation

The `split_contigs`, `contig_break` and `partial` categories are inferences from assembly
structure, so we tested three predictions that follow from IS insertion but not from incidental
fragmentation. First, position: the break should coincide with *mgrB*, compared against genomes
with intact *mgrB*. Second, repetitiveness: sequence at the *mgrB*-adjacent contig end should have
copies elsewhere in the same assembly. Because contig ends are repeat-enriched in general, this was
compared against randomly chosen other contig ends in the same assemblies, truncated to the same
length as each genome's *mgrB* flank. Third, recurrence: the same sequences should appear at *mgrB*
in unrelated isolates. Flanks were taken as up to 600 bp of contig end adjacent to *mgrB* (minimum
50 bp) and compared by blastn at ≥95% identity (within genome) or ≥90% (between genomes) over ≥80%
of the flank.

### 2.7 Validating the *mgrB* rule

Because three disruption categories are inferred from assembly structure, we tested whether the
association with phenotype could reflect assembly fragmentation by comparing contig counts and N50
between calls, and by stratifying the colistin association on assembly-quality quartile.

Because the decision to treat missense variants as resistance-predicting was made after observing
their enrichment in this cohort, we re-estimated that tier by leave-one-study-out: variants were
selected on training studies only (≥3 genomes, >50% resistant) and applied to the held-out study.
The disruption tier, being mechanism-specified and fitted to nothing, was evaluated in the same
framework for comparison.

### 2.8 Attributing performance to lineage

To test directly, rather than infer, how much performance is obtainable from lineage, we built a
predictor using sequence type alone: for each test genome it returns the majority label of that
genome's sequence type among training genomes, falling back to the overall training majority when
the type is unseen. It uses no resistance determinant. We also measured, for each split scheme,
the fraction of test genomes sharing a sequence type or a BioProject with any training genome.

### 2.9 Benchmark

Features were 226 binary indicators per genome: AMRFinderPlus gene families (allele suffixes
collapsed, OXA-48-like kept distinct from other OXA enzymes), 100 point-mutation calls, 33
AMRFinderPlus functional subclasses, and two *mgrB* features. Features present in fewer than five
genomes were dropped.

Three predictors were compared: the per-drug rule (for colistin including mechanism-specified
*mgrB* disruption, but not the missense tier, for the reason given in Section 3.5); L2 logistic regression [15]
(C = 1, class-balanced); and gradient boosting (HistGradientBoostingClassifier, 200 iterations, 15
leaves, class-balanced). Hyperparameters were fixed rather than tuned, since the object of study is
validation design rather than peak performance.

Four split schemes were used: random 5-fold; study-aware 5-fold grouped by BioProject;
lineage-aware 5-fold grouped by sequence type; and temporal. Genomes with missing group metadata
were assigned singleton groups so absent metadata can never merge unrelated isolates. Metrics are
computed on pooled held-out predictions with 95% Wilson intervals [16], reported both on all test labels
and with borderline-MIC labels removed. A split is scored only if its test set contains at least 10
resistant and 10 susceptible genomes. Differences between predictors were assessed by paired
bootstrap over held-out genomes (2,000 resamples), resampling once per iteration and recomputing
both predictors on the same resampled set.

Comparing every drug, split scheme and predictor pair yields a family of 54 tests, not one, so some
intervals will exclude zero by chance. Each comparison therefore also carries a two-sided bootstrap
p-value, floored at 1/2,000 since no finer resolution exists, and a Benjamini–Hochberg [17] q-value
computed across the whole family. Benjamini–Hochberg is our primary correction: the comparisons
share predictors, genomes and folds and are therefore positively dependent, and controlling the
false discovery rate is the appropriate goal for a benchmark mapping where differences exist rather
than testing a single hypothesis. We also report which comparisons survive the far stricter
Bonferroni threshold (p < 0.00093), and we treat only those as established beyond doubt.

### 2.10 External validation

Leave-one-study-out within BV-BRC holds out a study but not an aggregator: every such study was
curated the same way and was available to us when the rule was specified. To test the rule on
genomes whose phenotypes we had no access to, we assembled a cohort from NCBI Pathogen Detection,
joining submitted antimicrobial susceptibility results (`pdbrowser.ast`) to isolate records
(`pdbrowser.isolates`) for the assembly accession, since neither table carries both the MIC and
the accession. All *K. pneumoniae* colistin results with an assembly accession were
retrieved.

Because BV-BRC ingests from NCBI, overlap was removed before anything else: isolates matching our
cohort by BioSample accession, or by assembly accession ignoring the version suffix, were excluded.
Labels were re-derived from the submitted MIC under the same pinned EUCAST table used internally,
never taken from the submitted phenotype; isolates without a usable MIC were retained as a separate
stratum carrying the submitter's own interpretation. Assemblies were downloaded from NCBI Datasets
and verified against each archive's md5. AMRFinderPlus and Kleborate were run through the identical
pipeline code as the internal cohort, and the rules were imported unchanged from the analysis
module. Nothing was refitted, reselected or thresholded for this cohort, and the evaluation exposes
no parameter that could be.

### 2.11 Reproducibility

All analyses are deterministic for a fixed seed (20261024). Code, pinned environments, manifests,
decision records, reports and figures are at https://github.com/omarfaruqe/klebsiella-amr.

## 3. Results

### 3.1 Cohort

Of 6,817 labelled genomes, 6,652 passed quality control. Requiring agreement between the species
field and the genome name removed four genomes with internally inconsistent BV-BRC records, one
independently confirmed as *Citrobacter freundii*. Of the genomes that passed, 63% (4,190/6,652)
carry no NCBI assembly accession, which is why assemblies were retrieved by BV-BRC identifier. All
6,652 completed both genotyping tools (8.8 h wall time, ~140 CPU-hours, zero failures), and the
targeted *mgrB* search took 2 minutes for 6,652 assemblies on 16 threads.

The cohort comprised 6,652 QC-passing genomes with harmonised labels: meropenem 5,541 (1,854 R /
3,687 S), ciprofloxacin 4,640 (3,441 / 1,199), gentamicin 4,760 (2,113 / 2,647), colistin 1,528
broth-microdilution (367 / 1,161), tigecycline 948 (261 / 687). It is structurally non-random, as
expected: two BioProjects contribute 44% of genomes, the ten commonest sequence types account for
52%, and collection year is missing for 34% (2,266/6,652).

### 3.2 Known mechanisms explain four drugs well and colistin poorly

**Table 1.** Sensitivity and specificity of each gene screen per drug, on harmonised labels.

| Drug | AMRFinderPlus sensitivity | specificity | Kleborate sensitivity | specificity |
|---|---|---|---|---|
| Ciprofloxacin | 98.8% (98.4–99.1) | 88.2% | 98.7% | 88.4% |
| Meropenem | 96.7% (95.7–97.4) | 91.2% | 96.5% | 91.3% |
| Gentamicin | 85.1% (83.6–86.6) | 94.9% | 81.4% | 95.8% |
| **Colistin** | **28.3% (24.0–33.2)** | 97.8% | **39.5% (34.6–44.6)** | 96.3% |
| **Tigecycline** | **2.7% (1.3–5.4)** | 100.0% | **0.0%** | 100.0% |

The tools agree on 97.8–99.9% of genome–drug pairs for every drug except colistin (87.8%), the one
drug where they differ in what they look for: Kleborate's higher sensitivity derives from MgrB
truncation calls that AMRFinderPlus's point-mutation screen does not make. For colistin, 263 of 367
resistant genomes carried no AMRFinderPlus determinant at all.

### 3.3 A targeted *mgrB* check recovers part of the gap

Across all 6,652 genomes *mgrB* was intact in 6,056, missense in 129, and disrupted in 467, the
last dominated by the assembly signatures predicted in Section 1: split across contigs (162), partial
(153), premature stop (73), absent (34), frameshift or indel (27), contig break (13) and assembled
internal insertion (5).

**Table 2.** *mgrB* status against colistin phenotype among the 1,528 genomes with
broth-microdilution labels (24% resistant overall).

| *mgrB* status | Resistant | Susceptible |
|---|---|---|
| Disrupted | 153 | 41 |
| Missense | 70 | 13 |
| Intact | 144 | 1,107 |

Classification required comparing aligned residues directly. An initial implementation using
percentage identity called 6,048 genomes missense at exactly 97.9% identity, a single artefactual
mismatch arising because tblastn translates alternative start codons as V or L rather than M;
ignoring query position 1 removes it.

Adding mechanism-specified disruption raises sensitivity for both gene screens while costing
little specificity (Figure 1; exact values with intervals in Supplementary Table S1):

![Figure 1](../../reports/figures/colistin_mgrb_light.png)

**Figure 1.** Colistin sensitivity and specificity for each rule, with and without the targeted
*mgrB* disruption check. Dots are point estimates, bars 95% Wilson intervals; broth-microdilution
labels (n = 1,528; 367 resistant). The specificity axis begins at 85%.

The increment is large against AMRFinderPlus (+28.4 points) and modest against Kleborate (+6.0
points), which already detects MgrB truncation by sequence coverage. Kleborate is the appropriate
comparator for the state of the art, and the honest gain is therefore six points. Notably, the
*mgrB* check used alone, with no gene screen, slightly exceeds Kleborate alone.

**The association is not an assembly artefact.** Median contig counts were similar between calls
(disrupted 116, intact 104, missense 108), and stratifying on assembly quality left the association
intact:

**Table 3.** Colistin resistance by *mgrB* status within quartiles of assembly fragmentation.

| Contig quartile | *mgrB* intact, % resistant | *mgrB* disrupted, % resistant |
|---|---|---|
| Q1 (fewest contigs) | 4% (n=337) | 81% (n=32) |
| Q2 | 8% (n=315) | 69% (n=48) |
| Q3 | 17% (n=292) | 83% (n=65) |
| Q4 (most contigs) | 18% (n=307) | 82% (n=49) |

The gap is widest, not narrowest, among the least fragmented assemblies. Disruption calls inferred
from assembly structure behaved similarly to sequence-level calls (76% versus 86% resistant,
against 12% for intact *mgrB*), indicating that the assembly-break categories carry genuine signal.

### 3.4 The contig breaks carry insertion-sequence signatures

Among the 328 genomes whose *mgrB* disruption was inferred from assembly structure, the break
coincided with the gene: median distance from the *mgrB* fragment to the contig end was 79 bp, with
63% within 100 bp and 87% within 1 kb. In 328 matched genomes with intact *mgrB*, the gene lay a
median of 80,470 bp from any contig end, and 0% were within 100 bp. Whatever terminates these
assemblies is located at *mgrB* specifically.

The adjacent sequence behaved like a mobile element. Of 255 genomes with at least 50 bp of flank,
82% had a high-identity copy elsewhere in their own assembly and 72% had three or more (median 5),
against 66% and 35% (median 2) for length-matched random contig ends from the same assemblies. And
249 of 255 flanks (98%) matched a flank from a different genome, which idiosyncratic assembly
artefacts would not.

Together these support the interpretation that the assembly-inferred categories represent IS
insertions in *mgrB*. They do not identify which element is involved, nor confirm insertion in any
individual genome; both require an IS reference database or long-read assemblies of these isolates.

### 3.5 A plausible missense tier does not survive out-of-sample selection

Missense variants were enriched among resistant isolates — M27K (25 R / 3 S), T21N (9 / 0), N42Y
(7 / 0), K43I (7 / 0), C28W (4 / 1), C39G (3 / 0) — and C28 and C39 are the two cysteines of this
small lipoprotein, so a structural rationale was readily available. Counting missense as
resistance-predicting raised in-sample sensitivity of the *mgrB* check from 41.7% to 60.8%, an
apparent gain of 19.1 points.

Under leave-one-study-out selection the tier added **0.9 points** (47.0% to 47.9%). The apparent
benefit was almost entirely an artefact of selecting variants on the same data used to evaluate
them. We therefore exclude the missense tier from all reported rules and from the benchmark.

We report this negative result because the failure mode is instructive: the variants were
individually plausible, the structural rationale was real, and the in-sample improvement was large
and would have survived cross-validation that did not group by study.

### 3.6 Mechanism composition differs between collections

The disruption tier is fitted to nothing, so its leave-one-study-out performance measures transfer
rather than overfitting. It varied substantially (Figure 2; per-study values in Supplementary
Table S2):

![Figure 2](../../reports/figures/colistin_transfer_light.png)

**Figure 2.** Sensitivity of the *mgrB* disruption rule in each held-out study. The rule is
mechanism-specified and fitted to nothing, so the spread reflects differing mechanistic
composition between collections rather than overfitting.

In the largest collection, *mgrB* disruption accounts for under 40% of colistin resistance; in
another it accounts for over 70%. Since nothing is fitted, this is not an overfitting artefact but
a statement about the isolates: the mechanistic composition of colistin resistance genuinely
differs between collections, presumably reflecting differences in clonal background, geography and
selection history. Any fixed rule, and any model trained on one collection, inherits this limit.

### 3.7 The ordering of rules and models inverts with validation stringency

Balanced accuracy on held-out genomes reverses the ordering of rule and models as grouping
tightens (Figure 3; values in Supplementary Table S3):

![Figure 3](../../reports/figures/benchmark_splits_light.png)

**Figure 3.** Balanced accuracy of the per-drug rule, logistic regression and gradient boosting
under random, study-aware and lineage-aware splits. Rules are deterministic and therefore flat
across splits by construction; the comparison of interest is whether model lines converge on, or
fall below, the rule line as grouping tightens.

Paired bootstrap differences (rule minus model, 95% interval) confirm the reversal: under random
splitting the models are ahead by 0.085 (−0.113 to −0.059) and 0.083 (−0.111 to −0.054); under
study-aware splitting by 0.063 (−0.088 to −0.037) and 0.051 (−0.078 to −0.026); under lineage-aware
splitting the rule is ahead by 0.028 (+0.002 to +0.052) and 0.050 (+0.026 to +0.075). All six
intervals exclude zero, but they are drawn from a family of 54 comparisons and do not survive
correction equally. The rule's margin over gradient boosting under lineage-aware splitting holds
at p < 0.0005 (Benjamini–Hochberg q = 0.001) and survives Bonferroni; its margin over logistic
regression does not (p = 0.026, q = 0.045, Bonferroni threshold 0.00093). The reversal against
gradient boosting is therefore established, while the reversal against logistic regression should
be read as comparable performance rather than superiority. Across the family, 33 of 54 intervals
excluded zero uncorrected, 31 survived Benjamini–Hochberg and 24 survived Bonferroni; the two that
fell away under Benjamini–Hochberg were small ciprofloxacin study-aware differences (−0.008 and
−0.005) on which no claim here rests.

The models had access to the same *mgrB* features as the rule, and lost ground only when
sequence types were held out. The margin by which the rule leads under lineage-aware splitting is
small, and the practical claim is not that rules are superior but that the models' apparent
superiority is an artefact of validation design.

Across drugs, random-split optimism — balanced accuracy under random splitting minus the lower of
the two grouped splits — was 0.007–0.042 for meropenem, ciprofloxacin and gentamicin, but
0.113–0.133 for colistin and 0.136–0.149 for tigecycline. For tigecycline, models exceeded the rule
under every split including lineage-aware (differences −0.203 and −0.107), though at 0.62–0.72 all
predictors remain weak. Temporal validation was interpretable only for meropenem, ciprofloxacin and
gentamicin; colistin (34 test genomes, one resistant) and tigecycline (28, none) were not evaluable.

### 3.8 The leaked signal is finer-grained than sequence type

A predictor using sequence type alone, with no resistance determinant of any kind, reached 0.737
balanced accuracy for meropenem, 0.728 for ciprofloxacin and 0.713 for gentamicin under random
splitting, confirming that lineage carries substantial predictive signal in this cohort. For those
three drugs the models nonetheless reached ~0.95 and did not degrade when sequence types were held
out, so their performance rests on resistance determinants rather than lineage.

Colistin behaved differently, and not as we had expected. Sequence type alone reached only 0.554
under random splitting, despite 91% of test genomes sharing a type with training data. Colistin
resistance is therefore poorly predicted by sequence type: both phenotypes occur within the same
types. Yet the models still fell from 0.844 to 0.731 when types were held out. The signal they
lose is thus not sequence-type membership as such, but something finer-grained that co-varies with
it — plausibly accessory-genome similarity among closely related isolates within a type. We report
this because it corrects a natural but incorrect reading of the split results: the models were not
simply memorising lineage labels.

The two grouping schemes also proved not to be substitutes. Under study-aware splitting, 74–88% of
test genomes still shared a sequence type with training data; under lineage-aware splitting, 87–95%
still shared a BioProject. Each scheme leaves the other dimension almost entirely unconstrained, so
grouping on both would be stricter than either alone. Temporal splitting, often presented as the
realistic option, constrained neither: 73–78% of test genomes shared a sequence type and 76–82% a
BioProject with training data.

### 3.9 Much apparent genotype error is label artefact

**Meropenem.** Of 2,118 carbapenemase carriers, 326 tested susceptible, strongly enzyme-dependent:
52.9% of OXA-48-like carriers (182/344) and 75.0% of VIM carriers (33/44) versus 5.9% of KPC
(75/1,278) and 10.0% of NDM (28/280). Of the 240 with a recorded MIC, 191 (80%) lay at 1–2 mg/L and
224 (93%) exceeded the EUCAST carbapenemase screening cut-off of 0.125 mg/L. The genotype call is
correct in most cases; the clinical S category does not reflect carbapenemase carriage for these
families. Of 62 resistant non-carriers, 51 carried an ESBL or AmpC and 29 of those also had an
OmpK35/36 defect — but 624 susceptible non-carriers had the same combination, making it a risk
factor rather than a predictor.

**Gentamicin.** Of 299 unexplained resistant isolates with an MIC, 213 (71%) were at exactly 4
mg/L, one doubling dilution above the EUCAST breakpoint, whereas isolates with a detected mechanism
were mostly ≥8 mg/L. Among CLSI-tested records in this borderline group, 160 of 162 had been
reported susceptible by the originating laboratory. Excluding borderline labels raises gentamicin
rule balanced accuracy from 0.906 to 0.959, above every model under every split.

**A flagged study.** Of 135 gentamicin false positives, 65 carried *rmtB*; 63 came from a single
BioProject at MIC 1 mg/L, implausible for a 16S rRNA methyltransferase. That study's gentamicin
labels are flagged and excluded, with evidence recorded.

### 3.10 The rule transfers to an external collection, which is composed differently

Of 355 NCBI colistin results, 84 (23.7%) were already in our cohort and were excluded, leaving 271
held out. Kleborate's species check rejected six that NCBI's taxgroup labels *K. pneumoniae*: all
six are *Klebsiella aerogenes*, four of them complete single-contig genomes, so this is
misassignment rather than poor assembly. AMRFinderPlus processed all six without complaint. A
cohort built from NCBI taxgroups without an independent species check therefore inherits a 2.2%
contamination rate from a source that does not flag it.

Of the 270 held-out isolates with a usable label, six were excluded as *K. aerogenes* and one
intermediate isolate cannot inform sensitivity, leaving 263 (81 resistant, 182 susceptible). The
frozen rule reached 77.8% sensitivity (95% CI 67.6–85.5) at 99.5% specificity, against 56.7%
(51.6–61.7) internally; restricting to the 251 isolates with MIC-derived labels gave 78.3%
(67.2–86.4). The rule therefore transfers, and the external estimate is higher than the internal
one with non-overlapping intervals.

The reason is composition, not rule quality. The *mgrB* component performed almost identically in
the two collections (45.7% versus 41.7% alone), while AMRFinderPlus alone rose from 28.3% to
50.6%. Among the 81 resistant isolates, 26 were explained by an AMRFinderPlus determinant alone,
22 by *mgrB* disruption alone, 15 by both, and 18 (22.2%) by neither. The determinants were
predominantly chromosomal point mutations in AMRFinderPlus's catalogue rather than acquired genes:
*pmrB* T157P in 10 isolates, *mgrB* nonsense and frameshift alleles in 12, *crrB* substitutions in
4, and single substitutions in *pmrA* (two isolates), *pmrB* S85R, *phoQ* and *lapB*, with only
9 isolates carrying any *mcr* gene.
Correspondingly, the largest disrupted *mgrB* category here was premature stop (11 genomes), where
internally it was the assembly-structural categories of Section 3.4. Point substitutions are what
AMRFinderPlus detects and insertion-driven disruption is what it misses, so a collection weighted
toward the former leaves the off-the-shelf tool less to miss and the targeted check less to add.
The check still added 27.2 points over AMRFinderPlus alone, explaining 22 resistant isolates it
missed, so the contribution of Section 3.3 replicates out of sample.

Two secondary results follow. The missense tier added 1.2 points externally (77.8% to 79.0%, on
three isolates) against 10.6 points in-sample, independently confirming the decision in Section 3.5 to
discard it. And 22.2% of external resistant isolates remain unexplained by any determinant or
*mgrB* disruption, close to the internal picture, so the unexplained fraction is a property of the
mechanism catalogue rather than of one collection.

## 4. Discussion

Three findings here bear on how genomic AMR prediction should be evaluated, and they point in the
same direction: the limiting factor is not model capacity but the relationship between the training
collection and the target population.

The *mgrB* result is the smallest of the three but the most concrete. A two-minute check of a single
47-codon locus recovers 28 points of colistin sensitivity over AMRFinderPlus — though only six over
Kleborate, which already detects truncation, and we think the six-point figure is the one that
should be quoted. The mechanism is instructive even where the increment is modest: the categories
carrying most signal are not sequence features but assembly artefacts, contig boundaries that are
diagnostic because an IS element has landed in a very short gene. Three independent lines of
evidence support that reading (Section 3.4) — the breaks sit a median of 79 bp from *mgrB* against 80 kb
in intact genomes, the adjacent sequence is multi-copy within its own assembly beyond what
length-matched contig ends show, and 98% of these sequences recur in unrelated isolates. Screens that ask
"which genes are present?" cannot see this; a screen that asks "is this locus whole?" can. That
generalises to any resistance determinant whose loss, rather than presence, confers the phenotype.

The missense tier is a cautionary result we think worth reporting in full. The variants were
individually plausible, several fell on structurally critical cysteines, and the in-sample gain was
19 points. Out-of-sample it was one. Nothing about the biology was obviously wrong; what was wrong
was that the variants were chosen by looking at the same isolates used to score them. We suspect
this failure mode is common in curated resistance-determinant catalogues assembled from
observational association, and it is invisible to any validation scheme that does not group by
study.

The heterogeneity result (Section 3.6) is in our view the most important, and the least avoidable. With
nothing fitted, the same mechanism-specified rule achieved 39.5% sensitivity in one collection and
70.6% in another. Colistin resistance is simply constituted differently in different isolate sets.
This bounds what any method can achieve: a model trained where *mgrB* explains 70% of resistance
will underperform where it explains 40%, however well regularised. It also bears on the crossover
in Section 3.7. We initially read that crossover as the models memorising lineage, but the sequence-type
baseline in Section 3.8 shows otherwise: for colistin, type alone predicts barely better than chance.
What the models lose when types are held out is therefore a finer-grained resemblance between
related isolates, not type identity — a resource that likewise disappears on genuinely novel
lineages, but one that a reader should not mistake for simple label memorisation.

We want to be careful not to overstate the comparison with machine learning. The rule leads under
lineage-aware splitting by 0.028 and 0.050 balanced-accuracy points; both intervals exclude zero,
but after correcting across the 54-comparison family only the second survives a Bonferroni
threshold, and the first is better described as parity than as an advantage. For tigecycline the
models win under every split. The defensible
claim is not that rules beat models, but that the models' advantage under conventional validation
is attributable to information that will not be present in deployment. Yu and colleagues
demonstrated this across species using clade-based splits [7], and Hu and colleagues reported the
same pattern across 78 species–antibiotic datasets under homology-aware folds [9]; the
contribution here is to show it at the level of a single identifiable mechanism, where one can point to exactly what the rule uses
and what the model substitutes for it.

The label findings carry comparable methodological weight. Roughly a third of our usable labels did
not exist as S/I/R calls but as uninterpreted MICs. Conversely, a substantial fraction of apparent
genotype error was *created* by harmonisation: the gentamicin borderline effect moves rule balanced
accuracy by five points and depends entirely on whether one applies EUCAST or the CLSI revision
under which those isolates were originally reported. A benchmark that silently selects or mixes
standards reports a number whose provenance a reader cannot assess.

### 4.1 Limitations

Public surveillance data over-represent outbreaks and high-risk clones; two BioProjects supply 44%
of our genomes, and prevalences here are not population estimates. Residual label error remains
beyond the one flagged study. Colistin and tigecycline have modest numbers of resistant isolates
(367 and 261). Collection year is missing for 34% of genomes and dated genomes concentrate in
2010–2014, so temporal validation was possible for three drugs only. Features are presence/absence
calls from assemblies: copy number, expression, plasmid context and promoter-region insertions are
not represented, and the 144 colistin-resistant isolates with intact *mgrB* coding sequence likely
include promoter insertions we cannot see. Hyperparameters were fixed rather than tuned, so the
models are reasonable defaults rather than optimised systems. Finally, while the positional, repeat-content and
cross-genome evidence in Section 3.4 collectively supports the IS interpretation, we did not identify
which element is involved or confirm insertion in any individual genome; that requires an IS
reference database or long-read assemblies of these isolates, neither of which we generated. The
mechanistic label on those categories therefore remains an inference, albeit a well-supported one.
The external cohort did not advance this question: only one of its genomes fell into the
assembled-insertion category, so the element remains unidentified there too.

The external cohort has limitations of its own. It is small (263 isolates, 81 resistant), so its
sensitivity interval spans 18 points. Its isolates are those whose submitters chose to deposit
colistin susceptibility results, which over-represents collections curated for characterised
mechanisms; its higher sensitivity should therefore be read as an estimate on a mechanism-weighted
collection rather than a field estimate, and the internal and external figures bracket rather than
supersede one another. Thirteen isolates retain the submitter's own breakpoint interpretation,
though excluding them changed sensitivity by 0.5 points. For colistin specifically the CLSI and
EUCAST breakpoints effectively coincide, so harmonisation altered one label among the 355
records retrieved; it
remains necessary for the other drugs but did little work here.

Most importantly, this study measures discrimination, not clinical usefulness, and the two come
apart. Sensitivity and specificity are properties of the rule; predictive values also depend on how
common resistance is in the population tested, and 24% of our colistin-tested genomes are
resistant — a surveillance figure far above most clinical settings. Applying the internal rule's
56.7% sensitivity and 95.2% specificity at a 2% prevalence gives a positive predictive value of
19.4%: four in five genotype-positive isolates would not be resistant. At our own cohort's
prevalence the complementary error dominates, with a negative predictive value of 87.4%, meaning
12.6% of genotype-negative isolates are resistant. For a last-line agent these two errors are not
interchangeable — withholding colistin unnecessarily and administering it ineffectively carry
different costs — and we performed no decision analysis weighing them. Nothing here establishes
that genotypic prediction should inform colistin therapy, and we make no such claim: the rule's
role in this paper is as a mechanism-specified baseline against which to measure validation
design. Predictive values across a range of prevalences are tabulated in
`reports/review/predictive_values.md`.

### 4.2 Recommendations

1. **Report a mechanism-specified rule baseline per drug.** It is the only way to tell whether a
   model has learned a mechanism or a lineage.
2. **Select determinants out-of-sample.** A plausible variant enriched in a cohort is not a
   validated determinant; leave-one-study-out selection cost our missense tier 18 of its 19 points.
3. **State and harmonise the breakpoint standard, and report sensitivity to borderline labels.**
4. **Use lineage-aware splits and report the gap from random splitting.** The gap was under one
   point where rules work and 11–15 points where they do not.
5. **Group on study and lineage together, not either alone.** In our cohort, study-aware splits
   still shared a sequence type with training data for 74–88% of test genomes, and lineage-aware
   splits still shared a BioProject for 87–95%. Temporal splits constrained neither (73–78% and
   76–82%). Each single scheme leaves the other dimension almost entirely unconstrained.
6. **Report a lineage-only baseline.** It costs nothing and distinguishes a model that has learned
   a mechanism from one that has learned who its neighbours are.
7. **Report per-collection performance, not only pooled.** Pooled estimates concealed a two-fold
   range in our data, and a 21-point difference between our cohort and an external one.
8. **Verify species independently of the source's own labels.** Six of 270 isolates in NCBI's
   *K. pneumoniae* taxgroup were *K. aerogenes*, including complete genomes; Kleborate's species
   check caught them and AMRFinderPlus did not.

## 5. Conflicts of interest

The authors declare that there are no conflicts of interest.

## 6. Funding

This research received no external funding.

## 7. Author contributions

Following CRediT. **Md Omar Faruqe:** conceptualisation, methodology, software, formal analysis,
data curation, investigation, validation, visualisation, writing — original draft. **Md Humayun
Kabir:** conceptualisation, supervision, project administration, writing — review and editing.

## 8. Acknowledgements

Generative artificial intelligence (Claude, Anthropic) was used throughout this work for drafting
and editing support, including text generation, manuscript organisation and revision of tone, and
for assistance in writing the analysis code. Candidate interpretations suggested during the
analysis were treated as hypotheses and tested empirically; several were revised or discarded where
the data did not support them. All study design decisions, verification and conclusions are the
authors'. Every reported value was checked by the authors against the underlying data, the complete
analysis code and provenance manifests are openly available and archived (see Data Summary), and
the authors accept full responsibility for the content and integrity of the work.

## 9. References

All references were verified against the publisher record during preparation (references 1–6 on 2026-10-24, references 7–9 on 2026-10-25).

1. Antimicrobial Resistance Collaborators. Global burden of bacterial antimicrobial resistance in
   2019: a systematic analysis. *The Lancet*. 2022;399(10325):629–655.
   https://doi.org/10.1016/S0140-6736(21)02724-0
2. Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, Hoffmann M,
   Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. AMRFinderPlus and the Reference Gene
   Catalog facilitate examination of the genomic links among antimicrobial resistance, stress
   response, and virulence. *Scientific Reports*. 2021;11(1). https://doi.org/10.1038/s41598-021-91456-0
3. Lam MMC, Wick RR, Watts SC, Cerdeira LT, Wyres KL, Holt KE. A genomic surveillance framework and
   genotyping tool for *Klebsiella pneumoniae* and its related species complex. *Nature
   Communications*. 2021;12:4188. https://doi.org/10.1038/s41467-021-24448-3
4. Poirel L, Jayol A, Nordmann P. Polymyxins: antibacterial activity, susceptibility testing, and
   resistance mechanisms encoded by plasmids or chromosomes. *Clinical Microbiology Reviews*.
   2017;30(2):557–596. https://doi.org/10.1128/CMR.00064-16
5. Cannatelli A, Giani T, D'Andrea MM, Di Pilato V, Arena F, Conte V, Tryfinopoulou K,
   Vatopoulos A, Rossolini GM. MgrB inactivation is a common mechanism of colistin resistance in
   KPC-producing *Klebsiella pneumoniae* of clinical origin. *Antimicrobial Agents and
   Chemotherapy*. 2014;58:5696–5703. https://doi.org/10.1128/AAC.03110-14
6. Liu YY, Wang Y, Walsh TR, Yi LX, Zhang R, Spencer J, Doi Y, Tian G, Dong B, Huang X, Yu LF,
   Gu D, Ren H, Chen X, Lv L, He D, Zhou H, Liang Z, Liu JH, Shen J. Emergence of
   plasmid-mediated colistin resistance mechanism MCR-1 in animals and human beings in China: a
   microbiological and molecular biological study. *The Lancet Infectious Diseases*.
   2016;16(2):161–168. https://doi.org/10.1016/S1473-3099(15)00424-7
7. Yu Y, Wheeler NE, Barquist L. Biased sampling driven by bacterial population structure confounds
   machine learning prediction of antimicrobial resistance. *PLOS Biology*.
   2025;23(12):e3003539. https://doi.org/10.1371/journal.pbio.3003539
8. Nguyen M, Brettin T, Long SW, Musser JM, Olsen RJ, Olson R, Shukla M, Stevens RL, Xia F,
   Yoo H, Davis JJ. Developing an in silico minimum inhibitory concentration panel test for
   *Klebsiella pneumoniae*. *Scientific Reports*. 2018;8:421. https://doi.org/10.1038/s41598-017-18972-w
9. Hu K, Meyer F, Deng Z-L, Asgari E, Kuo T-H, Münch PC, McHardy AC. Assessing computational
   predictions of antimicrobial resistance phenotypes from microbial genomes. *Briefings in
   Bioinformatics*. 2024;25(3):bbae206. https://doi.org/10.1093/bib/bbae206
10. The European Committee on Antimicrobial Susceptibility Testing. Breakpoint tables for
    interpretation of MICs and zone diameters. Version 16.1, valid from 2026-06-24.
    https://www.eucast.org
11. Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. CheckM: assessing the quality of
    microbial genomes recovered from isolates, single cells, and metagenomes. *Genome Research*.
    2015;25(7):1043–1055. https://doi.org/10.1101/gr.186072.114
12. Olson RD, Assaf R, Brettin T, et al. Introducing the Bacterial and Viral Bioinformatics Resource
    Center (BV-BRC): a resource combining PATRIC, IRD and ViPR. *Nucleic Acids Research*.
    2023;51(D1):D678–D689. https://doi.org/10.1093/nar/gkac1003
13. Lam MMC, Wick RR, Judd LM, Holt KE, Wyres KL. Kaptive 2.0: updated capsule and
    lipopolysaccharide locus typing for the *Klebsiella pneumoniae* species complex. *Microbial
    Genomics*. 2022;8:000800. https://doi.org/10.1099/mgen.0.000800
14. Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+:
    architecture and applications. *BMC Bioinformatics*. 2009;10:421.
    https://doi.org/10.1186/1471-2105-10-421
15. Pedregosa F, Varoquaux G, Gramfort A, et al. Scikit-learn: machine learning in Python. *Journal
    of Machine Learning Research*. 2011;12:2825–2830.
    https://jmlr.org/papers/v12/pedregosa11a.html
16. Wilson EB. Probable inference, the law of succession, and statistical inference. *Journal of
    the American Statistical Association*. 1927;22(158):209–212.
    https://doi.org/10.1080/01621459.1927.10502953
17. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful
    approach to multiple testing. *Journal of the Royal Statistical Society: Series B
    (Statistical Methodology)*. 1995;57(1):289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
