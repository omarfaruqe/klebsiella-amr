# A targeted *mgrB* integrity check doubles genomic detection of colistin resistance in *Klebsiella pneumoniae* and outperforms machine learning under lineage-aware validation

Md Omar Faruqe

Montréal, Québec, Canada. Correspondence: faruqe.cse@gmail.com

Code and data: https://github.com/omarfaruqe/klebsiella-amr

---

## Abstract

**Background.** Genomic prediction of antimicrobial susceptibility works well for drugs whose
resistance is driven by acquired genes, but poorly for colistin, where resistance arises mainly
from loss of the chromosomal regulator *mgrB*. Insertion-sequence (IS) disruption of *mgrB*
frequently breaks short-read assemblies inside the gene, so gene-detection and point-mutation
screens report nothing.

**Methods.** We assembled 6,652 public *K. pneumoniae* genomes from BV-BRC that passed assembly
quality control and carried laboratory antimicrobial susceptibility testing (AST) results. All
minimum inhibitory concentrations (MICs) were re-interpreted under a single standard, EUCAST
v16.1, which recovered 30–50% more usable labels per drug than the deposited S/I/R calls and
resolved most apparent label conflicts. Genotypes were called with AMRFinderPlus 4.2.7 and
Kleborate 3.2.4. We added a targeted *mgrB* integrity check: tblastn of the reference MgrB
protein against each assembly, classifying the locus as intact, missense, or disrupted
(premature stop, frameshift, internal insertion, split across contigs, truncated at a contig
edge, or absent). We then benchmarked rule-based prediction against L2 logistic regression and
gradient boosting on 226 binary genomic features, under random, study-aware (BioProject-grouped),
lineage-aware (sequence-type-grouped) and temporal splits.

**Results.** Known mechanisms explained resistance well for meropenem (sensitivity 96.7%),
ciprofloxacin (98.8%) and gentamicin (85.1%), but poorly for colistin (28.3%) and tigecycline
(2.7%). Among 1,528 genomes with broth-microdilution colistin labels (367 resistant), *mgrB* was
disrupted in 153 resistant versus 41 susceptible genomes, and carried a missense variant in 70
versus 13. Adding the *mgrB* check to AMRFinderPlus raised colistin sensitivity from 28.3%
(95% CI 24.0–33.2) to 67.3% (62.3–71.9), at a cost of 3.5 percentage points of specificity
(97.8% to 94.3%). The *mgrB* check alone reached 60.8% sensitivity, exceeding either standard
tool. Under lineage-aware validation, machine-learning sensitivity for colistin collapsed
(logistic regression 58.3%, gradient boosting 50.4%) while the mechanism-based rule was unchanged
at 67.3%; balanced accuracy was 0.808 for the rule versus 0.731 and 0.709 for the models.
Random splits overstated model balanced accuracy by 11–15 points for colistin and tigecycline,
but under one point for drugs where rules already work. Several apparent genotype errors were
label artefacts: 93% of carbapenemase-carrying but meropenem-susceptible isolates with an MIC
exceeded the EUCAST carbapenemase screening cut-off, and 71% of unexplained gentamicin-resistant
isolates sat exactly one doubling dilution above the breakpoint, where 160 of 162 CLSI-tested
records had been reported susceptible by the originating laboratory.

**Conclusions.** A small, mechanism-specific assembly check recovers more colistin resistance
than either standard genotyping tool or either machine-learning model, and unlike the models it
does not degrade when tested on unseen lineages. For resistance mechanisms that break assemblies,
targeted locus-integrity checks are a better investment than larger models. Evaluations of
genomic AMR prediction should report a per-drug rule baseline, harmonise breakpoints explicitly,
and use grouped validation splits.

**Keywords:** antimicrobial resistance, *Klebsiella pneumoniae*, colistin, *mgrB*, insertion
sequence, genomic prediction, data leakage, EUCAST

---

## 1. Introduction

Whole-genome sequencing is increasingly used to infer antimicrobial susceptibility. For most
drug–bug combinations in *Klebsiella pneumoniae*, resistance is conferred by acquired genes or by
well-characterised target mutations, and tools such as AMRFinderPlus and Kleborate detect them
reliably. Colistin is the conspicuous exception. Colistin is a last-line agent against
carbapenem-resistant *K. pneumoniae*, and resistance arises predominantly not from an acquired
gene but from loss of function of *mgrB*, a 47-codon chromosomal gene encoding a negative
regulator of the PhoPQ two-component system. Loss of MgrB de-represses PhoPQ, triggering lipid A
modification and reducing colistin binding.

The commonest route to *mgrB* loss is insertion-sequence transposition into the gene. This creates
a specific and under-appreciated detection problem: an IS element inserted into a 144-base-pair
coding sequence is typically longer than the repeat the assembler can resolve, so the assembler cannot bridge
it. The assembly breaks inside *mgrB*. A gene-presence screen finds no acquired colistin gene; a
point-mutation screen finds no codon change, because the reference alignment simply ends. The
genome looks wild type at the one locus that matters.

Two further issues complicate any attempt to measure how well genotype predicts phenotype on
public data. First, deposited susceptibility labels mix interpretive standards: the same MIC can
be reported susceptible under one CLSI breakpoint revision and resistant under EUCAST, so
unharmonised cohorts contain systematic label noise that is easily mistaken for genotype error.
Second, public genome collections are dominated by outbreak investigations and surveys of
high-risk clones, so closely related isolates are abundant and a random train/test split places
near-identical genomes on both sides. Yu and colleagues recently demonstrated, across more than
24,000 genomes from five species, that machine-learning models trained this way substantially
conflate lineage markers with resistance determinants, and that performance falls when test
isolates come from held-out clades [1]. They recommended, as future work, comparison against
simple rule-based baselines.

Here we take that recommendation as the design principle. We ask not whether leakage exists, but
what to do about it: whether a small, mechanism-specific correction to the rules outperforms
larger models once validation is honest. We assemble a 6,652-genome *K. pneumoniae* cohort with
harmonised EUCAST v16.1 labels, quantify where known mechanisms fail, add a targeted *mgrB*
integrity check designed around the assembly-break failure mode, and benchmark rules against
models under random, study-aware, lineage-aware and temporal splits.

---

## 2. Methods

### 2.1 Cohort assembly and provenance

All *K. pneumoniae* (taxon 573) AMR phenotype records were downloaded from the BV-BRC
`genome_amr` endpoint (snapshot 2026-10-22; 1,840,793 records, 39,598 genomes). Retrieval was
re-run the following day and reproduced byte-identical output (SHA-256 d0e561a6…1505126),
confirming deterministic ingestion. Every data pull writes a manifest recording the query,
retrieval timestamp, record count and file checksum.

Only records with `evidence = "Laboratory Method"` were used as labels (85,291 records across
7,276 genomes); computational predictions deposited in the same database were excluded, as using
them as ground truth would be circular.

### 2.2 Label harmonisation

Susceptibility labels follow a single interpretive standard: EUCAST Clinical Breakpoint Tables
v16.1, valid from 2026-06-24 [2]. For each record carrying an MIC, the MIC was re-interpreted
(S if MIC ≤ s_le; R if MIC > r_gt; otherwise I), rather than accepting the deposited call.
Censored measurements that cannot be resolved against the breakpoint (for example "≤4" against
an S breakpoint of 2) yield no label.

This recovered substantial data. Records deposited without an S/I/R call were not missing data
but uninterpreted measurements: re-interpretation added 1,921 labelled genomes for meropenem,
959 for ciprofloxacin, 598 for colistin and 560 for gentamicin. It also resolved most apparent
contradictions — meropenem label conflicts fell from 60 genomes to 8 — because the majority were
not conflicting laboratory results but the same genome interpreted under different standards by
different studies. Where both a deposited call and a re-interpreted call existed, they agreed in
92.9–97.2% of records.

Three label policies were fixed in advance and are recorded as architecture decision records in
the repository:

- Intermediate and conflicting labels are excluded from the primary analysis.
- Colistin labels are restricted to broth microdilution, the only method EUCAST accepts for
  colistin MICs. Disk diffusion, agar dilution and method-not-reported colistin records are
  reported separately.
- Study-level label sets judged implausible on review are flagged and excluded, with the
  evidence recorded. One such flag was applied (Section 3.5).

Tigecycline has no EUCAST breakpoint for *K. pneumoniae* (breakpoints cover *E. coli* and
*C. koseri* only), so deposited calls were used and no borderline analysis was possible.

### 2.3 Assembly quality control

Genome metadata were retrieved for all labelled genomes. Assemblies were required to be
*K. pneumoniae* by both the species field and the genome name, 5.0–6.5 Mb, ≤500 contigs, N50
≥20 kb, CheckM completeness ≥95% and contamination ≤5%. Of 6,817 labelled genomes, 6,652 passed
(2.4% removed, predominantly for contig count). Requiring agreement between the species field
and the genome name removed four genomes whose BV-BRC records were internally inconsistent; one
was independently confirmed as *Citrobacter freundii* by Kleborate's species module.

### 2.4 Genotyping

Assemblies were downloaded from BV-BRC by genome identifier — 63% of QC-passing genomes (4,190/6,652) have no
NCBI assembly accession, so NCBI-based retrieval would have lost most of the cohort — and
validated for FASTA format and length consistency with deposited metadata.

AMRFinderPlus 4.2.7 (database 2026-08-07.1) was run in nucleotide mode with
`--organism Klebsiella_pneumoniae --plus`, and Kleborate 3.2.4 with Kaptive 3.2.2 using the
`kpsc` preset. All 6,652 genomes completed both tools (8.8 h wall time, ~140 CPU-hours, zero
failures). Tool versions, database version and container-equivalent environment files are pinned;
a cross-machine verification script reproduces identical calls on macOS arm64 and Ubuntu x86_64.

### 2.5 Rule-based prediction

Per-drug rules were defined from tool output before any phenotype comparison:

| Drug | Resistance predicted if |
|---|---|
| Meropenem | acquired carbapenemase (Subclass CARBAPENEM); porin mutations alone excluded |
| Ciprofloxacin | QRDR mutation (*gyrA*, *parC*, *parE*) or acquired quinolone gene |
| Gentamicin | gentamicin-modifying enzyme or 16S rRNA methyltransferase |
| Colistin | *mcr* genes or *mgrB*/*pmrAB*/*phoPQ* variants |
| Tigecycline | acquired *tet(X)* or *tmexCD-toprJ* |

Two rules required correction after an initial run, both identified as tool-annotation artefacts
rather than biology. AMRFinderPlus assigns the intrinsic chromosomal *oqxAB* efflux pump to
tigecycline as well as quinolones; because essentially all *K. pneumoniae* carry *oqxAB*, this
produced 88 of 91 tigecycline false positives in a 524-genome pilot. *oqxAB* is now excluded for
all drugs. Separately, efflux-regulator variants (*oqxR*, *ramR*, *acrR*) are largely lineage
polymorphisms and produced 21 of 26 pilot ciprofloxacin false positives; these are reported as a
separate rule tier ("+regulators") rather than folded into the core rule.

### 2.6 Targeted *mgrB* integrity check

The MgrB reference protein (WP_002911375.1, 47 aa) was extracted from the AMRFinderPlus database
and searched against each assembly with tblastn (e-value ≤1e-5, hits ≥10 aa). The best hit was
classified as:

- **intact** — full-length, no residue differences;
- **missense** — full-length with amino-acid substitutions;
- **premature_stop** — stop codon within the aligned protein;
- **frameshift_or_indel** — two alignment fragments on one contig, <300 bp apart;
- **insertion_assembled** — two fragments on one contig, >300 bp apart (an assembled insertion);
- **split_contigs** — fragments on different contigs;
- **contig_break** — alignment terminating at a contig edge;
- **partial** — a fragment only, not at a contig edge;
- **absent** — no hit.

All categories except intact and missense are counted as disrupted. The split_contigs and
contig_break categories are precisely the assembly signature of IS insertion described in
Section 1.

Residue differences are enumerated by direct comparison of aligned query and subject sequences,
ignoring query position 1: tblastn translates alternative start codons (GTG, TTG) as V or L
rather than M. An initial implementation that used percentage identity classified 6,048 genomes
as missense at exactly 97.9% identity — a single artefactual mismatch at the start codon. The
corrected method returns 6,056 intact genomes, which is consistent with *mgrB* being conserved in
most isolates. The whole search takes 2 minutes for 6,652 assemblies on 16 threads.

### 2.7 Benchmark

Features were 226 binary indicators per genome: AMRFinderPlus gene families (allele suffixes
collapsed, e.g. *bla*KPC-2 → *bla*KPC, but OXA-48-like carbapenemases kept distinct from other
OXA enzymes), 100 point-mutation calls retained in full, 33 AMRFinderPlus functional subclasses,
and two *mgrB* features. Features present in fewer than five genomes were dropped.

Three predictors were compared: the per-drug rule (for colistin, including disrupted or missense
*mgrB*); L2 logistic regression (C = 1, class-balanced); and gradient boosting
(HistGradientBoostingClassifier, 200 iterations, 15 leaves, class-balanced). Hyperparameters were
fixed, not tuned, because the object of study is validation design rather than peak model
performance.

Four split schemes were used: random 5-fold; study-aware 5-fold grouped by BioProject;
lineage-aware 5-fold grouped by sequence type; and temporal (train on earlier genomes, test on
the most recent ~20% of dated genomes). Genomes with missing group metadata were assigned
singleton groups so that absent metadata can never merge unrelated isolates. Metrics are computed
on pooled held-out predictions, with 95% Wilson intervals, and are reported both on all test
labels and with borderline-MIC labels removed. A split is scored only if its test set contains at
least 10 resistant and 10 susceptible genomes.

### 2.8 Reproducibility

All analyses are deterministic for a fixed seed (20261024); repeat runs reproduced identical
results. Code, pinned environments, provenance manifests, decision records and all generated
reports are available at https://github.com/omarfaruqe/klebsiella-amr. The benchmark completes in
52 seconds; the *mgrB* check in 2 minutes; the full genotyping pipeline in approximately 9 hours
on 32 cores.

---

## 3. Results

### 3.1 Cohort

The final cohort comprised 6,652 QC-passing *K. pneumoniae* genomes with harmonised EUCAST v16.1
laboratory labels: meropenem 5,541 (1,854 R / 3,687 S), ciprofloxacin 4,640 (3,441 / 1,199),
gentamicin 4,760 (2,113 / 2,647), colistin 1,528 broth-microdilution (367 / 1,161), tigecycline
948 (261 / 687).

The cohort is structurally non-random, as expected of public data. Two BioProjects contribute 44%
of genomes, and the ten commonest sequence types account for 52%, led by ST307 (11.5%), ST258
(10.1%) and ST15 (7.3%). Collection year is missing for 34% of genomes (2,266/6,652). These properties motivate
the grouped validation in Section 3.4 and constrain the temporal analysis.

### 3.2 Known mechanisms explain four drugs well and colistin poorly

Rule-based prediction against harmonised phenotype across the full cohort:

| Drug | AMRFinderPlus sensitivity | specificity | Kleborate sensitivity | specificity |
|---|---|---|---|---|
| Ciprofloxacin | 98.8% (98.4–99.1) | 88.2% | 98.7% | 88.4% |
| Meropenem | 96.7% (95.7–97.4) | 91.2% | 96.5% | 91.3% |
| Gentamicin | 85.1% (83.6–86.6) | 94.9% | 81.4% | 95.8% |
| **Colistin** | **28.3% (24.0–33.2)** | 97.8% | **39.5%** | 96.3% |
| **Tigecycline** | **2.7% (1.3–5.4)** | 100.0% | **0.0%** | 100.0% |

The two tools agree with each other on 97.8–99.9% of genome–drug pairs for meropenem,
ciprofloxacin, gentamicin and tigecycline, but on only 87.8% for colistin — the one drug where
they differ substantially in what they look for. Kleborate's higher colistin sensitivity derives
from its truncation calls, which AMRFinderPlus's point-mutation screen does not make.

For colistin, 263 of 367 resistant genomes carried no AMRFinderPlus colistin determinant at all.

### 3.3 A targeted *mgrB* check more than doubles colistin sensitivity

Across all 6,652 genomes, *mgrB* was intact in 6,056, carried missense variants in 129, and was
disrupted in 467. The disrupted categories were dominated by the assembly signatures predicted in
Section 1: split across contigs (162), partial (153), premature stop (73), absent (34),
frameshift or indel (27), contig break (13) and assembled internal insertion (5).

*mgrB* status was strongly associated with colistin phenotype among the 1,528 genomes with
broth-microdilution labels (overall 24% resistant):

| *mgrB* status | Resistant | Susceptible |
|---|---|---|
| Disrupted (all categories) | 153 | 41 |
| Missense | 70 | 13 |
| Intact | 144 | 1,107 |

Adding *mgrB* status to the rules produced the principal result (Figure 1):

| Rule | Sensitivity (95% CI) | Specificity (95% CI) |
|---|---|---|
| AMRFinderPlus | 28.3% (24.0–33.2) | 97.8% (96.7–98.5) |
| AMRFinderPlus + *mgrB* disrupted | 56.7% (51.6–61.7) | 95.2% (93.8–96.3) |
| **AMRFinderPlus + *mgrB* disrupted or missense** | **67.3% (62.3–71.9)** | **94.3% (92.8–95.5)** |
| Kleborate | 39.5% (34.6–44.6) | 96.3% (95.0–97.2) |
| Kleborate + *mgrB* disrupted or missense | 61.6% (56.5–66.4) | 94.5% (93.0–95.7) |
| *mgrB* disrupted or missense alone | 60.8% (55.7–65.6) | 95.3% (94.0–96.4) |

![Figure 1](../../reports/figures/colistin_mgrb_light.png)

**Figure 1.** Colistin sensitivity and specificity for each rule, with and without the targeted
*mgrB* check. Dots are point estimates; bars are 95% Wilson intervals. Broth-microdilution labels
(n = 1,528; 367 resistant). Note the specificity axis begins at 85%.

Sensitivity rose from 28.3% to 67.3%, explaining 143 additional resistant isolates, for 40
additional false positives — 3.5 percentage points of specificity. The *mgrB* check used alone,
with no gene screen at all, reached 60.8% sensitivity, exceeding both standard tools.

The missense tier contributed materially and is reported separately rather than merged, because
not all substitutions are functional. Variants enriched among resistant isolates included M27K
(25 R / 3 S), T21N (9 / 0), N42Y (7 / 0), K43I (7 / 0), C28W (4 / 1) and C39G (3 / 0). C28 and
C39 are the two cysteine residues of this small lipoprotein, consistent with substitutions there
being structurally disruptive.

A total of 144 resistant isolates retained an intact *mgrB* coding sequence. These are candidates
for *pmrB*, *phoQ* or *crrB* variants, or for IS insertion in the *mgrB* promoter region, which a
coding-sequence search does not interrogate.

### 3.4 Under lineage-aware validation the *mgrB* rule beats both models

For colistin, the three predictors behave very differently as validation becomes stricter:

| Split | Rule + *mgrB* | Logistic regression | Gradient boosting |
|---|---|---|---|
| Random | 0.808 | 0.844 | 0.842 |
| Study-aware | 0.808 | 0.822 | 0.810 |
| Lineage-aware | **0.808** | 0.731 | 0.709 |

(balanced accuracy on held-out genomes)

![Figure 2](../../reports/figures/benchmark_splits_light.png)

**Figure 2.** Balanced accuracy of the per-drug rule, logistic regression and gradient boosting
under random, study-aware and lineage-aware splits, for each drug. Lines that fall from left to
right indicate performance that a random split attributed to the model but that does not survive
grouping. Rules are deterministic and therefore flat by construction.

Under random splitting, both models appear to beat the rule. Under lineage-aware splitting both
fall well below it. The collapse is concentrated in sensitivity: logistic regression falls from
77.7% to 58.3% and gradient boosting from 81.2% to 50.4%, while the rule is 67.3% by construction
under every scheme, since it depends only on the genome in front of it.

The pattern generalises across drugs. Random-split optimism — balanced accuracy under random
splitting minus the lower of the two grouped splits — was 0.007–0.042 for meropenem,
ciprofloxacin and gentamicin, but 0.113–0.133 for colistin and 0.136–0.149 for tigecycline. The
apparent advantage of machine learning was largest precisely where it was least real. For
tigecycline, models appeared to improve on the rule from 0.513 to 0.769 under random splitting,
but reached only 0.620–0.626 under study-aware splitting.

Where rules already perform well, models added 1–2 points under random splitting and little under
grouped splitting; for gentamicin the rule (0.906) exceeded both models under the study-aware
split (0.882–0.886).

Temporal validation was interpretable only for meropenem, ciprofloxacin and gentamicin
(546–569 test genomes). For colistin the post-2014 test set contained 34 genomes with a single
resistant isolate, and for tigecycline 28 genomes with none; these are reported as not evaluable
rather than scored.

### 3.5 Much apparent genotype error is label artefact

**Meropenem.** Of 2,118 carbapenemase carriers, 326 tested susceptible. These were strongly
enzyme-dependent: 52.9% of OXA-48-like carriers (182/344) and 75.0% of VIM carriers (33/44)
tested susceptible, against 5.9% of KPC (75/1,278) and 10.0% of NDM (28/280) carriers. Of the 240
such isolates with a recorded MIC, 191 (80%) sat at 1–2 mg/L and 224 (93%) exceeded the EUCAST
carbapenemase screening cut-off of 0.125 mg/L. The genotype call is therefore correct in most
cases; the clinical S category simply does not reflect carbapenemase carriage for these enzyme
families. Conversely, of 62 resistant non-carriers, 51 carried an ESBL or AmpC enzyme and 29 of
those also had an OmpK35/36 porin defect — but 624 susceptible non-carriers had the same
combination, so this route is a risk factor rather than a predictor.

**Gentamicin.** Of 299 unexplained resistant isolates with an MIC, 213 (71%) were at exactly
4 mg/L, one doubling dilution above the EUCAST breakpoint of 2, whereas resistant isolates with a
detected mechanism were mostly at ≥8 mg/L. Among the CLSI-tested records in this borderline
group, 160 of 162 had been reported susceptible by the originating laboratory. The rule agrees
with the original laboratory; EUCAST re-interpretation created the discordance. Excluding
borderline labels raises gentamicin rule balanced accuracy from 0.906 to 0.959, above every
model under every split.

**A flagged study.** Of 135 gentamicin false positives, 65 carried *rmtB*, a 16S rRNA
methyltransferase conferring high-level aminoglycoside resistance; 63 of these came from a single
BioProject and were reported at MIC 1 mg/L. Such values are not biologically plausible for *rmtB*
carriage. That study's gentamicin labels are flagged and excluded from the primary benchmark,
with the evidence recorded in the repository.

---

## 4. Discussion

The central finding is narrow and practical: for colistin in *K. pneumoniae*, a two-minute
tblastn check of a single 47-codon locus recovers more resistance than either widely used
genotyping tool, and more than either machine-learning model once those models are tested on
unseen lineages.

This is not a failure of the existing tools so much as a mismatch between their design and the
mechanism. AMRFinderPlus and Kleborate are built to detect genes that are present and mutations
that are callable. IS disruption of *mgrB* produces neither: it produces absence, and in
short-read assemblies, absence that manifests as a contig boundary. The categories that carried
most signal here — split contigs, partial alignments, contig breaks — are not sequence features
at all but assembly artefacts that happen to be diagnostic. Any screen that asks "which genes are
here?" will miss them; a screen that asks "is this specific locus whole?" will not.

The contrast with the machine-learning models is instructive. Given the same assemblies, models
had every opportunity to learn the same signal, and under random splitting they appeared to. But
what they largely learned was lineage. Colistin resistance in this cohort is concentrated in
particular clones, so a model can score well by recognising the clone — until the clone is held
out, at which point sensitivity falls by roughly 20–30 points. The rule, which inspects only the
mechanism, is unaffected. Yu and colleagues demonstrated the general phenomenon across five
species [1]; the present work shows what follows from it in a single well-characterised case, and
answers their call for rule-based baselines: the baseline wins.

The label findings carry a methodological message of similar weight. Roughly a third of our
usable labels did not exist as S/I/R calls at all — they were uninterpreted MICs that
harmonisation converted into data. Conversely, a substantial fraction of apparent genotype error
was created by harmonisation: the gentamicin borderline effect moves rule balanced accuracy by
five points, and depends entirely on whether one applies EUCAST or the CLSI revision under which
those isolates were originally reported. A benchmark that silently selects a standard, or mixes
standards, is reporting a number whose provenance the reader cannot assess. We therefore report
every result twice, with and without borderline labels, and treat the difference as label
uncertainty rather than model error.

The meropenem result is a reminder that "false positive" is a claim about the label as much as
the prediction. Three quarters of VIM carriers and half of OXA-48-like carriers in this cohort
test clinically susceptible to meropenem, and 93% of those with an MIC nonetheless exceed the
screening cut-off that EUCAST provides precisely because carbapenemase carriage and clinical
susceptibility routinely co-occur. Scoring a genomic tool against the clinical category penalises
it for being right about the enzyme.

### 4.1 Limitations

Public surveillance data over-represent outbreaks, carbapenem-resistant isolates and high-risk
clones; two BioProjects supply 44% of our genomes and the resistance prevalences here should not
be read as population estimates. Labels originate from many laboratories using several methods,
and residual label error certainly remains beyond the one study we flagged. Colistin and
tigecycline have modest numbers of resistant isolates (367 and 261), which bounds the precision
of those estimates. Collection year is missing for 34% of genomes and the dated genomes
concentrate in 2010–2014, so temporal validation was possible for only three drugs. Our features
are presence/absence calls from assemblies: gene copy number, expression, plasmid context and
promoter-region insertions are not represented, and the 144 colistin-resistant isolates with
intact *mgrB* coding sequence likely include promoter insertions we cannot see. Hyperparameters
were fixed rather than tuned, so the models here are reasonable defaults, not optimised systems;
a tuned model might close part of the lineage-aware gap, though tuning cannot manufacture a
mechanism signal the features do not contain. Finally, *mgrB* status was validated against
phenotype but not against long-read sequencing, which would be required to confirm that the
contig-break category corresponds to IS insertion in each individual genome.

### 4.2 Recommendations

For evaluations of genomic AMR prediction we suggest four practices, each inexpensive:

1. **Report a per-drug rule baseline.** It is the only way to tell whether a model has learned a
   mechanism or a lineage.
2. **State and harmonise the breakpoint standard, and report sensitivity to borderline labels.**
   The standard alone moved one of our drugs by five points.
3. **Use grouped splits — by study and by lineage — and report the gap from random splitting.**
   The gap is itself a useful statistic: it was under one point where rules work and 11–15 points
   where they do not.
4. **When a mechanism is poorly detected, examine the failure mode before adding model capacity.**
   Here the failure mode was an assembly artefact, and the fix was 100 lines of code.

---

## 5. Data and code availability

All code, pinned environments, provenance manifests (with retrieval queries, timestamps and
SHA-256 checksums), architecture decision records, generated reports and figures are available at
https://github.com/omarfaruqe/klebsiella-amr. Genome assemblies and AST records are from BV-BRC.
Derived tables are regenerated deterministically by the pipeline documented in the repository
README. A model card for the benchmark is provided at `docs/model-card.md`.

## 6. Competing interests

None declared.

## 7. References

All references below were verified against the publisher record during preparation
(2026-10-24).

1. Yu Y, Wheeler NE, Barquist L. Biased sampling driven by bacterial population structure
   confounds machine learning prediction of antimicrobial resistance. *PLOS Biology*. 2025;
   published 16 December 2025. doi:10.1371/journal.pbio.3003539
2. The European Committee on Antimicrobial Susceptibility Testing. Breakpoint tables for
   interpretation of MICs and zone diameters. Version 16.1, valid from 2026-06-24.
   https://www.eucast.org
3. Lam MMC, Wick RR, Watts SC, Cerdeira LT, Wyres KL, Holt KE. A genomic surveillance framework
   and genotyping tool for *Klebsiella pneumoniae* and its related species complex. *Nature
   Communications*. 2021;12:4188. doi:10.1038/s41467-021-24448-3
4. Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, Hoffmann M,
   Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. AMRFinderPlus and the Reference Gene
   Catalog facilitate examination of the genomic links among antimicrobial resistance, stress
   response, and virulence. *Scientific Reports*. 2021;11(1).
   doi:10.1038/s41598-021-91456-0
5. Olson RD, Assaf R, Brettin T, et al. Introducing the Bacterial and Viral Bioinformatics
   Resource Center (BV-BRC): a resource combining PATRIC, IRD and ViPR. *Nucleic Acids Research*.
   2023;51(D1):D678–D689. PMID:36350631
6. Pedregosa F, Varoquaux G, Gramfort A, et al. Scikit-learn: machine learning in Python.
   *Journal of Machine Learning Research*. 2011;12:2825–2830.
