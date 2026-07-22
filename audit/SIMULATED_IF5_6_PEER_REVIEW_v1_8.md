# Simulated IF 5-6 peer review: v1.8

**Manuscript assessed:** `RCC_TKI_ICI_WES_HLA_manuscript_english_v1_8.md`  
**Review framework:** editor, biostatistics, RCC immunotherapy, immunogenomics/peptidomics, and adversarial reviewer perspectives.  
**Decision:** Major revision before submission. The evidence package is now materially more honest and reproducible than v1.7, but the manuscript is not yet submission-ready because it lacks a permanent code/data deposit, author declarations, and a target-journal-specific reference style.

## Editorial assessment

The paper is not a clinical biomarker paper. Its defensible contribution is an unusually transparent reanalysis that shows why a striking discovery-cohort subgroup cannot identify a mutation-by-HLA interaction, maps discordant RCC cohorts, and removes circular peptide evidence. That is potentially publishable as an exploratory computational/reanalysis article only if the journal accepts negative or boundary-defining analyses.

The title, abstract, and conclusion are now aligned with this framing. However, an IF 5-6 oncology/immunotherapy journal may still consider the clinical advance limited because the most direct discovery association has three exposed patients, the Miao cohort lacks HLA variables, the MSK result is post hoc and uses OS/ICI rather than PFS/TKI-ICI, and no patient-specific immunopeptidomics or T-cell experiment is available.

**Editorial recommendation:** Do not submit as a biomarker-discovery paper. Submit only after the mandatory operational fixes below, to a journal that accepts rigorous secondary analyses and explicitly allows negative results. A realistic submission sequence is Cancer Immunology, Immunotherapy or OncoImmunology only after a final scope check; the more conservative IF 3-5 route remains more probable.

## Reviewer 1: biostatistics and reproducibility

### Major comments

1. **Post hoc subgroup inference remains the central limitation.** The manuscript correctly explains that the n=3 discovery state is identical to HLA-A homozygosity and has an empty low-burden/homozygous cell. The log-rank and fixed-size permutation results must remain descriptive diagnostics. They cannot be presented as evidence for a biomarker or an interaction. This is addressed in the Results and Discussion, but the manuscript must preserve the same restraint in a cover letter and graphical abstract.

2. **The MSK conditional result is a multiple-testing-sensitive secondary observation.** The 20-test BH family is stated and the conditional result survives it. However, the combination of TMB stratum, histology scope, and HLA construct was searched. The result should be labeled `post hoc conditional association` in the figure panel itself and in the abstract. The manuscript must not imply that leave-one-out robustness corrects selection.

3. **Miao covariate handling requires a source-data dictionary.** First-line status and nivolumab dose need coded definitions, missing-data handling, and a stated sample count for the adjusted model in either Methods or Supplementary Table S3. The source-data workbook includes the patient-level fields, but a formal variable dictionary must be added before submission.

4. **Peptide event-level dependence is unresolved.** Sixty-three events are nested within seven patients. The Fisher and binomial tests describe event-level allocation and should not be read as seven independent patients. Add one sentence in Methods and Discussion that p values are event-level descriptive tests and do not provide patient-level inference.

### Minor comments

- State the exact IQR calculation convention for mutation count (empirical 75th minus 25th percentile).
- Report whether there were missing covariate values in the Miao and MSK adjusted models.
- Use `multiplicity-adjusted FDR` rather than an unqualified `FDR` where a panel contains both nominal and adjusted values.

## Reviewer 2: RCC immunotherapy and clinical interpretation

### Major comments

1. **Treatment context is a strength only if the comparisons are disciplined.** The manuscript should not call Braun, Miao, JAVELIN, Chowell, Nature Cancer, and MSK a validation series. It appropriately calls them contextual cohorts, but Figure 4 needs a conspicuous endpoint/treatment table or caption reminder. The final figure does this partially; retain it in any journal-resized version.

2. **Prognostic versus predictive is correctly separated but should appear earlier.** Add a sentence in the Abstract or Introduction: all analyses concern outcomes among treated patients and cannot establish treatment-predictive utility without non-ICI comparators.

3. **Clinical actionability is absent.** The conclusion should explicitly say the results should not influence treatment selection. The current final Discussion sentence does this and should be retained.

### Minor comments

- Use one term consistently: `protein-altering mutation count` for Lee/Miao and `TMB` only where a study provides mutations per megabase.
- Avoid `high mutation burden` in the Lee section; say `at or above the cohort median mutation count`.

## Reviewer 3: HLA immunogenomics and peptide biology

### Major comments

1. **The correction of the 29-event selection is the strongest scientific contribution.** It must be foregrounded in the abstract and Figure 5 legend as `complete archived retained universe`, not as a complete neoantigen universe from the raw data. The retained events are already downstream of an earlier pipeline.

2. **Predicted terminal grammar is not evidence of natural presentation.** PXD017149 should be described only as ligandome context. The text meets this standard. Do not use `validates`, `confirms patient peptides`, `immunogenic`, or `T-cell recognition` for this result.

3. **HLA genotype provenance is appropriately repaired.** The published HLA calls must be used consistently across all primary figures. The package documents the LEN07 A*34:01 correction. A final automated scan should ensure that the superseded A*34:05 call does not appear in the manuscript, legends, tables, or code comments.

## Reviewer 4: novelty and literature positioning

### Major comments

1. **Novelty is methodological and falsification-oriented, not discovery-oriented.** The paper adds little if framed as `HLA diversity predicts response`, because the source article and Chowell literature already make that point. The publishable increment is source reconstruction, non-identifiability, cross-cohort heterogeneity, and correction of circular peptide selection.

2. **The paper needs a concise prior-knowledge versus present-work table.** Add Supplementary Table S6 or a boxed paragraph that states what Lee, Chowell, SNiP-RCC, and the present reanalysis respectively test. This will prevent reviewers from confusing a non-replication with a contradiction.

3. **Reference formatting had substantive errors in the draft.** These have been corrected against Crossref for references 7, 8, 9, 16, and 18. The remaining journal-specific author formatting must be applied only after the target journal is selected.

## Reviewer 5: adversarial review

### Primary challenge

The manuscript could be criticized as an internally generated negative result built from heterogeneous public datasets. The defense is not to inflate the conditional MSK association. The defensible response is transparency: show every search family, preserve all null boundaries, state directness, keep the data dictionary and code public, and remove all language that lets readers infer clinical validation.

### Residual vulnerabilities that cannot be repaired with editing

- n=24 discovery cohort and n=3 HLA-A-homozygous cases.
- No low-mutation/HLA-A-homozygous discovery case; no interaction can be identified.
- No external cohort with matched TKI-ICI, PFS, HLA genotype/HED, mutation burden, and harmonized peptide data.
- No direct immunopeptidomics or T-cell functional validation.
- Heterogeneous definitions, assay modalities, endpoints, and treatments prevent pooled effect estimation.

These are evidence ceilings, not wording defects. The manuscript must retain them as limitations.

## Mandatory actions before formal submission

1. Create a permanent public repository/DOI for code and permitted aggregate source data; add the identifier to Code availability.
2. Add an Miao/MSK variable dictionary: definitions, coding, missingness, and adjusted-model sample counts.
3. Add an event-clustering caveat for Figure 5 statistics and label the universe as `complete archived retained` everywhere.
4. Add a prior-knowledge versus present-work supplementary table.
5. Complete all author-entered declarations: authorship, affiliations, corresponding author, funding, conflicts, ethics status, and journal-compliant AI disclosure.
6. Choose the target journal and convert references, title page, figure file types, word count, and declaration wording to its exact instructions.

## Actions already verified in v1.8

- SBS/APOBEC labels were removed after PASS-SNV bootstrap instability.
- The historical retained-record table is not labeled TMB.
- Official Lee HLA calls and locked 16-event/8-censored PFS provenance are used.
- The circular 29-event peptide subset was replaced by the complete archived 63-event retained universe.
- MSK patient-level data are not redistributed under the stated CC BY-NC-ND constraint.
- Every cited DOI in the 23-reference v1.8 manuscript resolves through Crossref; five incorrect reference records were corrected.

## Submission judgment

**Current scientific status:** coherent exploratory reanalysis with a rigorous falsification component.  
**Current operational status:** pre-submission only; the six mandatory actions remain.  
**Recommended claim ceiling:** treatment-context-dependent associations and a testable replication boundary, not a validated RCC biomarker.  
**Expected outcome if submitted without the mandatory actions:** high risk of editorial rejection or major revision.  
**Expected outcome after mandatory actions:** plausible external review at an IF 3-5 journal; a selective IF 5-6 submission is defensible but remains a higher-risk attempt.
