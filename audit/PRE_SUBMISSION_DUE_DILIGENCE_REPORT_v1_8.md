# RCC TKI-ICI WES/HLA v1.8: pre-submission due-diligence report

**Audit date:** 2026-07-17  
**Package status:** evidence-frozen pre-submission draft; not yet ready for formal journal upload.

## What v1.8 changes materially

1. **SBS/APOBEC was removed from the scientific claim.** The old signature matrix was not PASS-filtered in a way that supports the historical subgroup. Under PASS biallelic SNVs, COSMIC SBS v3.2, and 100 bootstrap fits, EIS01, LEN08, and LEN15 have 0/100 detection for SBS2 and SBS2+SBS13. Any APOBEC-based clinical conclusion is unsupported.
2. **The local historical `TMB.tab` is no longer called TMB.** It is a retained-record count normalized by a fixed 50 Mb denominator and has weak rank correlation with the official Lee MAF protein-altering count (Spearman rho=0.166, P=0.438). v1.8 uses the official count and calls it `protein-altering mutation count`.
3. **Discovery PFS was reconciled.** It contains 16 events and 8 censored patients. The conflicting historical event field was excluded.
4. **The discovery state was reclassified.** All three HLA-A-homozygous cases are above the mutation-count median and have early PFS events, but the low-count/HLA-A-homozygous cell is empty. Mutation and HLA-A effects cannot be separated; no discovery interaction is estimable.
5. **The paired-peptide result was de-circularized.** The old 29/29 result conditioned on `partner_only` selection. The complete archived retained 63-event universe has no global partner advantage (29/63; P=0.615); C-terminal residue class instead determines predicted allocation.
6. **External evidence was graded by directness.** Miao is mutation-only; PXD017149 is ligandome context only; MSK is a post hoc conditional OS association, not direct HLA-A/TKI-ICI/PFS replication. Braun, JAVELIN, Chowell, and Nature Cancer define non-replication boundaries.

## Claim-to-evidence ceiling

| Claim | Evidence status | Permitted wording | Prohibited wording |
|---|---|---|---|
| Lee n=3 observation | Record-level, post hoc | `post hoc discovery-cohort association` | `validated biomarker` |
| Miao n=35 | Independent mutation-only association | `mutation-only treated-cohort association` | `HLA/HED validation` |
| MSK n=138 | Conditional FDR-retained association; interaction FDR null | `post hoc conditional high-TMB subgroup association` | `direct replication` or `interaction proof` |
| 63-event audit | Computational event-level association | `terminal-residue-dependent predicted compatibility` | `presentation rescue`, natural presentation, or T-cell recognition |
| PXD017149 | Unlinked real-ligand background | `ligandome context` | clinical, patient-specific, or peptide validation |

The full claim matrix is `CLAIM_EVIDENCE_MATRIX.csv`.

## Deliverables verified

- Main figures: `figures/Figure1.png` through `Figure5.png`.
- Supplementary figures: `supplementary_figures/Supplementary_Figure_S1.png` through `Supplementary_Figure_S6.png`.
- Figure visual review: all five main figures and Supplementary Figure S5 were inspected at full resolution. Labels, legends, and panel alignment are readable; the visual message matches the permitted claims above.
- Source data workbook: `RCC_TKI_ICI_WES_HLA_v1_8_source_data.xlsx`, 24 sheets including README, all figure data, variable dictionary, and prior-knowledge table.
- Supplementary tables: S1-S8 in `supplementary_tables/`.
- English and Chinese Word files: generated with five embedded main figures each; supplementary Word file has six embedded supplementary figures.
- Reference audit: 23/23 DOI records resolve through Crossref; the ledger is `Reference_metadata_audit_v1_8.csv`.
- Reproducibility package: `analysis_code/` includes the analysis scripts and run-order manifest.

## Remaining hard blockers

1. **Permanent data/code deposit:** create a public repository and DOI for code plus permitted aggregate source data, then replace the Code availability placeholder. Do not deposit modified patient-level MSK records.
2. **Author-supplied declarations:** authors, affiliations, corresponding author, funding, conflicts, CRediT roles, ethics determination, and journal-specific AI disclosure cannot be inferred and are blank by design.
3. **Target-journal conversion:** journal-specific reference style, title page, word count, abstract format, figure format, supplementary naming, and data/AI wording remain target-dependent.
4. **Word visual render:** DOCX structural checks passed, but full page-image QA could not run because this workstation lacks both LibreOffice (`soffice`) and a usable Word automation server. Perform one final Word/PDF visual review on a machine with either application before upload.

## Scientific residual risks (not fixable by editing)

- The primary Lee observation is n=3 and post hoc.
- No low-count/HLA-A-homozygous discovery participant exists.
- No external cohort simultaneously provides matched TKI-ICI, PFS, HLA/HED, harmonized mutation data, and peptide-level readouts.
- Peptide predictions are not immunopeptidomics or T-cell assays.
- Cross-cohort heterogeneity prevents a pooled causal or predictive estimate.

## Submission recommendation

The v1.8 manuscript should be positioned as a **transparent multi-cohort reanalysis and falsification-aware boundary analysis**, not a biomarker study. Its realistic acceptance range remains IF 3-5. A selective IF 5-6 attempt is defensible only after the hard blockers are completed and only with the claim ceiling preserved. No version of the current data supports an IF 8+ biomarker claim.
