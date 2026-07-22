# Treatment-context dependence of mutation burden and HLA class I zygosity in renal cell carcinoma immunotherapy

Version 1.0.0 (2026-07-22)

## Purpose

This repository accompanies an exploratory multi-cohort reanalysis of mutation burden, HLA-I zygosity, and allele-level peptide prediction in renal cell carcinoma immunotherapy. It is a reproducibility and evidence-boundary package, **not** a clinically validated biomarker release.

## What is included

- Final English and Chinese Markdown manuscripts.
- Main and supplementary PNG figures.
- Aggregate figure source data, supplementary aggregate tables, claim-evidence matrix, and reference audit.
- Sanitized archival analysis scripts. Local absolute paths were replaced with environment variables listed in `scripts/input_environment_map.csv`.
- A SHA-256 checksum manifest.

## What is deliberately excluded

- Modified patient-level MSK records and any derived patient-level MSK merge table.
- Patient-level Miao/DFCI table, discovery patient table, event-level paired peptide table, and other per-patient/per-event derivative files.
- Raw WES/VCF/BAM files, source-study supplementary workbooks, and MHCflurry input caches.

These exclusions are intentional: the omitted records are third-party or human participant data, require source-specific permissions, or are not appropriate for unrestricted redistribution. The required source studies and accession routes are listed in `DATA_AVAILABILITY.md`.

## Reproducing the public release

The included aggregate CSV files reproduce the reported estimates and evidence maps. Full patient-level reanalysis requires independent retrieval of the listed public/controlled source data and setting the environment variables in `scripts/input_environment_map.csv`. The scripts are archival audit scripts; they do not grant redistribution rights for source datasets.

## Scientific claim boundary

- The Lee three-case state is a post hoc discovery-cohort observation, not an identifiable mutation-by-HLA interaction.
- Miao/DFCI is a mutation-only treated-cohort association.
- The MSK result is a post hoc conditional high-TMB subgroup association; its interaction family does not survive FDR correction.
- The 63-event audit supports terminal-residue-dependent predicted compatibility, not natural peptide presentation or T-cell recognition.
- PXD017149 is ligandomics context only.

See `audit/CLAIM_EVIDENCE_MATRIX.csv` and `audit/PRE_SUBMISSION_DUE_DILIGENCE_REPORT_v1_8.md` before reusing any result.
