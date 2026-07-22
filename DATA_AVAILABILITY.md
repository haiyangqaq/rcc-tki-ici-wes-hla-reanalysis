# Data and code availability

## Repository content

This repository contains code, figures, aggregate source data, and analysis outputs generated for the v1.8 reanalysis. It does not contain raw sequencing data or redistributed patient-level third-party data.

## Reused source data

| Source | Identifier / route | Role in this work | Redistribution in this repository |
|---|---|---|---|
| Lee lenvatinib-pembrolizumab RCC study | BioProject PRJNA610643; Lee et al., Molecular Cancer Research 2021, DOI 10.1158/1541-7786.MCR-21-0053 | Discovery WES/HLA/PFS reconstruction | No patient-level redistribution |
| Miao/DFCI ccRCC | cBioPortal `ccrcc_dfci_2019`; Miao et al., Science 2018, DOI 10.1126/science.aan5951 | Mutation-only treated-cohort analysis | No patient-level redistribution |
| MSK RCC ICI and MSK50K HLA | cBioPortal `tmb_mskcc_2018` and `msk_impact_50k_2026` | Conditional/interaction sensitivity analysis | No patient-level redistribution; source licenses apply |
| ccRCC HLA ligandome | ProteomeXchange PXD017149 | Ligandomics context | No peptide-level redistribution |

## Code

Code is provided in `scripts/`. Before running scripts requiring source inputs, set the documented environment variables and obtain data from the original source under its terms. No credential, token, raw clinical record, or modified MSK patient-level table is included.

## Citation

Use the DOI assigned to this release once published. Cite the original studies and accessions above for all reused data.
