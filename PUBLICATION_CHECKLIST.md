# GitHub and Zenodo release checklist

## Release identity

- Repository name: `rcc-tki-ici-wes-hla-reanalysis`
- Visibility: **Public**. This is required for the standard Zenodo-GitHub release archive workflow.
- Release tag: `v1.0.2`
- Release title: `v1.0.2 - Verified Zenodo archival release`

## Publish to GitHub

Create an empty public GitHub repository with the repository name above. Do not initialize it with a README, licence, or `.gitignore`; all three are already supplied here.

From this directory, run:

```powershell
git remote add origin https://github.com/<GITHUB_ACCOUNT>/rcc-tki-ici-wes-hla-reanalysis.git
git push -u origin main
git push origin v1.0.2
```

Then create a GitHub Release from tag `v1.0.2`. Use the release notes below:

> This frozen release contains code, aggregate source data, figures, and audit materials for an exploratory multi-cohort RCC immunogenomics reanalysis. Patient-level third-party data, modified MSK-derived records, raw sequencing files, and local source caches are intentionally excluded. See `DATA_AVAILABILITY.md` and `audit/CLAIM_EVIDENCE_MATRIX.csv` before reusing results.

## Archive with Zenodo

1. Sign in to Zenodo using the same GitHub account.
2. In Zenodo account settings, enable the GitHub integration and activate this repository.
3. Create the GitHub Release `v1.0.2` after activation. Zenodo will automatically archive it and mint a version-specific DOI. The preceding releases were published before the final checksum verification workflow was completed.
4. Review Zenodo metadata imported from `.zenodo.json` before publishing. Confirm the creator name, title, version, description, and licence.
5. Record both the concept DOI and version DOI. Cite the version DOI in the manuscript and include the concept DOI where an always-latest DOI is useful.

## Required post-publication edits

Replace the DOI placeholder in the manuscript Code Availability section with the version DOI. Confirm that the GitHub release and Zenodo file lists contain no patient-level record, raw VCF/BAM/CRAM/FASTQ file, private path, credential, or third-party file outside the documented aggregate outputs.
