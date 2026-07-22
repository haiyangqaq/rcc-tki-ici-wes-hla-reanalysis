"""Compare published and internal HLA genotype versions for the discovery cohort."""

from pathlib import Path
import os

import pandas as pd


OUT = Path(__file__).resolve().parent
OFFICIAL = Path(os.environ["RCC_INPUT_AUDIT_HLA_VERSION_PROVENANCE_01"])
INTERNAL = Path(os.environ["RCC_INPUT_AUDIT_HLA_VERSION_PROVENANCE_02"])
ANALYSIS = Path(os.environ["RCC_INPUT_AUDIT_HLA_VERSION_PROVENANCE_03"])


def split_genotype(value: str) -> dict:
    alleles = [part.strip() for part in str(value).split(",")]
    return {
        "A": ";".join(sorted(alleles[0:2])),
        "B": ";".join(sorted(alleles[2:4])),
        "C": ";".join(sorted(alleles[4:6])),
    }


def main() -> None:
    official = pd.read_excel(OFFICIAL, sheet_name="lenvatinib pembrolizumab cohort")
    official_parts = pd.DataFrame([split_genotype(value) for value in official["HLA class-I type"]]).add_prefix("official_")
    official = pd.concat([official[["Patient", "HLA class-I type", "Mean HED"]], official_parts], axis=1)

    internal = pd.read_csv(INTERNAL)
    internal_parts = pd.DataFrame([split_genotype(value) for value in internal.HLA_Alleles]).add_prefix("internal_")
    internal = pd.concat([internal[["PatientID", "HLA_Alleles", "HED_A", "HED_B", "HED_C", "Mean_HED"]], internal_parts], axis=1).rename(columns={"PatientID": "Patient"})

    audit = official.merge(internal, on="Patient", validate="one_to_one")
    analysis = pd.read_csv(ANALYSIS)[["Patient", "Mean_HED"]].rename(columns={"Mean_HED": "analysis_Mean_HED"})
    audit = audit.merge(analysis, on="Patient", validate="one_to_one")
    for locus in "ABC":
        audit[f"{locus}_match"] = audit[f"official_{locus}"].eq(audit[f"internal_{locus}"])
    audit["all_loci_match"] = audit[["A_match", "B_match", "C_match"]].all(axis=1)
    audit["mean_HED_difference_internal_minus_official"] = audit.Mean_HED - audit["Mean HED"]
    audit["mean_HED_difference_analysis_minus_official"] = audit.analysis_Mean_HED - audit["Mean HED"]
    audit.to_csv(OUT / "discovery_hla_version_comparison.csv", index=False)

    discordant = audit.loc[~audit.all_loci_match]
    large_hed = audit.loc[audit.mean_HED_difference_analysis_minus_official.abs().gt(0.1)]
    report = f"""# Discovery HLA-version provenance audit

- Published and internal HLA genotype pairs match at all A/B/C loci for {int(audit.all_loci_match.sum())}/24 patients.
- Genotype-discordant patients: {', '.join(discordant.Patient) if len(discordant) else 'none'}.
- Patients with absolute mean-HED difference >0.1 in the historical analysis table: {', '.join(large_hed.Patient) if len(large_hed) else 'none'}.
- LEN15 has a concordant published/internal HLA genotype but the historical analysis table reports mean HED 5.724 instead of the published 3.967, showing that genotype concordance alone does not guarantee HED-version concordance.
- LEN07 published HLA-A is A*02:01/A*34:01; internal OptiType-derived analysis used A*02:01/A*34:05.
- The three HLA-A-homozygous early-progression cases (EIS01, LEN08, LEN15) have concordant published/internal HLA-A genotypes.

## Submission rule

Use the published Supplementary Data 2 genotype as the primary manuscript genotype. Any analysis inherited from the internal HED/OptiType table must be version-audited. LEN07 peptide-HLA predictions must be re-scored and labeled with published HLA-A*34:01; the completed re-prediction is documented separately.
"""
    (OUT / "HLA_VERSION_PROVENANCE_AUDIT_REPORT.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
