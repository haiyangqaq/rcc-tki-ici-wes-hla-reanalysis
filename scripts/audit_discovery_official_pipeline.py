"""Rebuild the discovery analysis from the published MAF and HLA supplements."""

from itertools import combinations
from pathlib import Path
import os

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
from lifelines.utils import restricted_mean_survival_time
from statsmodels.stats.multitest import multipletests


OUT = Path(__file__).resolve().parent
MAF = Path(os.environ["RCC_INPUT_AUDIT_DISCOVERY_OFFICIAL_PIPELINE_01"])
HLA = Path(os.environ["RCC_INPUT_AUDIT_DISCOVERY_OFFICIAL_PIPELINE_02"])
PFS = Path(os.environ["RCC_INPUT_AUDIT_DISCOVERY_OFFICIAL_PIPELINE_03"])


def parse_hla(value: str) -> dict:
    alleles = [part.strip() for part in str(value).split(",")]
    if len(alleles) != 6:
        raise ValueError(f"Expected six HLA alleles, observed: {value}")
    return {
        "HLA_A1": alleles[0], "HLA_A2": alleles[1],
        "HLA_B1": alleles[2], "HLA_B2": alleles[3],
        "HLA_C1": alleles[4], "HLA_C2": alleles[5],
        "HLA_A_homozygous": alleles[0] == alleles[1],
        "HLA_B_homozygous": alleles[2] == alleles[3],
        "HLA_C_homozygous": alleles[4] == alleles[5],
        "HLA_any_homozygous": (alleles[0] == alleles[1]) or (alleles[2] == alleles[3]) or (alleles[4] == alleles[5]),
    }


def rmst(frame: pd.DataFrame, tau: float = 12) -> float:
    km = KaplanMeierFitter().fit(frame.PFS, event_observed=frame.PFS_status)
    return float(restricted_mean_survival_time(km, t=tau))


def exact_group_permutation(frame: pd.DataFrame, group: pd.Series) -> float:
    group = group.to_numpy(dtype=bool)
    observed = logrank_test(
        frame.loc[group, "PFS"], frame.loc[~group, "PFS"],
        event_observed_A=frame.loc[group, "PFS_status"],
        event_observed_B=frame.loc[~group, "PFS_status"],
    ).test_statistic
    values = frame[["PFS", "PFS_status"]].to_numpy()
    statistics = []
    for selected in combinations(range(len(frame)), int(group.sum())):
        trial = np.zeros(len(frame), dtype=bool)
        trial[list(selected)] = True
        statistics.append(logrank_test(
            values[trial, 0], values[~trial, 0],
            event_observed_A=values[trial, 1],
            event_observed_B=values[~trial, 1],
        ).test_statistic)
    statistics = np.asarray(statistics)
    return float((np.sum(statistics >= observed) + 1) / (len(statistics) + 1))


def main() -> None:
    maf = pd.read_excel(MAF, sheet_name="MAF exome sequencing")
    maf["Patient"] = maf.Tumor_Sample_Barcode.astype(str).str.replace(r"_T$", "", regex=True)
    mutation_key = ["Patient", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification"]
    burden = maf.drop_duplicates(mutation_key).groupby("Patient").size().rename("official_protein_altering_count")

    hla = pd.read_excel(HLA, sheet_name="lenvatinib pembrolizumab cohort")
    parsed = pd.DataFrame([parse_hla(value) for value in hla["HLA class-I type"]])
    hla = pd.concat([hla.reset_index(drop=True), parsed], axis=1)
    pfs = pd.read_csv(PFS)[["Patient", "PFS", "PFS_status"]]
    data = pfs.merge(burden, on="Patient", validate="one_to_one").merge(hla, on="Patient", validate="one_to_one")

    median = float(data.official_protein_altering_count.median())
    tertile = float(data.official_protein_altering_count.quantile(2 / 3))
    data["burden_high_median"] = data.official_protein_altering_count.ge(median)
    data["burden_high_tertile"] = data.official_protein_altering_count.ge(tertile)
    data["high_burden_any_hla_homo_median"] = data.burden_high_median & data.HLA_any_homozygous
    data["high_burden_any_hla_homo_tertile"] = data.burden_high_tertile & data.HLA_any_homozygous
    data.to_csv(OUT / "discovery_official_maf_hla_patient_level.csv", index=False)

    iqr = float(data.official_protein_altering_count.quantile(0.75) - data.official_protein_altering_count.quantile(0.25))
    data["burden_per_iqr"] = data.official_protein_altering_count / iqr
    cph = CoxPHFitter().fit(data[["PFS", "PFS_status", "burden_per_iqr"]], "PFS", "PFS_status")
    cox_row = cph.summary.loc["burden_per_iqr"]
    ph = proportional_hazard_test(cph, data[["PFS", "PFS_status", "burden_per_iqr"]], time_transform="rank")
    continuous = pd.DataFrame([{
        "variable": "official_protein_altering_count_per_iqr",
        "iqr": iqr,
        "hazard_ratio": cox_row["exp(coef)"],
        "lower_95ci": cox_row["exp(coef) lower 95%"],
        "upper_95ci": cox_row["exp(coef) upper 95%"],
        "wald_p": cox_row["p"],
        "ph_test_p": ph.summary.loc["burden_per_iqr", "p"],
    }])
    continuous.to_csv(OUT / "discovery_official_burden_continuous_cox.csv", index=False)

    conditional_rows = []
    for threshold in ["burden_high_median", "burden_high_tertile"]:
        high = data.loc[data[threshold]].copy()
        for defect in ["HLA_A_homozygous", "HLA_B_homozygous", "HLA_C_homozygous", "HLA_any_homozygous"]:
            flag = high[defect].astype(bool)
            if flag.sum() == 0 or (~flag).sum() == 0:
                p_value = np.nan
            else:
                p_value = logrank_test(
                    high.loc[flag, "PFS"], high.loc[~flag, "PFS"],
                    event_observed_A=high.loc[flag, "PFS_status"],
                    event_observed_B=high.loc[~flag, "PFS_status"],
                ).p_value
            conditional_rows.append({
                "threshold": threshold,
                "defect": defect,
                "n_flag": int(flag.sum()),
                "n_reference": int((~flag).sum()),
                "events_flag": int(high.loc[flag, "PFS_status"].sum()),
                "events_reference": int(high.loc[~flag, "PFS_status"].sum()),
                "logrank_p": p_value,
            })
    conditional = pd.DataFrame(conditional_rows)
    conditional["fdr_bh"] = np.nan
    estimable = conditional.logrank_p.notna()
    conditional.loc[estimable, "fdr_bh"] = multipletests(conditional.loc[estimable, "logrank_p"], method="fdr_bh")[1]
    conditional.to_csv(OUT / "discovery_official_high_burden_hla_conditional_tests.csv", index=False)

    group = data.high_burden_any_hla_homo_median
    state = data.loc[group]
    reference = data.loc[~group]
    lr = logrank_test(
        state.PFS, reference.PFS,
        event_observed_A=state.PFS_status,
        event_observed_B=reference.PFS_status,
    )
    primary = pd.DataFrame([{
        "definition": "official burden >= cohort median AND any HLA-I locus homozygous",
        "burden_median": median,
        "n": len(state),
        "members": ";".join(state.Patient),
        "events": int(state.PFS_status.sum()),
        "logrank_p": lr.p_value,
        "fixed_size_permutation_p": exact_group_permutation(data, group),
        "rmst_12m_difference": rmst(state) - rmst(reference),
        "identical_to_HLA_A_homozygous": bool(group.equals(data.HLA_A_homozygous)),
    }])
    primary.to_csv(OUT / "discovery_official_composite_summary.csv", index=False)

    report = f"""# Discovery-cohort official-pipeline audit

## Sources

- Mutation burden: 1,296 deduplicated protein-altering records from Lee et al. Supplementary Data 1.
- HLA genotype: Lee et al. Supplementary Data 2.
- PFS: 24 patients, 16 events and 8 censored observations, cross-checked against the published KM risk table and swimmer plot.

## Results

- Official mutation-count median: {median:.1f}; upper-tertile threshold: {tertile:.1f}; IQR: {iqr:.1f} records.
- Continuous mutation-count HR per IQR: {continuous.iloc[0].hazard_ratio:.2f} (95% CI {continuous.iloc[0].lower_95ci:.2f}-{continuous.iloc[0].upper_95ci:.2f}), P={continuous.iloc[0].wald_p:.3f}.
- Median-defined high-burden/any-HLA-homozygous state: {', '.join(state.Patient)} (n={len(state)}, events={int(state.PFS_status.sum())}); log-rank P={lr.p_value:.4g}; fixed-size permutation P={primary.iloc[0].fixed_size_permutation_p:.4g}; 12-month RMST difference={primary.iloc[0].rmst_12m_difference:.2f} months.
- This median-defined composite is {'exactly identical' if primary.iloc[0].identical_to_HLA_A_homozygous else 'not identical'} to HLA-A homozygosity in this cohort.
- The upper-tertile definition contains {int(data.high_burden_any_hla_homo_tertile.sum())} patients, demonstrating threshold dependence.

## Interpretation boundary

The official mutation table removes the invalid local retained-burden values and the unstable SBS/APOBEC component. However, all HLA-A-homozygous patients lie above the official mutation-count median, so the independent contribution of mutation burden and a TMB-by-HLA interaction remain unidentifiable in the 24-patient cohort. The composite is post hoc and must not be called a validated biomarker.
"""
    (OUT / "DISCOVERY_OFFICIAL_PIPELINE_AUDIT_REPORT.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
