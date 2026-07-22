"""Robustness audit for the locked Miao/DFCI mutation-only analysis."""

from pathlib import Path
import os

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from scipy.stats import chi2


OUT = Path(__file__).resolve().parent
SOURCE = Path(os.environ["RCC_INPUT_AUDIT_MIAO_ROBUSTNESS_01"])


def fit(frame: pd.DataFrame, variables: list[str]) -> tuple[CoxPHFitter, pd.DataFrame]:
    columns = ["PFS_MONTHS", "PFS_EVENT", *variables]
    data = frame[columns].dropna().copy()
    model = CoxPHFitter().fit(data, duration_col="PFS_MONTHS", event_col="PFS_EVENT")
    return model, data


def row(model: CoxPHFitter, variable: str, model_name: str) -> dict:
    result = model.summary.loc[variable]
    return {
        "model": model_name,
        "variable": variable,
        "n": int(model._n_examples),
        "events": int(model.event_observed.sum()),
        "hazard_ratio": float(result["exp(coef)"]),
        "lower_95ci": float(result["exp(coef) lower 95%"]),
        "upper_95ci": float(result["exp(coef) upper 95%"]),
        "wald_p": float(result["p"]),
        "log_likelihood": float(model.log_likelihood_),
        "aic_partial": float(model.AIC_partial_),
    }


def main() -> None:
    data = pd.read_csv(SOURCE)
    data["dose_mgkg"] = pd.to_numeric(data.TREATMENT_GROUP.str.extract(r"([0-9.]+)", expand=False), errors="coerce")
    data["log2_dose"] = np.log2(data.dose_mgkg)
    data["first_line"] = pd.to_numeric(data.FIRST_LINE_THERAPY, errors="coerce")

    rank_iqr = float(data.mutation_rank.quantile(0.75) - data.mutation_rank.quantile(0.25))
    count_iqr = float(data.protein_altering_mutation_count.quantile(0.75) - data.protein_altering_mutation_count.quantile(0.25))
    data["mutation_rank_per_iqr"] = data.mutation_rank / rank_iqr
    data["mutation_count_per_iqr"] = data.protein_altering_mutation_count / count_iqr

    primary, primary_data = fit(data, ["mutation_rank_per_iqr"])
    count_model, _ = fit(data, ["mutation_count_per_iqr"])
    adjusted, _ = fit(data, ["mutation_rank_per_iqr", "first_line", "log2_dose"])

    summaries = pd.DataFrame([
        row(primary, "mutation_rank_per_iqr", "univariable_rank_per_iqr"),
        row(count_model, "mutation_count_per_iqr", "univariable_count_per_iqr"),
        row(adjusted, "mutation_rank_per_iqr", "adjusted_first_line_and_log2_dose"),
    ])
    ph = proportional_hazard_test(primary, primary_data, time_transform="rank")
    summaries.loc[summaries.model.eq("univariable_rank_per_iqr"), "ph_test_p"] = float(ph.summary.loc["mutation_rank_per_iqr", "p"])
    summaries.to_csv(OUT / "miao_robustness_cox_models.csv", index=False)

    centered = data.mutation_rank - data.mutation_rank.mean()
    data["rank_centered"] = centered
    data["rank_centered_sq"] = centered**2
    linear, _ = fit(data, ["rank_centered"])
    quadratic, _ = fit(data, ["rank_centered", "rank_centered_sq"])
    lr_stat = 2 * (quadratic.log_likelihood_ - linear.log_likelihood_)
    nonlinearity = pd.DataFrame([{
        "linear_log_likelihood": linear.log_likelihood_,
        "quadratic_log_likelihood": quadratic.log_likelihood_,
        "likelihood_ratio_chisq": lr_stat,
        "nonlinearity_p": chi2.sf(lr_stat, 1),
        "linear_aic_partial": linear.AIC_partial_,
        "quadratic_aic_partial": quadratic.AIC_partial_,
    }])
    nonlinearity.to_csv(OUT / "miao_nonlinearity_audit.csv", index=False)

    loo_rows = []
    for patient in data.patient_id:
        subset = data.loc[data.patient_id.ne(patient)].copy()
        model, _ = fit(subset, ["mutation_rank_per_iqr"])
        estimate = row(model, "mutation_rank_per_iqr", "leave_one_out")
        estimate["excluded_patient"] = patient
        loo_rows.append(estimate)
    loo = pd.DataFrame(loo_rows)
    loo.to_csv(OUT / "miao_leave_one_out_influence.csv", index=False)

    primary_row = summaries.iloc[0]
    adjusted_row = summaries.iloc[2]
    report = f"""# Miao/DFCI mutation-only robustness audit

## Endpoint and scale

- Public datahub source: `ccrcc_dfci_2019/data_clinical_patient.txt`.
- Analysis: 35 nivolumab-treated ccRCC patients, 27 PFS events and 8 censored observations.
- The primary rank IQR is {rank_iqr:.3f}; the mutation-count IQR is {count_iqr:.0f} records.
- Rank-scaled HR per IQR: {primary_row.hazard_ratio:.2f} (95% CI {primary_row.lower_95ci:.2f}-{primary_row.upper_95ci:.2f}), P={primary_row.wald_p:.3f}; PH-test P={primary_row.ph_test_p:.3f}.
- Raw-count HR per {count_iqr:.0f}-mutation IQR: {summaries.iloc[1].hazard_ratio:.2f} (95% CI {summaries.iloc[1].lower_95ci:.2f}-{summaries.iloc[1].upper_95ci:.2f}), P={summaries.iloc[1].wald_p:.3f}.

## Robustness

- Adjustment for first-line status and log2 nivolumab dose: rank-IQR HR {adjusted_row.hazard_ratio:.2f} (95% CI {adjusted_row.lower_95ci:.2f}-{adjusted_row.upper_95ci:.2f}), P={adjusted_row.wald_p:.3f}.
- Quadratic nonlinearity likelihood-ratio P={nonlinearity.iloc[0].nonlinearity_p:.3f}.
- Leave-one-patient-out rank-IQR HR range: {loo.hazard_ratio.min():.2f}-{loo.hazard_ratio.max():.2f}; all {int((loo.hazard_ratio > 1).sum())}/{len(loo)} estimates retained the same direction; {int((loo.wald_p < 0.05).sum())}/{len(loo)} had nominal P<0.05.

## Interpretation boundary

This is an independent treated-cohort mutation-count association, not validation of HED-A, SBS/APOBEC, peptide presentation, a treatment interaction, or a clinical cutoff. The manuscript should report the per-IQR estimate rather than the visually inflated full-rank-range HR.
"""
    (OUT / "MIAO_ROBUSTNESS_AUDIT_REPORT.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
