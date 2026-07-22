"""Re-score the paired audit after replacing LEN07 A*34:05 by official A*34:01."""

from pathlib import Path
import os

import numpy as np
import pandas as pd
from mhcflurry import Class1PresentationPredictor


OUT = Path(__file__).resolve().parent
SOURCE = Path(os.environ["RCC_INPUT_AUDIT_LEN07_OFFICIAL_HLA_REPREDICTION_01"])
LEGACY = Path(os.environ["RCC_INPUT_AUDIT_LEN07_OFFICIAL_HLA_REPREDICTION_02"])


def main() -> None:
    data = pd.read_csv(SOURCE)
    official = data.copy()
    mask = official.Patient.eq("LEN07") & official.score_role.eq("partner")
    official.loc[mask, "allele"] = "HLA-A*34:01"
    official.loc[official.Patient.eq("LEN07"), "partner_hla"] = "HLA-A*34:01"
    official.loc[official.Patient.eq("LEN07"), "partner_allele"] = "A3401"

    predictor = Class1PresentationPredictor.load()
    allele_panel = {allele: [allele] for allele in official.allele.drop_duplicates()}
    prediction = predictor.predict(
        peptides=official.peptide.fillna("").astype(str).tolist(),
        n_flanks=official.n_flank.fillna("").astype(str).tolist(),
        c_flanks=official.c_flank.fillna("").astype(str).tolist(),
        alleles=allele_panel,
        sample_names=official.allele.tolist(),
        include_affinity_percentile=True,
        verbose=0,
        throw=False,
    ).reset_index(drop=True)

    score_columns = [
        "affinity",
        "affinity_percentile",
        "processing_score",
        "presentation_score",
        "presentation_percentile",
    ]
    scored = official.reset_index(drop=True).copy()
    for column in score_columns:
        scored[f"mhcflurry_{column}"] = prediction[column]
    scored["mhcflurry_best_allele"] = prediction.best_allele
    scored.to_csv(OUT / "paired_predictions_official_len07_a3401.csv", index=False)

    value_columns = [f"mhcflurry_{column}" for column in score_columns]
    pivot = scored.pivot(index="event_id", columns="score_role", values=value_columns)
    pivot.columns = [f"{metric}_{role}" for metric, role in pivot.columns]
    pivot = pivot.reset_index()
    meta = scored.drop_duplicates("event_id")[[
        "event_id", "Patient", "partner_allele", "partner_hla", "route_class", "peptide"
    ]]
    events = meta.merge(pivot, on="event_id", validate="one_to_one")
    events["partner_better_affinity"] = events.mhcflurry_affinity_partner < events.mhcflurry_affinity_a0201
    events["partner_better_affinity_percentile"] = events.mhcflurry_affinity_percentile_partner < events.mhcflurry_affinity_percentile_a0201
    events["partner_better_presentation"] = events.mhcflurry_presentation_score_partner > events.mhcflurry_presentation_score_a0201
    events["partner_better_presentation_percentile"] = events.mhcflurry_presentation_percentile_partner < events.mhcflurry_presentation_percentile_a0201
    flags = [
        "partner_better_affinity",
        "partner_better_affinity_percentile",
        "partner_better_presentation",
        "partner_better_presentation_percentile",
    ]
    events["partner_better_all"] = events[flags].all(axis=1)
    events["presentation_score_margin_partner"] = events.mhcflurry_presentation_score_partner - events.mhcflurry_presentation_score_a0201
    events.to_csv(OUT / "paired_events_official_len07_a3401.csv", index=False)

    legacy = pd.read_csv(LEGACY)
    legacy_scores = legacy[value_columns].reset_index(drop=True)
    unchanged = ~official.Patient.eq("LEN07")
    max_non_len07_delta = float(np.nanmax(np.abs(scored.loc[unchanged, value_columns].to_numpy() - legacy_scores.loc[unchanged, value_columns].to_numpy())))
    len07 = events.loc[events.Patient.eq("LEN07")]
    report = f"""# LEN07 official-HLA re-prediction audit

## HLA provenance conflict

- OptiType on both LEN07 WES runs: HLA-A*02:01/HLA-A*34:05.
- POLYSOLVER on both runs: HLA-A*34:01/HLA-A*34:01.
- Published Lee et al. Supplementary Data 2: HLA-A*02:01/HLA-A*34:01.
- The legacy paired audit used HLA-A*34:05. This audit substitutes the published HLA-A*34:01 call without changing peptides, flanks, or the comparison allele.

## Re-prediction result

- All paired events: {int(events.partner_better_all.sum())}/{len(events)} satisfy all four partner-better criteria after substitution.
- LEN07 events: {int(len07.partner_better_all.sum())}/{len(len07)} satisfy all four criteria.
- LEN07 median presentation-score margin: {len07.presentation_score_margin_partner.median():.4f}.
- Maximum absolute re-run delta for non-LEN07 score fields versus the archived MHCflurry output: {max_non_len07_delta:.3g}.

## Manuscript consequence

Any retained LEN07 panel must use the published HLA-A*34:01 label and disclose the cross-tool allele discrepancy. The re-prediction addresses prediction-label validity but does not resolve the underlying genotype discordance or turn the case audit into clinical validation.
"""
    (OUT / "LEN07_OFFICIAL_HLA_REPREDICTION_REPORT.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
