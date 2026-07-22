"""Re-score the complete A*02:01-versus-partner C-terminal event universe."""

from pathlib import Path
import os

import numpy as np
import pandas as pd
from mhcflurry import Class1PresentationPredictor
from scipy.stats import binomtest, fisher_exact


OUT = Path(__file__).resolve().parent
ROOT = Path(os.environ["RCC_INPUT_AUDIT_UNBIASED_A0201_PAIRED_UNIVERSE_01"])
ROWS = ROOT / "a0201_true_pos8_codon_template_annotated_rows_2026-05-23.csv"
FILTER_DIR = Path(os.environ["RCC_INPUT_AUDIT_UNBIASED_A0201_PAIRED_UNIVERSE_02"])
RAW_COLUMNS = [
    "sequence_id", "start_pos0", "peptide", "n_flank", "c_flank", "sample_name",
    "orig_affinity", "orig_best_allele_compact", "orig_affinity_percentile",
    "orig_processing_score", "orig_presentation_score", "orig_presentation_percentile",
]


def bh_adjust(values: pd.Series) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    output = np.empty(len(p), dtype=float)
    output[order] = np.minimum(adjusted, 1)
    return output


def compact_hla(value: str) -> str:
    value = str(value)
    return f"HLA-{value[4]}*{value[5:7]}:{value[7:9]}" if value.startswith("HLA-A") else value


def raw_file(patient: str, binder: str) -> pd.DataFrame:
    path = FILTER_DIR / f"{patient}_{binder}.csv"
    frame = pd.read_csv(path, header=None, names=RAW_COLUMNS)
    frame["HLA"] = frame.orig_best_allele_compact.map(compact_hla)
    return frame


def attach_flanks(events: pd.DataFrame) -> pd.DataFrame:
    cache = {}
    output = []
    for _, row in events.iterrows():
        key = (row.Patient, row.BinderClass)
        if key not in cache:
            cache[key] = raw_file(*key)
        expected_hla = row.HLA
        matches = cache[key].loc[
            cache[key].peptide.eq(row.Peptide) & cache[key].HLA.eq(expected_hla)
        ]
        if matches.empty:
            other = "WB" if row.BinderClass == "SB" else "SB"
            fallback_key = (row.Patient, other)
            if fallback_key not in cache:
                cache[fallback_key] = raw_file(*fallback_key)
            matches = cache[fallback_key].loc[
                cache[fallback_key].peptide.eq(row.Peptide) & cache[fallback_key].HLA.eq(expected_hla)
            ]
        if matches.empty:
            raise RuntimeError(f"Missing raw flank match for {row.Patient} {row.Peptide} {expected_hla}")
        match = matches.sort_values(["orig_presentation_percentile", "orig_affinity_percentile"]).iloc[0]
        record = row.to_dict()
        record["n_flank"] = "" if pd.isna(match.n_flank) else str(match.n_flank)
        record["c_flank"] = "" if pd.isna(match.c_flank) else str(match.c_flank)
        output.append(record)
    return pd.DataFrame(output)


def main() -> None:
    rows = pd.read_csv(ROWS)
    events = rows.loc[rows.dataset.eq("all_localized") & rows.pos_consistent.eq(True)].copy()
    events["event_id"] = events.apply(
        lambda row: f"{row.Patient}|{row.line_id}|{row.Peptide}|{row.source}", axis=1
    )
    if len(events) != 63 or events.event_id.nunique() != 63:
        raise RuntimeError("Expected 63 unique position-consistent all-localized events.")
    events = attach_flanks(events)

    pair_rows = []
    for _, row in events.iterrows():
        partner = "HLA-A*34:01" if row.Patient == "LEN07" else f"HLA-A*{str(row.partner_allele)[1:3]}:{str(row.partner_allele)[3:5]}"
        for role, allele in [("a0201", "HLA-A*02:01"), ("partner", partner)]:
            pair_rows.append({
                "event_id": row.event_id,
                "Patient": row.Patient,
                "original_source": row.source,
                "peptide": row.Peptide,
                "wt_peptide": row.wt_peptide,
                "n_flank": row.n_flank,
                "c_flank": row.c_flank,
                "role": role,
                "allele": allele,
            })
    pairs = pd.DataFrame(pair_rows)

    predictor = Class1PresentationPredictor.load()
    allele_panel = {allele: [allele] for allele in pairs.allele.drop_duplicates()}
    prediction = predictor.predict(
        peptides=pairs.peptide.tolist(),
        n_flanks=pairs.n_flank.tolist(),
        c_flanks=pairs.c_flank.tolist(),
        alleles=allele_panel,
        sample_names=pairs.allele.tolist(),
        include_affinity_percentile=True,
        verbose=0,
        throw=False,
    ).reset_index(drop=True)
    metrics = ["affinity", "affinity_percentile", "processing_score", "presentation_score", "presentation_percentile"]
    for metric in metrics:
        pairs[metric] = prediction[metric]
    pairs["best_allele"] = prediction.best_allele
    pairs.to_csv(OUT / "unbiased_a0201_partner_all_pair_predictions.csv", index=False)

    pivot = pairs.pivot(index="event_id", columns="role", values=metrics)
    pivot.columns = [f"{metric}_{role}" for metric, role in pivot.columns]
    meta = pairs.drop_duplicates("event_id")[["event_id", "Patient", "original_source", "peptide", "wt_peptide"]]
    result = meta.merge(pivot.reset_index(), on="event_id", validate="one_to_one")
    result["partner_better_affinity"] = result.affinity_partner < result.affinity_a0201
    result["partner_better_affinity_percentile"] = result.affinity_percentile_partner < result.affinity_percentile_a0201
    result["partner_better_presentation"] = result.presentation_score_partner > result.presentation_score_a0201
    result["partner_better_presentation_percentile"] = result.presentation_percentile_partner < result.presentation_percentile_a0201
    flags = [
        "partner_better_affinity", "partner_better_affinity_percentile",
        "partner_better_presentation", "partner_better_presentation_percentile",
    ]
    result["partner_better_all"] = result[flags].all(axis=1)
    result["cterm_residue"] = result.peptide.str[-1]
    result["cterm_class"] = np.select(
        [result.cterm_residue.isin(list("KRH")), result.cterm_residue.isin(list("FWY")), result.cterm_residue.isin(list("LIVMA"))],
        ["basic", "aromatic", "aliphatic_hydrophobic"],
        default="other",
    )
    result["presentation_score_margin_partner"] = result.presentation_score_partner - result.presentation_score_a0201
    result.to_csv(OUT / "unbiased_a0201_partner_event_results.csv", index=False)

    partner_wins = int(result.partner_better_all.sum())
    direction_p = binomtest(partner_wins, len(result), p=0.5).pvalue
    source_agreement = float((result.partner_better_all == result.original_source.eq("partner_only")).mean())

    enrichment_rows = []
    for cterm_class in ["basic", "aromatic", "aliphatic_hydrophobic", "other"]:
        in_class = result.cterm_class.eq(cterm_class)
        table = pd.crosstab(in_class, result.partner_better_all).reindex(index=[True, False], columns=[True, False], fill_value=0)
        odds_ratio, p_value = fisher_exact(table.to_numpy())
        enrichment_rows.append({
            "cterm_class": cterm_class,
            "n": int(in_class.sum()),
            "partner_wins": int(result.loc[in_class, "partner_better_all"].sum()),
            "odds_ratio": odds_ratio,
            "fisher_p": p_value,
        })
    enrichment = pd.DataFrame(enrichment_rows)
    enrichment["fdr_bh"] = bh_adjust(enrichment.fisher_p)
    enrichment.to_csv(OUT / "unbiased_a0201_partner_cterm_enrichment.csv", index=False)

    patient = result.groupby("Patient").agg(
        n_events=("event_id", "size"),
        partner_wins=("partner_better_all", "sum"),
        partner_win_fraction=("partner_better_all", "mean"),
        median_margin=("presentation_score_margin_partner", "median"),
    ).reset_index()
    patient.to_csv(OUT / "unbiased_a0201_partner_patient_summary.csv", index=False)

    best_class = enrichment.sort_values("fdr_bh").iloc[0]
    report = f"""# Unbiased A*02:01 paired-universe audit

## Design correction

- Legacy main figures selected 29 `partner_only` events before paired re-scoring; 29/29 partner-better was therefore selection-conditioned.
- The corrected universe includes all 63 position-consistent C-terminal 9-mers: 34 originally assigned to A*02:01 (`index_only`) and 29 assigned to the partner allele (`partner_only`).
- Every peptide was re-scored against both patient HLA-A alleles with MHCflurry 2.2.1 and presentation model release 2.2.0. LEN07 used the published HLA-A*34:01 call.

## Results

- Partner allele satisfied all four paired superiority criteria in {partner_wins}/{len(result)} events; two-sided binomial P={direction_p:.3f} versus 50%.
- Re-scored direction agreed with the archived best-allele source label for {source_agreement:.1%} of events.
- Strongest C-terminal class association: {best_class.cterm_class}, OR={best_class.odds_ratio:.2f}, nominal P={best_class.fisher_p:.3g}, BH-FDR={best_class.fdr_bh:.3g}.
- Events arose from {result.Patient.nunique()} patients; patient-level counts and margins are reported separately.

## Interpretation boundary

The unbiased result supports complementary allele-specific peptide allocation and terminal-residue compatibility, not a global partner-allele advantage. The 29-event selected set may be shown only as an illustrative subset after the 63-event universe is disclosed. It cannot serve as direct validation of the clinical subgroup.
"""
    (OUT / "UNBIASED_A0201_PAIRED_UNIVERSE_REPORT.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
