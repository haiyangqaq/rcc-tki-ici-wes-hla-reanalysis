from __future__ import annotations

import shutil
from itertools import combinations
from pathlib import Path
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import restricted_mean_survival_time


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "submission_package_v1_8_due_diligence"
FIG = OUT / "figures"
SUPP = OUT / "supplementary_figures"
SOURCE = OUT / "source_data"
for directory in [OUT, FIG, SUPP, SOURCE]:
    directory.mkdir(parents=True, exist_ok=True)

DISCOVERY = ROOT / "discovery_official_maf_hla_patient_level.csv"
MIAO = Path(os.environ["RCC_INPUT_BUILD_V18_FIGURES_01"])
MSK = ROOT / "msk138_patient_level_audit_source.csv"
PAIRED = ROOT / "unbiased_a0201_partner_event_results.csv"
PAIR_PRED = ROOT / "unbiased_a0201_partner_all_pair_predictions.csv"
PAIR_ENRICH = ROOT / "unbiased_a0201_partner_cterm_enrichment.csv"
HLA_AUDIT = ROOT / "discovery_hla_version_comparison.csv"
SIGNATURE = ROOT / "signature_full_strict_bootstrap_summary.csv"
PASS_COUNTS = ROOT / "pass_snv_counts.csv"
OFFICIAL_MAF = Path(os.environ["RCC_INPUT_BUILD_V18_FIGURES_02"])
PXD_SOURCE = Path(os.environ["RCC_INPUT_BUILD_V18_FIGURES_03"])
PXD_TAILS = PXD_SOURCE / "supplementary_figure_S5_pxd017149_logo_tail_source_fontplus2_2026-06-22_2026-06-22.csv"
PXD_CLASSES = PXD_SOURCE / "supplementary_figure_S5_pxd017149_ctail_class_source_fontplus2_2026-06-22_2026-06-22.csv"

COLORS = {
    "teal": "#007C83",
    "orange": "#D95F02",
    "blue": "#3B6FB6",
    "gold": "#C89B00",
    "green": "#3A923A",
    "purple": "#7A5195",
    "red": "#C23B3B",
    "gray": "#73777A",
    "light": "#E8ECEF",
    "dark": "#202124",
}


def set_style() -> None:
    sns.set_theme(style="white", context="talk")
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 13,
        "axes.labelsize": 15,
        "axes.titlesize": 16,
        "axes.titleweight": "bold",
        "axes.labelweight": "bold",
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 11,
        "figure.dpi": 160,
        "savefig.dpi": 300,
        "axes.linewidth": 1.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def panel(ax, label: str) -> None:
    ax.text(-0.10, 1.06, label, transform=ax.transAxes, fontsize=20, fontweight="bold", va="top")


def save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def km_step(ax, data, duration, event, group, palette, xlim, title=None):
    for name, subset in data.groupby(group, observed=False):
        km = KaplanMeierFitter().fit(subset[duration], subset[event], label=str(name))
        km.plot_survival_function(ax=ax, ci_show=False, color=palette[str(name)], linewidth=2.6, show_censors=True,
                                  censor_styles={"marker": "|", "ms": 8, "mew": 1.8})
    ax.set(xlim=xlim, ylim=(0, 1.04), xlabel="Time (months)", ylabel="Progression-free survival")
    if title:
        ax.set_title(title, loc="left")
    ax.legend(frameon=False, loc="lower left")
    sns.despine(ax=ax)


def risk_table(ax, data, duration, event, group, times, palette):
    levels = list(data[group].drop_duplicates())
    ax.set_xlim(min(times), max(times))
    ax.set_ylim(-0.6, len(levels) - 0.4)
    for y, level in enumerate(levels[::-1]):
        subset = data.loc[data[group].eq(level)]
        for t in times:
            ax.text(t, y, str(int((subset[duration] >= t).sum())), ha="center", va="center", fontsize=10,
                    color=palette[str(level)], fontweight="bold")
        ax.text(min(times) - (max(times) - min(times)) * 0.03, y, str(level), ha="right", va="center", fontsize=10)
    ax.set_xticks(times)
    ax.set_yticks([])
    ax.set_xlabel("Time (months)")
    ax.set_title("Number at risk", loc="left", fontsize=11)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)


def fig1(discovery: pd.DataFrame) -> None:
    order = discovery.sort_values("official_protein_altering_count", ascending=False).reset_index(drop=True)
    x = np.arange(len(order))
    core = order.HLA_A_homozygous.astype(bool)
    fig = plt.figure(figsize=(17, 12))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.75, 1.25], hspace=0.35, wspace=0.36)

    ax = fig.add_subplot(gs[0, 0])
    panel(ax, "A")
    ax.axis("off")
    boxes = [
        (0.01, 0.48, 0.19, 0.34, "Published\nRCC cohort\nn = 24"),
        (0.27, 0.48, 0.19, 0.34, "Official MAF,\nHLA and PFS\n16/8"),
        (0.53, 0.48, 0.19, 0.34, "External RCC\ntreated cohorts\nlinked outcomes"),
        (0.79, 0.48, 0.19, 0.34, "Unbiased paired\npeptide audit\n63 events"),
    ]
    for i, (x0, y0, w, h, text) in enumerate(boxes):
        color = [COLORS["teal"], COLORS["blue"], COLORS["orange"], COLORS["purple"]][i]
        ax.add_patch(mpl.patches.FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0.015,rounding_size=0.015",
                                                facecolor="white", edgecolor=color, linewidth=2.2))
        ax.text(x0 + w / 2, y0 + h / 2, text, ha="center", va="center", fontsize=10.5, fontweight="bold")
        if i < len(boxes) - 1:
            ax.annotate("", xy=(boxes[i + 1][0] - 0.01, y0 + h / 2), xytext=(x0 + w + 0.01, y0 + h / 2),
                        arrowprops=dict(arrowstyle="-|>", color=COLORS["gray"], lw=1.8))
    ax.text(0.01, 0.25, "Aim: distinguish reproducible associations from treatment-context effects,\nand global allele rescue from terminal-residue compatibility.",
            fontsize=10.5, va="top")
    ax.set_title("Audit design", loc="left")

    ax = fig.add_subplot(gs[0, 1])
    panel(ax, "B")
    colors = np.where(core, COLORS["orange"], COLORS["teal"])
    ax.bar(x, order.official_protein_altering_count, color=colors, width=0.78)
    ax.axhline(order.official_protein_altering_count.median(), color=COLORS["dark"], linestyle="--", linewidth=1.5,
               label=f"Median = {order.official_protein_altering_count.median():.1f}")
    for i in np.where(core)[0]:
        ax.text(i, order.loc[i, "official_protein_altering_count"] + 2.5, order.loc[i, "Patient"], rotation=90,
                ha="center", va="bottom", fontsize=10, fontweight="bold", color=COLORS["orange"])
    ax.set_xticks(x)
    ax.set_xticklabels(order.Patient, rotation=90)
    ax.set_ylabel("Protein-altering mutation count")
    ax.set_title("Published mutation calls and HLA-A zygosity", loc="left")
    handles = [mpl.patches.Patch(color=COLORS["orange"], label="HLA-A homozygous"),
               mpl.patches.Patch(color=COLORS["teal"], label="HLA-A heterozygous")]
    ax.legend(handles=handles + [mpl.lines.Line2D([], [], color=COLORS["dark"], ls="--", label="Cohort median")],
              frameon=False, ncol=2, loc="upper right")
    sns.despine(ax=ax)

    ax = fig.add_subplot(gs[1, 0])
    panel(ax, "C")
    matrix = order[["HLA_A_homozygous", "HLA_B_homozygous", "HLA_C_homozygous"]].T.astype(int)
    sns.heatmap(matrix, cmap=mpl.colors.ListedColormap(["#F2F4F5", COLORS["orange"]]), cbar=False,
                linewidths=0.7, linecolor="white", ax=ax)
    ax.set_yticklabels(["HLA-A homozygous", "HLA-B homozygous", "HLA-C homozygous"], rotation=0)
    ax.set_xticklabels(order.Patient, rotation=90)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Locus-specific HLA-I zygosity", loc="left")

    ax = fig.add_subplot(gs[1, 1])
    panel(ax, "D")
    timeline = discovery.sort_values("PFS", ascending=True).reset_index(drop=True)
    y = np.arange(len(timeline))
    c = np.where(timeline.HLA_A_homozygous, COLORS["orange"], COLORS["teal"])
    ax.hlines(y, 0, timeline.PFS, color=c, linewidth=3)
    events = timeline.PFS_status.eq(1)
    ax.scatter(timeline.loc[events, "PFS"], y[events], marker="o", s=54, color=c[events], edgecolor="white", linewidth=0.8, zorder=3)
    ax.scatter(timeline.loc[~events, "PFS"], y[~events], marker=">", s=70, color=c[~events], edgecolor="white", linewidth=0.8, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(timeline.Patient)
    ax.set_xlabel("Progression-free survival (months)")
    ax.set_title("Patient-level PFS with censoring", loc="left")
    ax.legend(handles=[mpl.lines.Line2D([], [], marker="o", ls="", color=COLORS["dark"], label="Progression event"),
                       mpl.lines.Line2D([], [], marker=">", ls="", color=COLORS["dark"], label="Censored")],
              frameon=False, loc="lower right")
    sns.despine(ax=ax)
    save(fig, FIG / "Figure1.png")
    order.to_csv(SOURCE / "figure1_discovery_landscape.csv", index=False)


def rmst(frame, flag, tau=12):
    a = KaplanMeierFitter().fit(frame.loc[flag, "PFS"], frame.loc[flag, "PFS_status"])
    b = KaplanMeierFitter().fit(frame.loc[~flag, "PFS"], frame.loc[~flag, "PFS_status"])
    return float(restricted_mean_survival_time(a, t=tau) - restricted_mean_survival_time(b, t=tau))


def fig2(discovery: pd.DataFrame) -> None:
    data = discovery.copy()
    data["State"] = np.where(data.HLA_A_homozygous, "HLA-A homozygous (n=3)", "Other patients (n=21)")
    palette = {"HLA-A homozygous (n=3)": COLORS["orange"], "Other patients (n=21)": COLORS["teal"]}
    fig = plt.figure(figsize=(17, 13))
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.30)

    ax = fig.add_subplot(gs[0, 0])
    panel(ax, "A")
    sns.boxplot(data=data, x="HLA_A_homozygous", y="official_protein_altering_count", color="white", width=0.48,
                showfliers=False, ax=ax)
    sns.stripplot(data=data, x="HLA_A_homozygous", y="official_protein_altering_count", hue="HLA_A_homozygous",
                  palette={False: COLORS["teal"], True: COLORS["orange"]}, size=8, jitter=0.15, ax=ax, legend=False)
    ax.axhline(data.official_protein_altering_count.median(), ls="--", lw=1.4, color=COLORS["gray"])
    ax.set_xticklabels(["Heterozygous\n(n = 21)", "Homozygous\n(n = 3)"])
    ax.set_xlabel("HLA-A zygosity")
    ax.set_ylabel("Official protein-altering mutation count")
    ax.set_title("All HLA-A-homozygous cases are above the burden median", loc="left")
    sns.despine(ax=ax)

    inner = gs[0, 1].subgridspec(2, 1, height_ratios=[4.0, 1.1], hspace=0.08)
    ax = fig.add_subplot(inner[0])
    panel(ax, "B")
    km_step(ax, data, "PFS", "PFS_status", "State", palette, (0, 42), "Post hoc discovery-state PFS")
    state = data.HLA_A_homozygous
    p = logrank_test(data.loc[state, "PFS"], data.loc[~state, "PFS"],
                     event_observed_A=data.loc[state, "PFS_status"], event_observed_B=data.loc[~state, "PFS_status"]).p_value
    ax.text(0.98, 0.94, f"Log-rank P = {p:.4f}\n12-month RMST difference = {rmst(data, state):.2f} months",
            transform=ax.transAxes, ha="right", va="top", fontsize=11,
            bbox=dict(facecolor="white", edgecolor=COLORS["light"], pad=5))
    ax2 = fig.add_subplot(inner[1], sharex=ax)
    risk_table(ax2, data, "PFS", "PFS_status", "State", [0, 10, 20, 30, 40], palette)
    ax.tick_params(labelbottom=False)
    ax.set_xlabel("")

    ax = fig.add_subplot(gs[1, 0])
    panel(ax, "C")
    rows = []
    definitions = {
        "HLA-A homozygous": data.HLA_A_homozygous.astype(bool),
        "High burden + any HLA homo\n(median cutoff)": data.high_burden_any_hla_homo_median.astype(bool),
        "High burden + any HLA homo\n(upper-tertile cutoff)": data.high_burden_any_hla_homo_tertile.astype(bool),
        "Any HLA-I homozygous": data.HLA_any_homozygous.astype(bool),
    }
    for name, flag in definitions.items():
        lr = logrank_test(data.loc[flag, "PFS"], data.loc[~flag, "PFS"],
                          event_observed_A=data.loc[flag, "PFS_status"], event_observed_B=data.loc[~flag, "PFS_status"])
        rows.append({"definition": name, "n": int(flag.sum()), "events": int(data.loc[flag, "PFS_status"].sum()),
                     "rmst_difference": rmst(data, flag), "logrank_p": lr.p_value})
    sens = pd.DataFrame(rows)
    yy = np.arange(len(sens))
    ax.axvline(0, color=COLORS["gray"], lw=1.2)
    ax.scatter(sens.rmst_difference, yy, s=85, color=[COLORS["orange"], COLORS["orange"], COLORS["gold"], COLORS["blue"]])
    for i, row in sens.iterrows():
        ax.text(row.rmst_difference + 0.18, i, f"n={row.n}; P={row.logrank_p:.3g}", va="center", fontsize=11)
    ax.set_yticks(yy)
    ax.set_yticklabels(sens.definition)
    ax.invert_yaxis()
    ax.set_xlabel("12-month RMST difference vs remainder (months)")
    ax.set_title("Definition and threshold sensitivity", loc="left")
    sns.despine(ax=ax)

    ax = fig.add_subplot(gs[1, 1])
    panel(ax, "D")
    table = pd.crosstab(data.burden_high_median.map({False: "Low", True: "High"}),
                        data.HLA_A_homozygous.map({False: "Heterozygous", True: "Homozygous"})).reindex(
                            index=["Low", "High"], columns=["Heterozygous", "Homozygous"], fill_value=0)
    sns.heatmap(table, annot=True, fmt="d", cmap=mpl.colors.LinearSegmentedColormap.from_list("cells", ["#FFFFFF", "#B8D8D8"]),
                cbar=False, linewidths=2, linecolor="white", annot_kws={"fontsize": 22, "fontweight": "bold"}, ax=ax)
    ax.set_xlabel("HLA-A zygosity")
    ax.set_ylabel("Official mutation burden")
    ax.set_title("The low-burden/homozygous cell is empty", loc="left")
    ax.text(0.5, -0.25, "Mutation and HLA-A effects cannot be separated in n = 24;\nno interaction estimate is identifiable.",
            transform=ax.transAxes, ha="center", fontsize=12, color=COLORS["red"], fontweight="bold")
    save(fig, FIG / "Figure2.png")
    sens.to_csv(SOURCE / "figure2_definition_sensitivity.csv", index=False)


def fig3(miao: pd.DataFrame) -> None:
    d = miao.copy()
    d["PFS_EVENT"] = pd.to_numeric(d.PFS_EVENT)
    d["mutation_count"] = d.protein_altering_mutation_count
    d = d.sort_values("mutation_count", ascending=False).reset_index(drop=True)
    fig = plt.figure(figsize=(17, 13))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.28)

    ax = fig.add_subplot(gs[0, 0])
    panel(ax, "A")
    x = np.arange(len(d))
    ax.bar(x, d.mutation_count, color=np.where(d.PFS_EVENT.eq(1), COLORS["orange"], COLORS["teal"]), width=0.78)
    ax.set_xticks(x)
    ax.set_xticklabels(d.patient_id, rotation=90)
    ax.set_ylabel("Protein-altering mutation count")
    ax.set_title("Miao/DFCI nivolumab cohort (n = 35)", loc="left")
    ax.legend(handles=[mpl.patches.Patch(color=COLORS["orange"], label="PFS event"),
                       mpl.patches.Patch(color=COLORS["teal"], label="Censored")], frameon=False)
    sns.despine(ax=ax)

    ax = fig.add_subplot(gs[0, 1])
    panel(ax, "B")
    models = pd.read_csv(ROOT / "miao_robustness_cox_models.csv")
    models["label"] = models.model.map({
        "univariable_rank_per_iqr": "Rank per IQR",
        "univariable_count_per_iqr": "Count per 34-mutation IQR",
        "adjusted_first_line_and_log2_dose": "Rank per IQR, adjusted",
    })
    y = np.arange(len(models))
    ax.axvline(1, color=COLORS["gray"], ls="--", lw=1.3)
    ax.errorbar(models.hazard_ratio, y,
                xerr=[models.hazard_ratio - models.lower_95ci, models.upper_95ci - models.hazard_ratio],
                fmt="o", color=COLORS["orange"], ecolor=COLORS["dark"], capsize=4, markersize=8)
    for i, row in models.iterrows():
        ax.text(row.upper_95ci * 1.04, i, f"P={row.wald_p:.3f}", va="center", fontsize=11)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(models.label)
    ax.invert_yaxis()
    ax.set_xlabel("Hazard ratio for shorter PFS (95% CI)")
    ax.set_title("Continuous mutation burden association", loc="left")
    sns.despine(ax=ax)

    ax = fig.add_subplot(gs[1, 0])
    panel(ax, "C")
    loo = pd.read_csv(ROOT / "miao_leave_one_out_influence.csv").sort_values("hazard_ratio")
    y = np.arange(len(loo))
    ax.axvline(1, color=COLORS["gray"], ls="--", lw=1.2)
    ax.hlines(y, loo.lower_95ci, loo.upper_95ci, color="#AAB2B8", linewidth=1)
    ax.scatter(loo.hazard_ratio, y, c=np.where(loo.wald_p.lt(0.05), COLORS["orange"], COLORS["gray"]), s=38)
    ax.set_yticks(y[::3])
    ax.set_yticklabels(loo.excluded_patient.iloc[::3])
    ax.set_xscale("log")
    ax.set_xlabel("Leave-one-out HR per mutation-rank IQR")
    ax.set_ylabel("Excluded patient")
    ax.set_title("Direction retained in all 35 exclusions", loc="left")
    sns.despine(ax=ax)

    inner = gs[1, 1].subgridspec(2, 1, height_ratios=[4.0, 1.1], hspace=0.08)
    ax = fig.add_subplot(inner[0])
    panel(ax, "D")
    median = d.mutation_count.median()
    d["Burden"] = np.where(d.mutation_count.ge(median), "At/above median", "Below median")
    pal = {"At/above median": COLORS["orange"], "Below median": COLORS["teal"]}
    km_step(ax, d, "PFS_MONTHS", "PFS_EVENT", "Burden", pal, (0, 48), "Median split is descriptive")
    flag = d.Burden.eq("At/above median")
    p = logrank_test(d.loc[flag, "PFS_MONTHS"], d.loc[~flag, "PFS_MONTHS"],
                     event_observed_A=d.loc[flag, "PFS_EVENT"], event_observed_B=d.loc[~flag, "PFS_EVENT"]).p_value
    ax.text(0.97, 0.94, f"Log-rank P = {p:.3f}", transform=ax.transAxes, ha="right", va="top", fontsize=11)
    ax2 = fig.add_subplot(inner[1], sharex=ax)
    risk_table(ax2, d, "PFS_MONTHS", "PFS_EVENT", "Burden", [0, 12, 24, 36, 48], pal)
    ax.tick_params(labelbottom=False)
    ax.set_xlabel("")
    save(fig, FIG / "Figure3.png")
    d.to_csv(SOURCE / "figure3_miao_patient_level.csv", index=False)
    models.to_csv(SOURCE / "figure3_miao_cox_models.csv", index=False)
    loo.to_csv(SOURCE / "figure3_miao_leave_one_out.csv", index=False)


def fig4(msk: pd.DataFrame) -> None:
    d = msk.copy()
    high = d.loc[d.is_ccrcc.astype(bool) & d.TMB_high_median.eq(1)].copy()
    high["HLA-I zygosity"] = np.where(high.HLA_any_homozygous.eq(1), "Any locus homozygous", "All loci heterozygous")
    pal = {"Any locus homozygous": COLORS["orange"], "All loci heterozygous": COLORS["teal"]}
    fig = plt.figure(figsize=(17, 13))
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.30)

    ax = fig.add_subplot(gs[0, 0])
    panel(ax, "A")
    for label, subset in high.groupby("HLA-I zygosity"):
        km = KaplanMeierFitter().fit(subset.OS_MONTHS, subset.OS_EVENT, label=f"{label} (n={len(subset)})")
        km.plot_survival_function(ax=ax, ci_show=False, color=pal[label], linewidth=2.6, show_censors=True,
                                  censor_styles={"marker": "|", "ms": 8, "mew": 1.8})
    ax.set(xlim=(0, 70), ylim=(0, 1.04), xlabel="Overall survival from ICI start (months)", ylabel="Overall survival")
    ax.set_title("MSK ccRCC, high-TMB subgroup", loc="left")
    ax.legend(frameon=False, loc="lower left")
    sns.despine(ax=ax)

    ax = fig.add_subplot(gs[0, 1])
    panel(ax, "B")
    cond = pd.read_csv(ROOT / "msk138_high_tmb_conditional_hla_tests.csv")
    view = cond.loc[(cond.cohort.eq("ccRCC")) & (cond.tmb_threshold.eq("TMB_high_median"))].copy()
    label_map = {"HLA_A_homozygous": "HLA-A homozygous", "HLA_any_homozygous": "Any HLA-I homozygous",
                 "HLA_A_LOH": "HLA-A LOH", "HLA_any_LOH": "Any HLA-I LOH", "HED_A_zero": "HED-A = 0"}
    view = view.reset_index(drop=True)
    view["label"] = view.apply(lambda row: f"{label_map[row.defect]}\nFDR={row.cox_fdr_bh:.3f}", axis=1)
    y = np.arange(len(view))
    ax.axvline(1, color=COLORS["gray"], ls="--", lw=1.2)
    ax.errorbar(view.hazard_ratio, y, xerr=[view.hazard_ratio-view.lower_95ci, view.upper_95ci-view.hazard_ratio],
                fmt="o", color=COLORS["orange"], ecolor=COLORS["dark"], capsize=4, markersize=8)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(view.label)
    ax.invert_yaxis()
    ax.set_xlabel("Adjusted OS hazard ratio (95% CI)")
    ax.set_title("Conditional tests within high TMB", loc="left")
    sns.despine(ax=ax)

    ax = fig.add_subplot(gs[1, 0])
    panel(ax, "C")
    inter = pd.read_csv(ROOT / "msk138_tmb_hla_interaction_models.csv")
    inter = inter.loc[inter.analysis.str.startswith("binary")].copy()
    inter = inter.reset_index(drop=True)
    inter["label"] = inter.apply(lambda row: f"{label_map[row.analysis.split(':')[1]]}\nFDR={row.fdr_bh:.3f}", axis=1)
    y = np.arange(len(inter))
    ax.axvline(1, color=COLORS["gray"], ls="--", lw=1.2)
    ax.errorbar(inter.hazard_ratio, y, xerr=[inter.hazard_ratio-inter.lower_95ci, inter.upper_95ci-inter.hazard_ratio],
                fmt="o", color=COLORS["blue"], ecolor=COLORS["dark"], capsize=4, markersize=8)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(inter.label)
    ax.invert_yaxis()
    ax.set_xlabel("TMB-high × HLA interaction HR (95% CI)")
    ax.set_title("No interaction survives multiplicity correction", loc="left")
    sns.despine(ax=ax)

    ax = fig.add_subplot(gs[1, 1])
    panel(ax, "D")
    evidence = pd.DataFrame([
        ["Lee lenvatinib+pembrolizumab", 24, "PFS", "Post hoc", "Adverse", "Composite = HLA-A homozygosity"],
        ["Miao nivolumab", 35, "PFS", "Mutation only", "Adverse", "No linked HLA"],
        ["Braun nivolumab", 261, "PFS/OS", "Direct", "Null", "TMB, neoantigen and zygosity null"],
        ["JAVELIN avelumab+axitinib", 333, "PFS", "Direct proxy", "Null", "High-TMB HLA-A single vs multi"],
        ["Chowell multi-ICI RCC", 91, "PFS/OS", "Direct proxy", "Null", "High-TMB/low-HED null"],
        ["SNiP-RCC nivolumab", 222, "Benefit/CSS", "HLA only", "Protective", "High HED-B associated with benefit"],
        ["Nature Cancer IO/IO", 128, "Response", "HLA only", "Null", "Mean HED ER vs PD null"],
        ["MSK RCC ICI", 138, "OS", "Post hoc", "Adverse", "Conditional FDR; interaction FDR null"],
    ], columns=["cohort", "n", "endpoint", "directness", "direction", "boundary"])
    color_map = {"Adverse": COLORS["orange"], "Null": COLORS["gray"], "Protective": COLORS["green"]}
    y = np.arange(len(evidence))
    ax.scatter(np.zeros(len(evidence)), y, s=130, color=evidence.direction.map(color_map), edgecolor="white", linewidth=1)
    for i, row in evidence.iterrows():
        ax.text(-0.08, i, row.cohort, ha="right", va="center", fontsize=10.5, fontweight="bold")
        ax.text(0.10, i, f"n={row.n}; {row.endpoint}; {row.boundary}", ha="left", va="center", fontsize=9.7)
    ax.set_xlim(-0.75, 1.75)
    ax.set_ylim(len(evidence)-0.4, -0.6)
    ax.axis("off")
    ax.set_title("Cross-cohort evidence is heterogeneous", loc="left")
    ax.legend(handles=[mpl.lines.Line2D([], [], marker="o", ls="", color=v, label=k, markersize=9)
                       for k, v in color_map.items()], frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.08))
    save(fig, FIG / "Figure4.png")
    high.to_csv(SOURCE / "figure4_msk_high_tmb_ccrcc.csv", index=False)
    view.to_csv(SOURCE / "figure4_msk_conditional_tests.csv", index=False)
    inter.to_csv(SOURCE / "figure4_msk_interaction_tests.csv", index=False)
    evidence.to_csv(SOURCE / "figure4_cross_cohort_evidence_matrix.csv", index=False)


def fig5(events: pd.DataFrame) -> None:
    d = events.copy()
    class_order = ["basic", "aromatic", "aliphatic_hydrophobic", "other"]
    d["cterm_class"] = pd.Categorical(d.cterm_class, class_order, ordered=True)
    d = d.sort_values(["cterm_class", "presentation_score_margin_partner"], ascending=[True, False]).reset_index(drop=True)
    fig = plt.figure(figsize=(17, 13))
    gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.30)

    ax = fig.add_subplot(gs[0, 0])
    panel(ax, "A")
    x = np.arange(len(d))
    bar_colors = np.where(d.partner_better_all, COLORS["orange"], COLORS["teal"])
    ax.bar(x, d.presentation_score_margin_partner, color=bar_colors, width=0.82)
    ax.axhline(0, color=COLORS["dark"], lw=1)
    ax.set_xticks(x[::3])
    ax.set_xticklabels(d.cterm_residue.iloc[::3], rotation=0)
    ax.set_xlabel("C-terminal residue (every third event labelled)")
    ax.set_ylabel("Presentation score: partner - A*02:01")
    ax.set_title("Complete 63-event paired universe", loc="left")
    ax.legend(handles=[mpl.patches.Patch(color=COLORS["orange"], label="Partner wins all four metrics"),
                       mpl.patches.Patch(color=COLORS["teal"], label="A*02:01/complementary allocation")], frameon=False)
    sns.despine(ax=ax)

    ax = fig.add_subplot(gs[0, 1])
    panel(ax, "B")
    sns.scatterplot(data=d, x="presentation_score_a0201", y="presentation_score_partner", hue="cterm_class",
                    palette={"basic": COLORS["orange"], "aromatic": COLORS["purple"],
                             "aliphatic_hydrophobic": COLORS["teal"], "other": COLORS["gray"]},
                    style="partner_better_all", s=90, ax=ax, legend=False)
    lim = [min(d.presentation_score_a0201.min(), d.presentation_score_partner.min()) - 0.05,
           max(d.presentation_score_a0201.max(), d.presentation_score_partner.max()) + 0.05]
    ax.plot(lim, lim, ls="--", color=COLORS["dark"], lw=1.2)
    ax.set(xlim=lim, ylim=lim, xlabel="A*02:01 presentation score", ylabel="Partner-allele presentation score")
    ax.set_title("Terminal class separates allele preference", loc="left")
    handles = [mpl.lines.Line2D([], [], marker="o", ls="", color=color, label=label, markersize=8)
               for label, color in [("Basic", COLORS["orange"]), ("Aromatic", COLORS["purple"]),
                                    ("Aliphatic hydrophobic", COLORS["teal"]), ("Other", COLORS["gray"])]]
    handles += [mpl.lines.Line2D([], [], marker="x", ls="", color=COLORS["dark"], label="Partner wins", markersize=8),
                mpl.lines.Line2D([], [], marker="o", ls="", color=COLORS["dark"], label="A*02:01/complementary", markersize=7)]
    ax.legend(handles=handles, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", title="C-terminal class")
    sns.despine(ax=ax)

    ax = fig.add_subplot(gs[1, 0])
    panel(ax, "C")
    enr = pd.read_csv(PAIR_ENRICH)
    enr["fraction"] = enr.partner_wins / enr.n
    labels = {"basic": "Basic (K/R/H)", "aromatic": "Aromatic (F/W/Y)",
              "aliphatic_hydrophobic": "Aliphatic hydrophobic", "other": "Other"}
    enr["label"] = enr.cterm_class.map(labels)
    y = np.arange(len(enr))
    ax.barh(y, enr.fraction, color=[COLORS["orange"], COLORS["purple"], COLORS["teal"], COLORS["gray"]])
    for i, row in enr.iterrows():
        ax.text(min(row.fraction + 0.03, 0.98), i, f"{int(row.partner_wins)}/{int(row.n)}; FDR={row.fdr_bh:.3g}",
                va="center", fontsize=11, color=COLORS["dark"])
    ax.set_yticks(y)
    ax.set_yticklabels(enr.label)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Fraction with partner winning all four metrics")
    ax.set_title("C-terminal compatibility, not global partner rescue", loc="left")
    sns.despine(ax=ax)

    ax = fig.add_subplot(gs[1, 1])
    panel(ax, "D")
    patient = d.groupby("Patient").agg(n_events=("event_id", "size"), partner_wins=("partner_better_all", "sum"),
                                        median_margin=("presentation_score_margin_partner", "median")).reset_index()
    patient = patient.sort_values("n_events", ascending=True)
    y = np.arange(len(patient))
    ax.barh(y, patient.n_events, color=COLORS["light"], label="All paired events")
    ax.barh(y, patient.partner_wins, color=COLORS["orange"], label="Partner wins")
    ax.set_yticks(y)
    ax.set_yticklabels(patient.Patient)
    ax.set_xlabel("Number of events")
    ax.set_title("Events span seven patients", loc="left")
    ax.legend(frameon=False, loc="lower right")
    sns.despine(ax=ax)
    save(fig, FIG / "Figure5.png")
    d.to_csv(SOURCE / "figure5_complete_paired_event_universe.csv", index=False)
    enr.to_csv(SOURCE / "figure5_cterm_enrichment.csv", index=False)
    patient.to_csv(SOURCE / "figure5_patient_summary.csv", index=False)


def supp1_oncoprint(discovery: pd.DataFrame) -> None:
    maf = pd.read_excel(OFFICIAL_MAF, sheet_name="MAF exome sequencing")
    maf["Patient"] = maf.Tumor_Sample_Barcode.astype(str).str.replace(r"_T$", "", regex=True)
    freq = maf.groupby("Hugo_Symbol").Patient.nunique().sort_values(ascending=False)
    canonical = ["VHL", "PBRM1", "SETD2", "BAP1", "KDM5C", "MTOR", "PTEN", "TP53", "TSC1", "TSC2", "PIK3CA", "ARID1A"]
    genes = []
    for gene in canonical + freq.index.tolist():
        if gene in set(maf.Hugo_Symbol) and gene not in genes:
            genes.append(gene)
        if len(genes) == 20:
            break
    patients = discovery.sort_values("official_protein_altering_count", ascending=False).Patient.tolist()
    classes = ["Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Other"]
    cmap = {"Missense_Mutation": COLORS["blue"], "Nonsense_Mutation": COLORS["red"], "Frame_Shift_Del": COLORS["orange"],
            "Frame_Shift_Ins": COLORS["gold"], "In_Frame_Del": COLORS["green"], "In_Frame_Ins": COLORS["purple"],
            "Splice_Site": "#8C564B", "Other": COLORS["gray"]}
    mat = pd.DataFrame("", index=genes, columns=patients)
    for (gene, patient), group in maf.loc[maf.Hugo_Symbol.isin(genes)].groupby(["Hugo_Symbol", "Patient"]):
        cls = group.Variant_Classification.iloc[0]
        mat.loc[gene, patient] = cls if cls in classes else "Other"
    fig, ax = plt.subplots(figsize=(17, 10))
    for iy, gene in enumerate(genes):
        for ix, patient in enumerate(patients):
            value = mat.loc[gene, patient]
            ax.add_patch(mpl.patches.Rectangle((ix, iy), 1, 1, facecolor=cmap.get(value, "#F1F3F4") if value else "#F1F3F4",
                                               edgecolor="white", linewidth=0.8))
    ax.set_xlim(0, len(patients)); ax.set_ylim(len(genes), 0)
    ax.set_xticks(np.arange(len(patients)) + 0.5); ax.set_xticklabels(patients, rotation=90)
    ax.set_yticks(np.arange(len(genes)) + 0.5); ax.set_yticklabels(genes, fontweight="bold")
    ax.set_title("Supplementary Figure S1 | Published RCC mutation landscape (20 selected recurrent/canonical genes)", loc="left")
    ax.legend(handles=[mpl.patches.Patch(color=v, label=k.replace("_", " ")) for k, v in cmap.items()],
              frameon=False, ncol=4, bbox_to_anchor=(0.5, -0.16), loc="upper center")
    for spine in ax.spines.values(): spine.set_visible(False)
    save(fig, SUPP / "Supplementary_Figure_S1.png")
    mat.to_csv(SOURCE / "supplementary_figure_S1_official_oncoprint.csv")


def supp2_hla() -> None:
    d = pd.read_csv(HLA_AUDIT)
    fig, axes = plt.subplots(1, 2, figsize=(17, 8), gridspec_kw={"width_ratios": [1.2, 1]})
    ax = axes[0]; panel(ax, "A")
    m = d.set_index("Patient")[["A_match", "B_match", "C_match"]].T.astype(int)
    sns.heatmap(m, cmap=mpl.colors.ListedColormap([COLORS["red"], COLORS["green"]]), cbar=False,
                linewidths=0.8, linecolor="white", ax=ax)
    ax.set_yticklabels(["HLA-A", "HLA-B", "HLA-C"], rotation=0)
    ax.set_title("Published vs internal genotype agreement", loc="left")
    ax.set_xlabel(""); ax.set_ylabel("")
    ax = axes[1]; panel(ax, "B")
    diff = d.sort_values("mean_HED_difference_analysis_minus_official")
    ax.axvline(0, color=COLORS["gray"], lw=1.2)
    ax.scatter(diff.mean_HED_difference_analysis_minus_official, np.arange(len(diff)),
               color=np.where(diff.mean_HED_difference_analysis_minus_official.abs().gt(0.1), COLORS["orange"], COLORS["teal"]), s=55)
    ax.set_yticks(np.arange(len(diff))); ax.set_yticklabels(diff.Patient)
    ax.set_xlabel("Historical analysis mean HED - published mean HED")
    ax.set_title("HED version differences", loc="left")
    sns.despine(ax=ax)
    save(fig, SUPP / "Supplementary_Figure_S2.png")
    d.to_csv(SOURCE / "supplementary_figure_S2_hla_version_audit.csv", index=False)


def supp3_sbs() -> None:
    s = pd.read_csv(SIGNATURE)
    c = pd.read_csv(PASS_COUNTS)
    fig, axes = plt.subplots(1, 2, figsize=(17, 7))
    ax = axes[0]; panel(ax, "A")
    sample_col = "Patient" if "Patient" in c.columns else "sample"
    c = c.sort_values("pass_snv_count")
    value_col = "pass_snv_count"
    ax.bar(np.arange(len(c)), c[value_col], color=COLORS["teal"])
    ax.set_xticks(np.arange(len(c))); ax.set_xticklabels(c[sample_col], rotation=90)
    ax.set_ylabel("PASS biallelic SNVs")
    ax.set_title("Filtered Mutect2 input", loc="left")
    sns.despine(ax=ax)
    ax = axes[1]; panel(ax, "B")
    if "signature" in s.columns:
        ap = s.loc[s.signature.isin(["SBS2", "SBS13"])].copy()
        metric = "bootstrap_detection_rate" if "bootstrap_detection_rate" in ap.columns else ap.select_dtypes("number").columns[-1]
        sns.barplot(data=ap, x="Patient", y=metric, hue="signature", palette=[COLORS["orange"], COLORS["purple"]], ax=ax)
        ax.set_ylabel("Bootstrap detection rate")
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=90)
    else:
        ap = s[["sample", "SBS2_detection_fraction", "APOBEC_detection_fraction"]].melt(
            id_vars="sample", var_name="signature", value_name="bootstrap_detection_rate"
        )
        ap["signature"] = ap.signature.map({"SBS2_detection_fraction": "SBS2", "APOBEC_detection_fraction": "SBS2+SBS13"})
        sns.barplot(data=ap, x="sample", y="bootstrap_detection_rate", hue="signature",
                    palette=[COLORS["orange"], COLORS["purple"]], ax=ax)
        ax.set_ylabel("Bootstrap detection rate")
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=90)
    ax.set_title("APOBEC signal fails strict stability testing", loc="left")
    if ax.get_legend(): ax.legend(frameon=False)
    save(fig, SUPP / "Supplementary_Figure_S3.png")
    s.to_csv(SOURCE / "supplementary_figure_S3_signature_stability.csv", index=False)
    c.to_csv(SOURCE / "supplementary_figure_S3_pass_snv_counts.csv", index=False)


def supp4_pair_metrics(events: pd.DataFrame) -> None:
    metrics = ["affinity", "affinity_percentile", "presentation_score", "presentation_percentile"]
    fig, axes = plt.subplots(2, 2, figsize=(17, 13))
    for i, (ax, metric) in enumerate(zip(axes.flat, metrics)):
        panel(ax, chr(65 + i))
        x = events[f"{metric}_a0201"]
        y = events[f"{metric}_partner"]
        sns.scatterplot(x=x, y=y, hue=events.cterm_class, palette={"basic": COLORS["orange"], "aromatic": COLORS["purple"],
                        "aliphatic_hydrophobic": COLORS["teal"], "other": COLORS["gray"]}, s=75, ax=ax, legend=i == 0)
        lim = [min(x.min(), y.min()), max(x.max(), y.max())]
        ax.plot(lim, lim, ls="--", color=COLORS["dark"], lw=1)
        ax.set_xlabel(f"A*02:01 {metric.replace('_', ' ')}")
        ax.set_ylabel(f"Partner {metric.replace('_', ' ')}")
        ax.set_title(metric.replace("_", " ").title(), loc="left")
        if ax.get_legend(): ax.legend(frameon=False)
        sns.despine(ax=ax)
    save(fig, SUPP / "Supplementary_Figure_S4.png")
    events.to_csv(SOURCE / "supplementary_figure_S4_complete_pair_metrics.csv", index=False)


def supp5_pxd() -> None:
    tails = pd.read_csv(PXD_TAILS)
    classes = pd.read_csv(PXD_CLASSES)
    classes = classes.loc[~classes.set_short.str.contains("Current selected", na=False)].copy()
    set_order = ["PXD017149 all HLA-I ligands", "PXD017149 tumor-exclusive HLA-I ligands"]
    aa_order = list("ACDEFGHIKLMNPQRSTVWY")
    matrices = {}
    for name in set_order:
        seqs = tails.loc[tails["set"].eq(name), "cterm_tail"].dropna().astype(str)
        matrix = pd.DataFrame(0.0, index=aa_order, columns=["-5", "-4", "-3", "-2", "-1"])
        for seq in seqs:
            if len(seq) != 5:
                continue
            for pos, aa in enumerate(seq):
                if aa in matrix.index:
                    matrix.iloc[matrix.index.get_loc(aa), pos] += 1
        matrix = matrix.div(matrix.sum(axis=0), axis=1)
        matrices[name] = matrix
    fig = plt.figure(figsize=(17, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.55], hspace=0.38, wspace=0.28)
    for i, name in enumerate(set_order):
        ax = fig.add_subplot(gs[0, i]); panel(ax, chr(65+i))
        sns.heatmap(matrices[name], cmap="viridis", vmin=0, vmax=max(x.values.max() for x in matrices.values()),
                    cbar=i == 1, cbar_kws={"label": "Residue frequency"}, ax=ax)
        ax.set_xlabel("Position from peptide C terminus")
        ax.set_ylabel("Amino acid")
        n = int(tails.loc[tails["set"].eq(name)].shape[0])
        ax.set_title(f"{name}\n(n = {n} ligands)", loc="left")
    ax = fig.add_subplot(gs[1, :]); panel(ax, "C")
    pivot = classes.pivot(index="set_short", columns="class", values="frac")
    pivot = pivot.reindex(columns=["Basic K/R/H", "Aromatic F/Y/W", "Other"])
    pivot.index = [x.replace("\n", " ") for x in pivot.index]
    left = np.zeros(len(pivot))
    for cls, color in zip(pivot.columns, [COLORS["orange"], COLORS["gold"], COLORS["gray"]]):
        values = pivot[cls].to_numpy()
        ax.barh(np.arange(len(pivot)), values, left=left, color=color, label=cls)
        for y, (lft, val) in enumerate(zip(left, values)):
            ax.text(lft + val / 2, y, f"{val:.1%}", ha="center", va="center", color="white", fontsize=11, fontweight="bold")
        left += values
    ax.set_yticks(np.arange(len(pivot))); ax.set_yticklabels(pivot.index)
    ax.set_xlim(0, 1); ax.set_xlabel("Fraction of ligands")
    ax.set_title("C-terminal residue-class composition", loc="left")
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.48))
    sns.despine(ax=ax)
    fig.suptitle("Supplementary Figure S5 | PXD017149 ccRCC HLA-ligandome context (not clinical validation)",
                 x=0.07, ha="left", fontsize=18, fontweight="bold")
    save(fig, SUPP / "Supplementary_Figure_S5.png")
    tails.loc[tails["set"].isin(set_order)].to_csv(SOURCE / "supplementary_figure_S5_pxd017149_tails.csv", index=False)
    classes.to_csv(SOURCE / "supplementary_figure_S5_pxd017149_classes.csv", index=False)


def supp6_msk() -> None:
    cond = pd.read_csv(ROOT / "msk138_high_tmb_conditional_hla_tests.csv")
    inter = pd.read_csv(ROOT / "msk138_tmb_hla_interaction_models.csv")
    fig, axes = plt.subplots(1, 2, figsize=(17, 8), gridspec_kw={"wspace": 0.42})
    ax = axes[0]; panel(ax, "A")
    plot = cond.copy(); plot["key"] = plot.cohort + " | " + plot.tmb_threshold.str.replace("TMB_high_", "") + " | " + plot.defect
    plot = plot.sort_values("hazard_ratio")
    y = np.arange(len(plot))
    ax.axvline(1, color=COLORS["gray"], ls="--")
    ax.hlines(y, plot.lower_95ci, plot.upper_95ci, color="#B0B7BC", lw=1)
    ax.scatter(plot.hazard_ratio, y, c=np.where(plot.cox_fdr_bh.lt(0.05), COLORS["orange"], COLORS["gray"]), s=42)
    ax.set_xscale("log"); ax.set_yticks(y); ax.set_yticklabels(plot.key, fontsize=8.5)
    ax.set_xlabel("Conditional OS HR (95% CI)"); ax.set_title("Twenty conditional tests", loc="left")
    sns.despine(ax=ax)
    ax = axes[1]; panel(ax, "B")
    plot = inter.sort_values("hazard_ratio")
    short = {
        "binary_interaction:HLA_any_homozygous": "Binary: any-HLA homo",
        "binary_interaction:HED_A_zero": "Binary: HED-A=0",
        "continuous_interaction:HLA_A_LOH": "Continuous: HLA-A LOH",
        "continuous_interaction:HLA_any_LOH": "Continuous: any-HLA LOH",
        "binary_interaction:HLA_A_homozygous": "Binary: HLA-A homo",
        "continuous_interaction:HLA_any_homozygous": "Continuous: any-HLA homo",
        "binary_interaction:HLA_A_LOH": "Binary: HLA-A LOH",
        "continuous_interaction:HED_A_zero": "Continuous: HED-A=0",
        "binary_interaction:HLA_any_LOH": "Binary: any-HLA LOH",
        "continuous_interaction:HLA_A_homozygous": "Continuous: HLA-A homo",
    }
    y = np.arange(len(plot))
    ax.axvline(1, color=COLORS["gray"], ls="--")
    ax.hlines(y, plot.lower_95ci, plot.upper_95ci, color="#B0B7BC", lw=1)
    ax.scatter(plot.hazard_ratio, y, c=COLORS["blue"], s=42)
    ax.set_xscale("log"); ax.set_yticks(y); ax.set_yticklabels(plot.analysis.map(short), fontsize=9)
    ax.set_xlabel("TMB × HLA interaction HR (95% CI)"); ax.set_title("Ten interaction tests; all FDR > 0.05", loc="left")
    sns.despine(ax=ax)
    save(fig, SUPP / "Supplementary_Figure_S6.png")
    cond.to_csv(SOURCE / "supplementary_figure_S6_msk_conditional_tests.csv", index=False)
    inter.to_csv(SOURCE / "supplementary_figure_S6_msk_interactions.csv", index=False)


def main() -> None:
    set_style()
    discovery = pd.read_csv(DISCOVERY)
    miao = pd.read_csv(MIAO)
    msk = pd.read_csv(MSK)
    paired = pd.read_csv(PAIRED)
    fig1(discovery)
    fig2(discovery)
    fig3(miao)
    fig4(msk)
    fig5(paired)
    supp1_oncoprint(discovery)
    supp2_hla()
    supp3_sbs()
    supp4_pair_metrics(paired)
    supp5_pxd()
    supp6_msk()
    print(f"Figures written to {OUT}")


if __name__ == "__main__":
    main()
