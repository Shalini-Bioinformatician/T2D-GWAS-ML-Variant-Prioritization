# =============================================================================
# PHASE 13 — FINAL RESULTS VISUALIZATION & PROJECT OUTPUTS
# =============================================================================

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# CONFIGURATION
# =============================================================================

OUTPUT_DIR = "../phase13_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 90)
print("PHASE 13 — FINAL RESULTS VISUALIZATION & PROJECT OUTPUTS")
print("=" * 90)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def minmax(series):
    """Min-max normalize a numeric pandas Series."""
    series = pd.to_numeric(series, errors="coerce")

    if series.max() == series.min():
        return pd.Series(np.ones(len(series)), index=series.index)

    return (series - series.min()) / (series.max() - series.min())


def save_plot(filename):
    path = os.path.join(OUTPUT_DIR, filename)
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved figure: {path}")


# =============================================================================
# 13.0 — LOCATE INPUT FILES
# =============================================================================

print("\n" + "=" * 90)
print("13.0 — LOCATING INPUT FILES")
print("=" * 90)

ranking_candidates = [
    "../T2D_phase11D_final_candidate_ranking.tsv",
    "T2D_phase11D_final_candidate_ranking.tsv",
    "../phase11C/T2D_phase11D_final_candidate_ranking.tsv"
]

ranking_file = None

for candidate in ranking_candidates:
    if os.path.exists(candidate):
        ranking_file = candidate
        break


if ranking_file is None:
    matches = glob.glob("../**/T2D_phase11D_final_candidate_ranking.tsv",
                        recursive=True)

    if matches:
        ranking_file = matches[0]


if ranking_file is None:
    raise FileNotFoundError(
        "Could not locate T2D_phase11D_final_candidate_ranking.tsv"
    )


print(f"Final ranking file found: {ranking_file}")


# =============================================================================
# LOAD FINAL CANDIDATE RANKING
# =============================================================================

df = pd.read_csv(ranking_file, sep="\t")

print("\nFinal ranking loaded successfully.")
print(f"Rows: {df.shape[0]}")
print(f"Columns: {df.shape[1]}")
print(f"Unique SNPs: {df['Index_SNV'].nunique()}")


# =============================================================================
# 13.1 — FINAL CANDIDATE RANKING PLOT
# =============================================================================

print("\n" + "=" * 90)
print("13.1 — FINAL CANDIDATE RANKING VISUALIZATION")
print("=" * 90)

plot_df = df.sort_values("Final_Rank").copy()

labels = (
    plot_df["Index_SNV"].astype(str)
    + "\n"
    + plot_df["Primary_gene"].astype(str)
)

plt.figure(figsize=(12, 8))

plt.barh(
    labels,
    plot_df["Integrated_Evidence_Score"]
)

plt.xlabel("Integrated Evidence Score")
plt.ylabel("Candidate SNP / Gene")
plt.title("Final Prioritized T2D Candidate Ranking")

plt.gca().invert_yaxis()

save_plot("Figure_1_Final_Candidate_Ranking.png")


# =============================================================================
# 13.2 — INTEGRATED EVIDENCE SCORE COMPONENTS
# =============================================================================

print("\n" + "=" * 90)
print("13.2 — INTEGRATED EVIDENCE SCORE COMPONENTS")
print("=" * 90)

component_df = df.sort_values("Final_Rank").head(10).copy()

x = np.arange(len(component_df))
width = 0.25

ml_component = minmax(
    component_df["Protein_altering_probability"]
)

external_component = minmax(
    component_df["External_GWAS_Evidence_Score"]
)

integrated_component = component_df[
    "Integrated_Evidence_Score"
]

plt.figure(figsize=(14, 8))

plt.bar(
    x - width,
    ml_component,
    width,
    label="ML / Protein-altering Evidence"
)

plt.bar(
    x,
    external_component,
    width,
    label="External GWAS Evidence"
)

plt.bar(
    x + width,
    integrated_component,
    width,
    label="Integrated Evidence"
)

plt.xticks(
    x,
    component_df["Index_SNV"].astype(str)
    + "\n"
    + component_df["Primary_gene"].astype(str),
    rotation=45,
    ha="right"
)

plt.ylabel("Normalized / Integrated Score")
plt.title("Comparison of Evidence Components for Top 10 Candidates")
plt.legend()

save_plot("Figure_2_Integrated_Evidence_Components.png")


# =============================================================================
# 13.3 — EXTERNAL GWAS EVIDENCE PLOT
# =============================================================================

print("\n" + "=" * 90)
print("13.3 — EXTERNAL GWAS EVIDENCE VISUALIZATION")
print("=" * 90)

external_df = df.sort_values(
    "External_GWAS_Evidence_Score",
    ascending=False
).head(10)

labels = (
    external_df["Index_SNV"].astype(str)
    + "\n"
    + external_df["Primary_gene"].astype(str)
)

plt.figure(figsize=(12, 8))

plt.barh(
    labels,
    external_df["External_GWAS_Evidence_Score"]
)

plt.xlabel("External GWAS Evidence Score")
plt.ylabel("Candidate SNP / Gene")
plt.title("Top Candidates Ranked by External GWAS Evidence")

plt.gca().invert_yaxis()

save_plot("Figure_3_External_GWAS_Evidence.png")


# =============================================================================
# 13.4 — MODEL VS BIOLOGICAL EVIDENCE DISCORDANCE
# =============================================================================

print("\n" + "=" * 90)
print("13.4 — MODEL VS BIOLOGICAL EVIDENCE DISCORDANCE")
print("=" * 90)

if "Validation_Group" in df.columns:

    discordance_counts = (
        df["Validation_Group"]
        .fillna("Unknown")
        .value_counts()
    )

    plt.figure(figsize=(10, 7))

    plt.bar(
        discordance_counts.index,
        discordance_counts.values
    )

    plt.xlabel("Validation Group")
    plt.ylabel("Number of Candidates")
    plt.title("Model and Biological Evidence Agreement / Discordance")

    plt.xticks(rotation=30, ha="right")

    save_plot("Figure_4_Model_Biology_Discordance.png")

else:
    print(
        "Validation_Group column not found. "
        "Skipping discordance count plot."
    )


# =============================================================================
# 13.5 — FINAL CANDIDATE EVIDENCE MATRIX / HEATMAP
# =============================================================================

print("\n" + "=" * 90)
print("13.5 — FINAL CANDIDATE EVIDENCE MATRIX")
print("=" * 90)

evidence_columns = [
    "Protein_altering_probability",
    "Direct_T2D_Associations",
    "Glycemic_Associations",
    "Metabolic_Associations",
    "Lipid_Cardiovascular_Associations",
    "External_GWAS_Evidence_Score",
    "Integrated_Evidence_Score"
]

available_columns = [
    col for col in evidence_columns
    if col in df.columns
]

matrix_df = (
    df.sort_values("Final_Rank")
    .set_index("Index_SNV")[available_columns]
    .copy()
)

matrix_normalized = matrix_df.copy()

for col in matrix_normalized.columns:
    matrix_normalized[col] = minmax(matrix_normalized[col])

plt.figure(
    figsize=(
        14,
        max(8, len(matrix_normalized) * 0.5)
    )
)

plt.imshow(
    matrix_normalized.values,
    aspect="auto"
)

plt.colorbar(label="Normalized Evidence Strength")

plt.xticks(
    np.arange(len(matrix_normalized.columns)),
    matrix_normalized.columns,
    rotation=45,
    ha="right"
)

plt.yticks(
    np.arange(len(matrix_normalized.index)),
    matrix_normalized.index
)

plt.xlabel("Evidence Type")
plt.ylabel("Candidate SNP")
plt.title("Final Candidate Evidence Matrix")

save_plot("Figure_5_Final_Evidence_Heatmap.png")


# =============================================================================
# 13.6 — FINAL PRIORITY DISTRIBUTION
# =============================================================================

print("\n" + "=" * 90)
print("13.6 — FINAL PRIORITY DISTRIBUTION")
print("=" * 90)

priority_counts = (
    df["Final_Priority"]
    .value_counts()
)

plt.figure(figsize=(9, 7))

plt.bar(
    priority_counts.index,
    priority_counts.values
)

plt.xlabel("Final Priority Category")
plt.ylabel("Number of Candidates")
plt.title("Distribution of Final Candidate Priorities")

plt.xticks(rotation=20, ha="right")

save_plot("Figure_6_Final_Priority_Distribution.png")


# =============================================================================
# 13.7 — SENSITIVITY ANALYSIS
# =============================================================================

print("\n" + "=" * 90)
print("13.7 — SENSITIVITY ANALYSIS")
print("=" * 90)

rank_columns = [
    col for col in df.columns
    if col.startswith("Rank_ML_")
]

score_columns = [
    col for col in df.columns
    if col.startswith("Score_ML_")
]


# Try locating Phase 12 output files if ranking columns are not present
if len(rank_columns) == 0:

    phase12_files = glob.glob(
        "../**/*phase12*.tsv",
        recursive=True
    )

    phase12_files += glob.glob(
        "../**/*Phase12*.tsv",
        recursive=True
    )

    phase12_files = list(set(phase12_files))

    print("\nSearching for Phase 12 sensitivity output...")

    for file in phase12_files:

        try:
            temp = pd.read_csv(file, sep="\t")

            possible_rank_columns = [
                col for col in temp.columns
                if col.startswith("Rank_ML_")
            ]

            if len(possible_rank_columns) > 0:

                print(f"Sensitivity file found: {file}")

                sensitivity_df = temp.copy()
                rank_columns = possible_rank_columns
                break

        except Exception:
            pass

else:
    sensitivity_df = df.copy()


if "sensitivity_df" in locals() and len(rank_columns) > 0:

    print("Sensitivity rank columns:")
    for col in rank_columns:
        print(f"  {col}")

    top_sensitivity = (
        sensitivity_df.sort_values("Final_Rank")
        .head(10)
        .copy()
    )

    plt.figure(figsize=(14, 8))

    for col in rank_columns:

        plt.plot(
            top_sensitivity["Index_SNV"],
            top_sensitivity[col],
            marker="o",
            label=col.replace("Rank_", "")
        )

    plt.gca().invert_yaxis()

    plt.xlabel("Candidate SNP")
    plt.ylabel("Rank")
    plt.title(
        "Sensitivity Analysis: Candidate Rank Stability "
        "Across ML / External Evidence Weighting Schemes"
    )

    plt.xticks(rotation=45, ha="right")
    plt.legend()

    save_plot("Figure_7_Sensitivity_Rank_Stability.png")

else:

    print(
        "No sensitivity rank columns found automatically. "
        "Skipping sensitivity plot."
    )


# =============================================================================
# 13.8 — CREATE FINAL PROJECT SUMMARY TABLE
# =============================================================================

print("\n" + "=" * 90)
print("13.8 — FINAL PROJECT SUMMARY TABLE")
print("=" * 90)

summary_columns = [
    "Final_Rank",
    "Index_SNV",
    "Primary_gene",
    "Most_severe_consequence",
    "Protein_altering_probability",
    "Has_protein_altering",
    "Predicted_label",
    "OR",
    "neg_log10_P",
    "MAF",
    "Direct_T2D_Associations",
    "Glycemic_Associations",
    "Metabolic_Associations",
    "Lipid_Cardiovascular_Associations",
    "External_GWAS_Evidence_Score",
    "External_Evidence_Tier",
    "Integrated_Evidence_Score",
    "Integrated_Evidence_Tier",
    "Biological_Interpretation",
    "Final_Priority"
]

available_summary_columns = [
    col for col in summary_columns
    if col in df.columns
]

final_summary = (
    df[available_summary_columns]
    .sort_values("Final_Rank")
    .copy()
)

summary_file = os.path.join(
    OUTPUT_DIR,
    "T2D_Final_Project_Candidate_Summary.tsv"
)

final_summary.to_csv(
    summary_file,
    sep="\t",
    index=False
)

print(f"Saved final summary table: {summary_file}")


# =============================================================================
# 13.9 — CREATE TOP 10 CANDIDATE TABLE
# =============================================================================

top10 = final_summary.head(10).copy()

top10_file = os.path.join(
    OUTPUT_DIR,
    "T2D_Top10_Final_Candidates.tsv"
)

top10.to_csv(
    top10_file,
    sep="\t",
    index=False
)

print(f"Saved Top 10 candidate table: {top10_file}")


# =============================================================================
# 13.10 — FINAL QC
# =============================================================================

print("\n" + "=" * 90)
print("13.10 — FINAL QC")
print("=" * 90)

print(f"Total final candidates: {len(df)}")
print(f"Unique SNPs: {df['Index_SNV'].nunique()}")

print(
    f"Missing Final_Rank: "
    f"{df['Final_Rank'].isna().sum()}"
)

print(
    f"Duplicate SNPs: "
    f"{df['Index_SNV'].duplicated().sum()}"
)

print("\nPriority distribution:")

if "Final_Priority" in df.columns:
    print(df["Final_Priority"].value_counts())

print("\nTop 10 candidates:")

display_columns = [
    col for col in [
        "Final_Rank",
        "Index_SNV",
        "Primary_gene",
        "Integrated_Evidence_Score",
        "Final_Priority"
    ]
    if col in df.columns
]

print(
    df.sort_values("Final_Rank")[display_columns]
    .head(10)
    .to_string(index=False)
)


# =============================================================================
# COMPLETE
# =============================================================================

print("\n" + "=" * 90)
print("PHASE 13 COMPLETE")
print("=" * 90)

print("\nGenerated outputs are saved in:")
print(OUTPUT_DIR)

print("\nExpected key outputs:")
print("1. Figure_1_Final_Candidate_Ranking.png")
print("2. Figure_2_Integrated_Evidence_Components.png")
print("3. Figure_3_External_GWAS_Evidence.png")
print("4. Figure_4_Model_Biology_Discordance.png")
print("5. Figure_5_Final_Evidence_Heatmap.png")
print("6. Figure_6_Final_Priority_Distribution.png")
print("7. Figure_7_Sensitivity_Rank_Stability.png")
print("8. T2D_Final_Project_Candidate_Summary.tsv")
print("9. T2D_Top10_Final_Candidates.tsv")

print("\nFINAL PROJECT ANALYSIS PIPELINE COMPLETE.")
