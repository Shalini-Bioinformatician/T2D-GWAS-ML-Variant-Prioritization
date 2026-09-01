import pandas as pd
import numpy as np
from sklearn.metrics import (
    confusion_matrix,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    balanced_accuracy_score
)

print("=" * 80)
print("PHASE 12 — FINAL MODEL VALIDATION & ROBUSTNESS ANALYSIS")
print("=" * 80)

# ============================================================
# 12.0 LOAD FINAL PHASE 11D DATA
# ============================================================

input_file = "../../phase11C/GWAS_Catalog/T2D_phase11D_final_candidate_ranking.tsv"

print("\nLoading final candidate ranking...")
df = pd.read_csv(input_file, sep="\t")

print(f"Rows: {len(df)}")
print(f"Unique SNPs: {df['Index_SNV'].nunique()}")

# ============================================================
# 12.0 INPUT QC
# ============================================================

required = [
    "Index_SNV",
    "Primary_gene",
    "Protein_altering_probability",
    "Has_protein_altering",
    "Predicted_label",
    "Direct_T2D_Associations",
    "External_GWAS_Evidence_Score",
    "Integrated_Evidence_Score",
    "Final_Rank",
    "Final_Priority"
]

missing = [c for c in required if c not in df.columns]

if missing:
    raise ValueError(f"Missing required columns: {missing}")

if df["Index_SNV"].duplicated().any():
    raise ValueError("Duplicate SNPs detected.")

print("\n" + "=" * 80)
print("12.0 INPUT QC")
print("=" * 80)

print("QC PASS")
print(f"SNP count: {len(df)}")
print(f"Unique SNPs: {df['Index_SNV'].nunique()}")
print(f"Missing required columns: {len(missing)}")

# ============================================================
# 12.1 DESCRIPTIVE MODEL PERFORMANCE
# ============================================================

print("\n" + "=" * 80)
print("12.1 MODEL PERFORMANCE — DESCRIPTIVE VALIDATION PANEL")
print("=" * 80)

y_true = df["Has_protein_altering"].astype(int)
y_pred = df["Predicted_label"].astype(int)

cm = confusion_matrix(y_true, y_pred, labels=[0, 1])

tn, fp, fn, tp = cm.ravel()

accuracy = accuracy_score(y_true, y_pred)
precision = precision_score(y_true, y_pred, zero_division=0)
recall = recall_score(y_true, y_pred, zero_division=0)
f1 = f1_score(y_true, y_pred, zero_division=0)
balanced_acc = balanced_accuracy_score(y_true, y_pred)

specificity = tn / (tn + fp) if (tn + fp) > 0 else np.nan

print("\nConfusion Matrix")
print("                 Predicted 0    Predicted 1")
print(f"Actual 0             {tn:3d}            {fp:3d}")
print(f"Actual 1             {fn:3d}            {tp:3d}")

print("\nPerformance metrics")
print(f"Accuracy:           {accuracy:.4f}")
print(f"Precision:          {precision:.4f}")
print(f"Recall/Sensitivity: {recall:.4f}")
print(f"Specificity:        {specificity:.4f}")
print(f"F1 score:           {f1:.4f}")
print(f"Balanced accuracy:  {balanced_acc:.4f}")

metrics = pd.DataFrame({
    "Metric": [
        "Accuracy",
        "Precision",
        "Recall",
        "Specificity",
        "F1",
        "Balanced_Accuracy"
    ],
    "Value": [
        accuracy,
        precision,
        recall,
        specificity,
        f1,
        balanced_acc
    ]
})

metrics.to_csv(
    "phase12_model_performance_metrics.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 12.2 PROBABILITY DISTRIBUTION BY BIOLOGICAL LABEL
# ============================================================

print("\n" + "=" * 80)
print("12.2 MODEL PROBABILITY DISTRIBUTION")
print("=" * 80)

prob_summary = (
    df.groupby("Final_Priority")
      .agg(
          N=("Index_SNV", "count"),
          Mean_Probability=("Protein_altering_probability", "mean"),
          Median_Probability=("Protein_altering_probability", "median"),
          Min_Probability=("Protein_altering_probability", "min"),
          Max_Probability=("Protein_altering_probability", "max")
      )
      .reset_index()
)

print(prob_summary.to_string(index=False))

prob_summary.to_csv(
    "phase12_probability_by_priority.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 12.3 ERROR ANALYSIS
# ============================================================

print("\n" + "=" * 80)
print("12.3 MODEL ERROR ANALYSIS")
print("=" * 80)

def classify_error(row):
    if row["Has_protein_altering"] == 1 and row["Predicted_label"] == 0:
        return "False_Negative"
    elif row["Has_protein_altering"] == 0 and row["Predicted_label"] == 1:
        return "False_Positive"
    elif row["Has_protein_altering"] == 1 and row["Predicted_label"] == 1:
        return "True_Positive"
    else:
        return "True_Negative"

df["Model_Error_Type"] = df.apply(classify_error, axis=1)

error_cols = [
    "Index_SNV",
    "Primary_gene",
    "Most_severe_consequence",
    "Protein_altering_probability",
    "Has_protein_altering",
    "Predicted_label",
    "Direct_T2D_Associations",
    "External_GWAS_Evidence_Score",
    "Integrated_Evidence_Score",
    "Final_Priority",
    "Model_Error_Type"
]

errors = df[df["Model_Error_Type"].isin(
    ["False_Positive", "False_Negative"]
)][error_cols].copy()

print("\nDiscordant candidates:")
print(errors.to_string(index=False))

errors.to_csv(
    "phase12_model_error_analysis.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 12.4 SENSITIVITY ANALYSIS
# ============================================================

print("\n" + "=" * 80)
print("12.4 INTEGRATED SCORE SENSITIVITY ANALYSIS")
print("=" * 80)

# Normalize external evidence to 0-1
ext_min = df["External_GWAS_Evidence_Score"].min()
ext_max = df["External_GWAS_Evidence_Score"].max()

if ext_max == ext_min:
    df["External_Normalized"] = 0.0
else:
    df["External_Normalized"] = (
        df["External_GWAS_Evidence_Score"] - ext_min
    ) / (ext_max - ext_min)

df["ML_Normalized"] = df["Protein_altering_probability"]

weight_schemes = {
    "ML_30_External_70": (0.30, 0.70),
    "ML_40_External_60": (0.40, 0.60),
    "ML_50_External_50": (0.50, 0.50),
    "ML_60_External_40": (0.60, 0.40),
    "ML_70_External_30": (0.70, 0.30)
}

rank_tables = []

for scheme, (ml_weight, ext_weight) in weight_schemes.items():

    score_col = f"Score_{scheme}"
    rank_col = f"Rank_{scheme}"

    df[score_col] = (
        ml_weight * df["ML_Normalized"] +
        ext_weight * df["External_Normalized"]
    )

    df[rank_col] = (
        df[score_col]
        .rank(method="min", ascending=False)
        .astype(int)
    )

    tmp = df[
        [
            "Index_SNV",
            "Primary_gene",
            score_col,
            rank_col
        ]
    ].copy()

    tmp["Weighting_Scheme"] = scheme

    rank_tables.append(tmp)

    print(f"\n{scheme}")
    print(
        tmp.sort_values(rank_col)
           .head(10)
           .to_string(index=False)
    )

sensitivity_long = pd.concat(rank_tables, ignore_index=True)

sensitivity_long.to_csv(
    "phase12_sensitivity_rankings.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 12.5 RANK STABILITY
# ============================================================

print("\n" + "=" * 80)
print("12.5 RANK STABILITY ANALYSIS")
print("=" * 80)

rank_cols = [f"Rank_{x}" for x in weight_schemes.keys()]

stability = df[
    [
        "Index_SNV",
        "Primary_gene",
        "Final_Rank"
    ] + rank_cols
].copy()

stability["Mean_Sensitivity_Rank"] = stability[rank_cols].mean(axis=1)
stability["Median_Sensitivity_Rank"] = stability[rank_cols].median(axis=1)
stability["Min_Sensitivity_Rank"] = stability[rank_cols].min(axis=1)
stability["Max_Sensitivity_Rank"] = stability[rank_cols].max(axis=1)
stability["Rank_SD"] = stability[rank_cols].std(axis=1)

# Smaller SD = more stable
def stability_class(sd):
    if sd <= 1.0:
        return "HIGHLY_STABLE"
    elif sd <= 2.5:
        return "STABLE"
    elif sd <= 4.0:
        return "MODERATELY_STABLE"
    else:
        return "SENSITIVE_RANKING"

stability["Rank_Stability"] = stability["Rank_SD"].apply(
    stability_class
)

stability = stability.sort_values(
    ["Mean_Sensitivity_Rank", "Rank_SD"]
)

print(
    stability[
        [
            "Index_SNV",
            "Primary_gene",
            "Final_Rank",
            "Mean_Sensitivity_Rank",
            "Min_Sensitivity_Rank",
            "Max_Sensitivity_Rank",
            "Rank_SD",
            "Rank_Stability"
        ]
    ].to_string(index=False)
)

stability.to_csv(
    "phase12_rank_stability.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 12.6 ROBUSTNESS CLASSIFICATION
# ============================================================

print("\n" + "=" * 80)
print("12.6 FINAL ROBUSTNESS CLASSIFICATION")
print("=" * 80)

merged = df.merge(
    stability[
        [
            "Index_SNV",
            "Mean_Sensitivity_Rank",
            "Min_Sensitivity_Rank",
            "Max_Sensitivity_Rank",
            "Rank_SD",
            "Rank_Stability"
        ]
    ],
    on="Index_SNV",
    how="left"
)

def robustness_label(row):

    if (
        row["Mean_Sensitivity_Rank"] <= 3
        and row["Rank_SD"] <= 2.5
    ):
        return "ROBUST_TOP_CANDIDATE"

    elif (
        row["Mean_Sensitivity_Rank"] <= 8
        and row["Rank_SD"] <= 4.0
    ):
        return "ROBUST_STRONG_CANDIDATE"

    elif row["Rank_SD"] <= 4.0:
        return "MODERATELY_STABLE"

    else:
        return "SENSITIVE_RANKING"

merged["Robustness_Class"] = merged.apply(
    robustness_label,
    axis=1
)

final_cols = [
    "Index_SNV",
    "Primary_gene",
    "Final_Rank",
    "Final_Priority",
    "Integrated_Evidence_Score",
    "Protein_altering_probability",
    "External_GWAS_Evidence_Score",
    "Mean_Sensitivity_Rank",
    "Min_Sensitivity_Rank",
    "Max_Sensitivity_Rank",
    "Rank_SD",
    "Rank_Stability",
    "Robustness_Class"
]

final_robustness = merged[
    final_cols
].sort_values(
    ["Mean_Sensitivity_Rank", "Rank_SD"]
)

print(
    final_robustness.to_string(index=False)
)

final_robustness.to_csv(
    "phase12_final_robustness_classification.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 12.7 TOP ROBUST CANDIDATES
# ============================================================

print("\n" + "=" * 80)
print("12.7 TOP ROBUST CANDIDATES")
print("=" * 80)

top_robust = final_robustness[
    final_robustness["Robustness_Class"].isin(
        [
            "ROBUST_TOP_CANDIDATE",
            "ROBUST_STRONG_CANDIDATE"
        ]
    )
].copy()

print(
    top_robust.to_string(index=False)
)

top_robust.to_csv(
    "phase12_top_robust_candidates.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 12.8 FINAL QC
# ============================================================

print("\n" + "=" * 80)
print("12.8 FINAL QC")
print("=" * 80)

checks = {
    "17 candidates retained": len(df) == 17,
    "17 unique SNPs": df["Index_SNV"].nunique() == 17,
    "No missing SNPs": df["Index_SNV"].notna().all(),
    "No duplicate SNPs": not df["Index_SNV"].duplicated().any(),
    "All sensitivity rankings generated":
        all(df[c].notna().all() for c in rank_cols),
    "All robustness classifications generated":
        merged["Robustness_Class"].notna().all()
}

for name, result in checks.items():
    print(f"{name}: {'PASS' if result else 'FAIL'}")

if not all(checks.values()):
    raise ValueError("One or more Phase 12 QC checks failed.")

print("\n" + "=" * 80)
print("PHASE 12 COMPLETE")
print("=" * 80)

print("\nSAVED FILES")
print("phase12_model_performance_metrics.tsv")
print("phase12_probability_by_priority.tsv")
print("phase12_model_error_analysis.tsv")
print("phase12_sensitivity_rankings.tsv")
print("phase12_rank_stability.tsv")
print("phase12_final_robustness_classification.tsv")
print("phase12_top_robust_candidates.tsv")
