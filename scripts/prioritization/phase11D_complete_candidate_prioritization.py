import pandas as pd
import numpy as np

print("=" * 80)
print("PHASE 11D — COMPLETE CANDIDATE INTERPRETATION & FINAL PRIORITIZATION")
print("=" * 80)

# ============================================================
# FILE PATHS
# ============================================================

integrated_file = "phase11C_integrated_candidate_evidence.tsv"
ranking_file = "phase11C_integrated_candidate_ranking.tsv"

# ============================================================
# LOAD DATA
# ============================================================

print("\nLoading integrated candidate evidence...")

df = pd.read_csv(
    integrated_file,
    sep="\t"
)

print("Rows:", len(df))
print("Unique SNPs:", df["Index_SNV"].nunique())

# Ranking file is loaded as an independent QC source
print("\nLoading integrated ranking...")

ranking = pd.read_csv(
    ranking_file,
    sep="\t"
)

print("Rows:", len(ranking))
print("Unique SNPs:", ranking["Index_SNV"].nunique())

# ============================================================
# 11D.0 — INPUT QC
# ============================================================

print("\n" + "=" * 80)
print("11D.0 — INPUT QC")
print("=" * 80)

if df["Index_SNV"].duplicated().any():
    raise ValueError("Duplicate SNPs detected in integrated evidence.")

if ranking["Index_SNV"].duplicated().any():
    raise ValueError("Duplicate SNPs detected in integrated ranking.")

if set(df["Index_SNV"]) != set(ranking["Index_SNV"]):
    raise ValueError(
        "Mismatch between integrated evidence and ranking SNP sets."
    )

required_columns = [
    "Index_SNV",
    "Primary_gene",
    "Most_severe_consequence",
    "Validation_Group",
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
    "Integrated_Evidence_Tier"
]

missing = [
    c for c in required_columns
    if c not in df.columns
]

if missing:
    raise ValueError(
        "Missing required columns: " + str(missing)
    )

print("QC PASS")
print("SNP count:", df["Index_SNV"].nunique())
print("Required columns present:", len(required_columns))

# ============================================================
# 11D.1 — FINAL CANDIDATE SHORTLIST
# ============================================================

print("\n" + "=" * 80)
print("11D.1 — FINAL CANDIDATE SHORTLIST")
print("=" * 80)

final_df = df.copy()

# Existing integrated score from Phase 11C.4 is retained.
final_df = final_df.sort_values(
    "Integrated_Evidence_Score",
    ascending=False
).reset_index(drop=True)

final_df["Final_Rank"] = range(
    1,
    len(final_df) + 1
)

# Priority based on the already established integrated evidence tier.
tier_map = {
    "Tier_1_High_Priority": "HIGH_PRIORITY",
    "Tier_2_Strong_Candidate": "STRONG_CANDIDATE",
    "Tier_3_Promising_Candidate": "PROMISING_CANDIDATE"
}

final_df["Final_Priority"] = (
    final_df["Integrated_Evidence_Tier"]
    .map(tier_map)
    .fillna("UNCLASSIFIED")
)

print("\nFinal priority counts:")
print(
    final_df["Final_Priority"]
    .value_counts()
)

print("\nTOP FINAL CANDIDATES")
print(
    final_df[
        [
            "Final_Rank",
            "Index_SNV",
            "Primary_gene",
            "Most_severe_consequence",
            "Protein_altering_probability",
            "Direct_T2D_Associations",
            "External_GWAS_Evidence_Score",
            "Integrated_Evidence_Score",
            "Integrated_Evidence_Tier",
            "Final_Priority"
        ]
    ]
    .to_string(index=False)
)

# Save shortlist
final_df.to_csv(
    "T2D_phase11D_final_candidate_shortlist.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 11D.2 — MODEL / BIOLOGY DISCORDANCE
# ============================================================

print("\n" + "=" * 80)
print("11D.2 — MODEL / BIOLOGY DISCORDANCE ANALYSIS")
print("=" * 80)

# False positives:
# model predicts protein-altering but biological label is 0
false_positive = final_df[
    (final_df["Has_protein_altering"] == 0) &
    (final_df["Predicted_label"] == 1)
].copy()

# False negatives:
# biological label is 1 but model predicts 0
false_negative = final_df[
    (final_df["Has_protein_altering"] == 1) &
    (final_df["Predicted_label"] == 0)
].copy()

# True positives
true_positive = final_df[
    (final_df["Has_protein_altering"] == 1) &
    (final_df["Predicted_label"] == 1)
].copy()

print("True positives:", len(true_positive))
print("False positives:", len(false_positive))
print("False negatives:", len(false_negative))

discordance = pd.concat(
    [
        false_positive.assign(
            Discordance_Type="False_Positive"
        ),
        false_negative.assign(
            Discordance_Type="False_Negative"
        )
    ],
    ignore_index=True
)

print("\nDISCORDANT CANDIDATES")

print(
    discordance[
        [
            "Index_SNV",
            "Primary_gene",
            "Most_severe_consequence",
            "Discordance_Type",
            "Protein_altering_probability",
            "Has_protein_altering",
            "Predicted_label",
            "Direct_T2D_Associations",
            "External_GWAS_Evidence_Score",
            "Integrated_Evidence_Score",
            "Integrated_Evidence_Tier"
        ]
    ]
    .sort_values(
        "Integrated_Evidence_Score",
        ascending=False
    )
    .to_string(index=False)
)

discordance.to_csv(
    "T2D_phase11D_model_biology_discordance.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 11D.3 — COMPLETE EVIDENCE MATRIX
# ============================================================

print("\n" + "=" * 80)
print("11D.3 — COMPLETE CANDIDATE EVIDENCE MATRIX")
print("=" * 80)

evidence_columns = [
    "Final_Rank",
    "Index_SNV",
    "Primary_gene",
    "Most_severe_consequence",
    "Validation_Group",
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
    "Final_Priority"
]

evidence_matrix = final_df[
    evidence_columns
].copy()

# ============================================================
# 11D.4 — BIOLOGICAL INTERPRETATION CLASSIFICATION
# ============================================================

print("\n" + "=" * 80)
print("11D.4 — BIOLOGICAL INTERPRETATION")
print("=" * 80)


def classify_candidate(row):

    protein = row["Has_protein_altering"]
    predicted = row["Predicted_label"]
    consequence = str(
        row["Most_severe_consequence"]
    )
    direct_t2d = row["Direct_T2D_Associations"]
    external_score = row[
        "External_GWAS_Evidence_Score"
    ]

    # Model-supported protein-altering candidates
    if protein == 1 and predicted == 1:

        if (
            "stop_gained" in consequence
            or "missense" in consequence
        ):
            return (
                "Protein_altering_model_supported"
            )

        return (
            "Protein_altering_model_supported"
        )

    # Biologically supported but missed by model
    if protein == 1 and predicted == 0:

        if direct_t2d > 0:
            return (
                "Model_missed_but_external_T2D_supported"
            )

        return (
            "Model_missed_protein_altering_candidate"
        )

    # Model-supported but not protein altering
    if protein == 0 and predicted == 1:

        if direct_t2d > 0:
            return (
                "Non_protein_altering_but_external_T2D_supported"
            )

        return (
            "Model_supported_non_protein_altering_candidate"
        )

    return "Other"


evidence_matrix["Biological_Interpretation"] = (
    evidence_matrix.apply(
        classify_candidate,
        axis=1
    )
)

print(
    evidence_matrix[
        [
            "Final_Rank",
            "Index_SNV",
            "Primary_gene",
            "Most_severe_consequence",
            "Biological_Interpretation",
            "Integrated_Evidence_Tier"
        ]
    ]
    .to_string(index=False)
)

# ============================================================
# 11D.5 — FINAL RANKING TABLE
# ============================================================

print("\n" + "=" * 80)
print("11D.5 — FINAL CANDIDATE RANKING")
print("=" * 80)

final_ranking = evidence_matrix[
    [
        "Final_Rank",
        "Index_SNV",
        "Primary_gene",
        "Most_severe_consequence",
        "Biological_Interpretation",
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
        "Final_Priority"
    ]
].copy()

print(
    final_ranking.to_string(index=False)
)

final_ranking.to_csv(
    "T2D_phase11D_final_candidate_ranking.tsv",
    sep="\t",
    index=False
)

# ============================================================
# SUMMARY STATISTICS
# ============================================================

print("\n" + "=" * 80)
print("PHASE 11D SUMMARY")
print("=" * 80)

print(
    "\nTotal candidates:",
    len(final_df)
)

print(
    "\nPriority distribution:"
)

print(
    final_df["Final_Priority"]
    .value_counts()
    .to_string()
)

print(
    "\nBiological interpretation distribution:"
)

print(
    evidence_matrix[
        "Biological_Interpretation"
    ]
    .value_counts()
    .to_string()
)

print(
    "\nHighest-ranked candidates:"
)

print(
    final_ranking[
        [
            "Final_Rank",
            "Index_SNV",
            "Primary_gene",
            "Integrated_Evidence_Score",
            "Final_Priority"
        ]
    ]
    .head(10)
    .to_string(index=False)
)

# ============================================================
# FINAL QC
# ============================================================

print("\n" + "=" * 80)
print("FINAL QC")
print("=" * 80)

checks = {
    "17 candidates retained":
        len(final_df) == 17,

    "17 unique SNPs":
        final_df["Index_SNV"].nunique() == 17,

    "No missing SNPs":
        final_df["Index_SNV"].notna().all(),

    "No duplicate SNPs":
        not final_df["Index_SNV"].duplicated().any(),

    "All candidates ranked":
        final_df["Final_Rank"].notna().all(),

    "All candidates prioritized":
        final_df["Final_Priority"].notna().all()
}

for name, result in checks.items():
    print(
        f"{name}: {'PASS' if result else 'FAIL'}"
    )

if not all(checks.values()):
    raise ValueError(
        "One or more final QC checks failed."
    )

# ============================================================
# OUTPUT FILES
# ============================================================

print("\n" + "=" * 80)
print("SAVED FILES")
print("=" * 80)

print(
    "T2D_phase11D_final_candidate_shortlist.tsv"
)

print(
    "T2D_phase11D_model_biology_discordance.tsv"
)

print(
    "T2D_phase11D_final_candidate_ranking.tsv"
)

print(
    "\n" + "=" * 80
)

print(
    "PHASE 11D COMPLETE"
)

print("=" * 80)
