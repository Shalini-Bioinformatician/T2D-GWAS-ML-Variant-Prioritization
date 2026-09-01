import pandas as pd
import numpy as np

print("=" * 70)
print("PHASE 11B.3 — COMPARATIVE BIOLOGICAL CANDIDATE ANALYSIS")
print("=" * 70)

# =========================================================
# FILES
# =========================================================

input_files = {
    "True_Positive":
        "T2D_phase11B_true_positive_validation.tsv",

    "High_Confidence_False_Positive":
        "T2D_phase11B_false_positive_validation.tsv",

    "False_Negative":
        "T2D_phase11B_false_negative_validation.tsv"
}

output_combined = (
    "T2D_phase11B_comparative_candidate_analysis.tsv"
)

output_summary = (
    "T2D_phase11B_group_summary.tsv"
)

output_consequence = (
    "T2D_phase11B_consequence_comparison.tsv"
)

# =========================================================
# LOAD VALIDATION GROUPS
# =========================================================

groups = []

for category, filename in input_files.items():

    print("\nLoading:", filename)

    temp = pd.read_csv(
        filename,
        sep="\t"
    )

    print("Rows:", len(temp))

    # Force the validation category
    temp["Validation_Group"] = category

    groups.append(temp)


# =========================================================
# COMBINE
# =========================================================

combined = pd.concat(
    groups,
    ignore_index=True
)

print("\n" + "=" * 70)
print("COMBINED VALIDATION PANEL")
print("=" * 70)

print("Total candidates:", len(combined))
print(
    "Unique SNPs:",
    combined["Index_SNV"].nunique()
)

print(
    "\nValidation groups:"
)

print(
    combined["Validation_Group"]
    .value_counts()
)


# =========================================================
# BASIC QC
# =========================================================

print("\n" + "=" * 70)
print("BASIC QC")
print("=" * 70)

if combined["Index_SNV"].duplicated().any():

    duplicated = combined[
        combined["Index_SNV"].duplicated(
            keep=False
        )
    ]

    print("\nWARNING: Duplicate SNPs detected:")
    print(
        duplicated[
            ["Index_SNV", "Validation_Group"]
        ].to_string(index=False)
    )

else:

    print("No duplicate SNPs.")

print(
    "\nMissing values:"
)

print(
    combined[
        [
            "Index_SNV",
            "Primary_gene",
            "Most_severe_consequence",
            "Protein_altering_probability",
            "OR",
            "neg_log10_P",
            "MAF"
        ]
    ]
    .isna()
    .sum()
)


# =========================================================
# GROUP-LEVEL NUMERICAL SUMMARY
# =========================================================

print("\n" + "=" * 70)
print("GROUP-LEVEL NUMERICAL COMPARISON")
print("=" * 70)

numeric_features = [
    "Protein_altering_probability",
    "OR",
    "neg_log10_P",
    "MAF"
]

summary_rows = []

for group, subset in combined.groupby(
    "Validation_Group"
):

    row = {
        "Validation_Group": group,
        "N": len(subset)
    }

    for feature in numeric_features:

        row[f"{feature}_mean"] = (
            subset[feature].mean()
        )

        row[f"{feature}_median"] = (
            subset[feature].median()
        )

        row[f"{feature}_min"] = (
            subset[feature].min()
        )

        row[f"{feature}_max"] = (
            subset[feature].max()
        )

    summary_rows.append(row)


group_summary = pd.DataFrame(
    summary_rows
)

print(
    group_summary.to_string(
        index=False
    )
)


# =========================================================
# FUNCTIONAL CONSEQUENCE COMPARISON
# =========================================================

print("\n" + "=" * 70)
print("FUNCTIONAL CONSEQUENCE COMPARISON")
print("=" * 70)

consequence_comparison = (
    pd.crosstab(
        combined["Most_severe_consequence"],
        combined["Validation_Group"]
    )
)

print(
    consequence_comparison.to_string()
)


# =========================================================
# CONSEQUENCE PERCENTAGES
# =========================================================

print("\n" + "=" * 70)
print("CONSEQUENCE DISTRIBUTION (%)")
print("=" * 70)

consequence_percent = (
    pd.crosstab(
        combined["Most_severe_consequence"],
        combined["Validation_Group"],
        normalize="columns"
    )
    * 100
)

print(
    consequence_percent.round(2).to_string()
)


# =========================================================
# MODEL PROBABILITY BANDS
# =========================================================

print("\n" + "=" * 70)
print("MODEL PROBABILITY BANDS")
print("=" * 70)

def probability_band(p):

    if p >= 0.90:
        return "Very_high_>=0.90"

    elif p >= 0.75:
        return "High_0.75-0.899"

    elif p >= 0.50:
        return "Moderate_0.50-0.749"

    else:
        return "Low_<0.50"


combined["Probability_band"] = (
    combined[
        "Protein_altering_probability"
    ].apply(probability_band)
)

probability_comparison = (
    pd.crosstab(
        combined["Probability_band"],
        combined["Validation_Group"]
    )
)

print(
    probability_comparison.to_string()
)


# =========================================================
# GWAS SIGNIFICANCE BANDS
# =========================================================

print("\n" + "=" * 70)
print("GWAS EVIDENCE BANDS")
print("=" * 70)

def pvalue_band(x):

    if x >= 50:
        return "Extreme_>=50"

    elif x >= 20:
        return "Very_strong_20-49.99"

    elif x >= 10:
        return "Strong_10-19.99"

    elif x >= 7.3:
        return "Genome_wide_7.3-9.99"

    else:
        return "Below_genome_wide_<7.3"


combined["GWAS_evidence_band"] = (
    combined["neg_log10_P"].apply(
        pvalue_band
    )
)

gwas_comparison = (
    pd.crosstab(
        combined["GWAS_evidence_band"],
        combined["Validation_Group"]
    )
)

print(
    gwas_comparison.to_string()
)


# =========================================================
# PROTEIN-ALTERING STATUS VS MODEL
# =========================================================

print("\n" + "=" * 70)
print("MODEL / BIOLOGICAL LABEL COMPARISON")
print("=" * 70)

comparison = (
    pd.crosstab(
        [
            combined["Has_protein_altering"],
            combined["Predicted_label"]
        ],
        combined["Validation_Group"]
    )
)

print(
    comparison.to_string()
)


# =========================================================
# INDIVIDUAL CANDIDATE RANKING
# =========================================================

print("\n" + "=" * 70)
print("ALL 17 CANDIDATES — RANKED BY MODEL PROBABILITY")
print("=" * 70)

ranking_columns = [
    "Index_SNV",
    "Primary_gene",
    "Most_severe_consequence",
    "Validation_Group",
    "Protein_altering_probability",
    "Has_protein_altering",
    "Predicted_label",
    "OR",
    "neg_log10_P",
    "MAF"
]

ranked = (
    combined[
        ranking_columns
    ]
    .sort_values(
        "Protein_altering_probability",
        ascending=False
    )
    .reset_index(drop=True)
)

print(
    ranked.to_string(index=False)
)


# =========================================================
# IDENTIFY MODEL DISCORDANCE
# =========================================================

print("\n" + "=" * 70)
print("MODEL-BIOLOGY DISCORDANCE")
print("=" * 70)

discordant = combined[
    combined["Has_protein_altering"]
    != combined["Predicted_label"]
].copy()

print(
    "Discordant candidates:",
    len(discordant)
)

print(
    discordant[
        ranking_columns
    ]
    .sort_values(
        "Protein_altering_probability",
        ascending=False
    )
    .to_string(index=False)
)


# =========================================================
# TRUE POSITIVE BIOLOGICAL PROFILE
# =========================================================

print("\n" + "=" * 70)
print("TRUE-POSITIVE PROFILE")
print("=" * 70)

tp = combined[
    combined["Validation_Group"]
    == "True_Positive"
]

print(
    "Mean ML probability:",
    round(
        tp[
            "Protein_altering_probability"
        ].mean(),
        4
    )
)

print(
    "Median ML probability:",
    round(
        tp[
            "Protein_altering_probability"
        ].median(),
        4
    )
)

print(
    "Mean -log10(P):",
    round(
        tp["neg_log10_P"].mean(),
        4
    )
)

print(
    "Median -log10(P):",
    round(
        tp["neg_log10_P"].median(),
        4
    )
)

print(
    "Missense variants:",
    (
        tp["Most_severe_consequence"]
        == "missense_variant"
    ).sum()
)

print(
    "Stop-gained variants:",
    (
        tp["Most_severe_consequence"]
        == "stop_gained"
    ).sum()
)


# =========================================================
# FALSE-POSITIVE PROFILE
# =========================================================

print("\n" + "=" * 70)
print("FALSE-POSITIVE PROFILE")
print("=" * 70)

fp = combined[
    combined["Validation_Group"]
    == "High_Confidence_False_Positive"
]

print(
    "Mean ML probability:",
    round(
        fp[
            "Protein_altering_probability"
        ].mean(),
        4
    )
)

print(
    "Median ML probability:",
    round(
        fp[
            "Protein_altering_probability"
        ].median(),
        4
    )
)

print(
    "Mean -log10(P):",
    round(
        fp["neg_log10_P"].mean(),
        4
    )
)

print(
    "Median -log10(P):",
    round(
        fp["neg_log10_P"].median(),
        4
    )
)


# =========================================================
# FALSE-NEGATIVE PROFILE
# =========================================================

print("\n" + "=" * 70)
print("FALSE-NEGATIVE PROFILE")
print("=" * 70)

fn = combined[
    combined["Validation_Group"]
    == "False_Negative"
]

print(
    "Mean ML probability:",
    round(
        fn[
            "Protein_altering_probability"
        ].mean(),
        4
    )
)

print(
    "Median ML probability:",
    round(
        fn[
            "Protein_altering_probability"
        ].median(),
        4
    )
)

print(
    "Mean -log10(P):",
    round(
        fn["neg_log10_P"].mean(),
        4
    )
)

print(
    "Median -log10(P):",
    round(
        fn["neg_log10_P"].median(),
        4
    )
)


# =========================================================
# SAVE FILES
# =========================================================

combined.to_csv(
    output_combined,
    sep="\t",
    index=False
)

group_summary.to_csv(
    output_summary,
    sep="\t",
    index=False
)

consequence_comparison.to_csv(
    output_consequence,
    sep="\t"
)


# =========================================================
# FINAL QC
# =========================================================

print("\n" + "=" * 70)
print("SAVED FILES")
print("=" * 70)

print(output_combined)
print(output_summary)
print(output_consequence)

print("\n" + "=" * 70)
print("PHASE 11B.3 COMPLETE")
print("=" * 70)
