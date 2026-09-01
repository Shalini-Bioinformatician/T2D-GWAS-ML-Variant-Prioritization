import pandas as pd

print("=" * 70)
print("PHASE 11C.3 — EXTERNAL GWAS EVIDENCE SCORING")
print("=" * 70)

# ---------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------

input_file = "phase11C_external_validation_summary.tsv"

df = pd.read_csv(
    input_file,
    sep="\t"
)

print("\nInput shape:", df.shape)
print("Unique SNPs:", df["Query_SNP"].nunique())


# ---------------------------------------------------------
# EXTERNAL EVIDENCE SCORE
# ---------------------------------------------------------
#
# Direct T2D evidence receives the highest weight.
# Related metabolic traits provide supporting evidence.
#
# We use log1p counts so that SNPs with hundreds/thousands
# of associations do not dominate the ranking.
# ---------------------------------------------------------

import numpy as np

df["Direct_T2D_score"] = np.log1p(
    df["Direct_T2D_Associations"]
)

df["Glycemic_score"] = np.log1p(
    df["Glycemic_Associations"]
)

df["Metabolic_score"] = np.log1p(
    df["Metabolic_Associations"]
)

df["Lipid_Cardiovascular_score"] = np.log1p(
    df["Lipid_Cardiovascular_Associations"]
)


# ---------------------------------------------------------
# WEIGHTED EXTERNAL GWAS SCORE
# ---------------------------------------------------------

df["External_GWAS_Evidence_Score"] = (
    5.0 * df["Direct_T2D_score"] +
    3.0 * df["Glycemic_score"] +
    2.0 * df["Metabolic_score"] +
    1.0 * df["Lipid_Cardiovascular_score"]
)


# ---------------------------------------------------------
# EVIDENCE TIER
# ---------------------------------------------------------

def assign_tier(row):

    direct = row["Direct_T2D_Associations"]
    glycemic = row["Glycemic_Associations"]
    metabolic = row["Metabolic_Associations"]

    if direct >= 3 and (glycemic >= 1 or metabolic >= 1):
        return "Tier_1_Strong_T2D_related"

    elif direct >= 2:
        return "Tier_2_Replicated_T2D"

    elif direct >= 1 and (glycemic >= 1 or metabolic >= 1):
        return "Tier_3_T2D_plus_related"

    elif direct >= 1:
        return "Tier_4_Direct_T2D_only"

    else:
        return "Tier_5_No_direct_T2D"


df["External_Evidence_Tier"] = df.apply(
    assign_tier,
    axis=1
)


# ---------------------------------------------------------
# SORT
# ---------------------------------------------------------

df = df.sort_values(
    [
        "External_Evidence_Tier",
        "External_GWAS_Evidence_Score",
        "Direct_T2D_Associations"
    ],
    ascending=[True, False, False]
).reset_index(drop=True)


# ---------------------------------------------------------
# RANK
# ---------------------------------------------------------

df["External_GWAS_Rank"] = (
    df["External_GWAS_Evidence_Score"]
    .rank(
        ascending=False,
        method="min"
    )
    .astype(int)
)


# ---------------------------------------------------------
# DISPLAY
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("EXTERNAL GWAS EVIDENCE RANKING")
print("=" * 70)

display_columns = [
    "External_GWAS_Rank",
    "Query_SNP",
    "Total_GWAS_Associations",
    "Direct_T2D_Associations",
    "Glycemic_Associations",
    "Metabolic_Associations",
    "Lipid_Cardiovascular_Associations",
    "Has_T2D_Related_Evidence",
    "External_GWAS_Evidence_Score",
    "External_Evidence_Tier"
]

print(
    df[display_columns]
    .sort_values("External_GWAS_Rank")
    .to_string(index=False)
)


# ---------------------------------------------------------
# TIER COUNTS
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("EVIDENCE TIER COUNTS")
print("=" * 70)

print(
    df["External_Evidence_Tier"]
    .value_counts()
    .sort_index()
)


# ---------------------------------------------------------
# SAVE
# ---------------------------------------------------------

output_file = (
    "phase11C_external_GWAS_evidence_ranking.tsv"
)

df.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("\nSaved:")
print(output_file)

print("\n" + "=" * 70)
print("PHASE 11C.3 COMPLETE")
print("=" * 70)
