import pandas as pd
import re

print("=" * 70)
print("PHASE 11C.2 — GWAS CATALOG TRAIT CLASSIFICATION")
print("=" * 70)

# ---------------------------------------------------------
# LOAD GWAS CATALOG ASSOCIATIONS
# ---------------------------------------------------------

input_file = "phase11C_all_gwas_associations.tsv"

df = pd.read_csv(
    input_file,
    sep="\t",
    dtype=str
)

print("\nInput shape:", df.shape)
print("Unique SNPs:", df["Query_SNP"].nunique())

# ---------------------------------------------------------
# COMBINE TRAIT INFORMATION
# ---------------------------------------------------------

df["TRAIT_TEXT"] = (
    df["DISEASE/TRAIT"].fillna("") + " | " +
    df["MAPPED_TRAIT"].fillna("")
).str.lower()

# ---------------------------------------------------------
# TRAIT CLASSIFICATION RULES
# ---------------------------------------------------------

direct_t2d = [
    r"\btype 2 diabetes\b",
    r"\btype 2 diabetes mellitus\b",
    r"\bt2d\b",
    r"\bnon[- ]insulin[- ]dependent diabetes\b"
]

glycemic_insulin = [
    r"\bglucose\b",
    r"\bglycaemia\b",
    r"\bglycemia\b",
    r"\bfasting glucose\b",
    r"\bfasting insulin\b",
    r"\binsulin\b",
    r"\bhba1c\b",
    r"\bhaemoglobin a1c\b",
    r"\bhemoglobin a1c\b",
    r"\bproinsulin\b",
    r"\binsulin secretion\b",
    r"\binsulin sensitivity\b",
    r"\binsulin resistance\b"
]

metabolic_obesity = [
    r"\bobesity\b",
    r"\bbody mass index\b",
    r"\bbmi\b",
    r"\bwaist circumference\b",
    r"\bbody fat\b",
    r"\badiposity\b",
    r"\bmetabolic syndrome\b"
]

lipid_cardiovascular = [
    r"\bhdl\b",
    r"\bldl\b",
    r"\btriglyceride",
    r"\bcholesterol\b",
    r"\blipid\b",
    r"\bcoronary artery disease\b",
    r"\bcardiovascular\b"
]

# ---------------------------------------------------------
# CLASSIFICATION FUNCTION
# ---------------------------------------------------------

def contains_any(text, patterns):

    return any(
        re.search(pattern, text, flags=re.IGNORECASE)
        for pattern in patterns
    )


def classify_trait(text):

    # Direct T2D takes priority
    if contains_any(text, direct_t2d):
        return "Direct_T2D"

    if contains_any(text, glycemic_insulin):
        return "Glycemic_Insulin"

    if contains_any(text, metabolic_obesity):
        return "Metabolic_Obesity"

    if contains_any(text, lipid_cardiovascular):
        return "Lipid_Cardiovascular"

    return "Other"


df["Trait_Category"] = df["TRAIT_TEXT"].apply(
    classify_trait
)

# ---------------------------------------------------------
# BINARY EVIDENCE FLAGS
# ---------------------------------------------------------

df["Is_Direct_T2D"] = (
    df["Trait_Category"] == "Direct_T2D"
).astype(int)

df["Is_Glycemic"] = (
    df["Trait_Category"] == "Glycemic_Insulin"
).astype(int)

df["Is_Metabolic"] = (
    df["Trait_Category"] == "Metabolic_Obesity"
).astype(int)

df["Is_Lipid_Cardiovascular"] = (
    df["Trait_Category"] == "Lipid_Cardiovascular"
).astype(int)

# ---------------------------------------------------------
# TRAIT CATEGORY QC
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("TRAIT CATEGORY COUNTS")
print("=" * 70)

print(
    df["Trait_Category"]
    .value_counts()
    .to_string()
)

print("\n" + "=" * 70)
print("TRAIT CATEGORY PERCENTAGES")
print("=" * 70)

print(
    (df["Trait_Category"].value_counts(normalize=True) * 100)
    .round(2)
    .to_string()
)

# ---------------------------------------------------------
# SNP-LEVEL EXTERNAL EVIDENCE SUMMARY
# ---------------------------------------------------------

summary = (
    df.groupby("Query_SNP")
    .agg(
        Total_GWAS_Associations=("Query_SNP", "size"),

        Direct_T2D_Associations=(
            "Is_Direct_T2D",
            "sum"
        ),

        Glycemic_Associations=(
            "Is_Glycemic",
            "sum"
        ),

        Metabolic_Associations=(
            "Is_Metabolic",
            "sum"
        ),

        Lipid_Cardiovascular_Associations=(
            "Is_Lipid_Cardiovascular",
            "sum"
        )
    )
    .reset_index()
)

# ---------------------------------------------------------
# SNP-LEVEL EVIDENCE FLAGS
# ---------------------------------------------------------

summary["Has_Direct_T2D_Evidence"] = (
    summary["Direct_T2D_Associations"] > 0
).astype(int)

summary["Has_Glycemic_Evidence"] = (
    summary["Glycemic_Associations"] > 0
).astype(int)

summary["Has_Metabolic_Evidence"] = (
    summary["Metabolic_Associations"] > 0
).astype(int)

summary["Has_Lipid_Cardiovascular_Evidence"] = (
    summary["Lipid_Cardiovascular_Associations"] > 0
).astype(int)

# ---------------------------------------------------------
# OVERALL BIOLOGICAL EVIDENCE FLAG
# ---------------------------------------------------------

summary["Has_T2D_Related_Evidence"] = (
    (
        summary["Has_Direct_T2D_Evidence"] == 1
    ) |
    (
        summary["Has_Glycemic_Evidence"] == 1
    ) |
    (
        summary["Has_Metabolic_Evidence"] == 1
    ) |
    (
        summary["Has_Lipid_Cardiovascular_Evidence"] == 1
    )
).astype(int)

# ---------------------------------------------------------
# SORT
# ---------------------------------------------------------

summary = summary.sort_values(
    [
        "Has_Direct_T2D_Evidence",
        "Has_T2D_Related_Evidence",
        "Total_GWAS_Associations"
    ],
    ascending=[False, False, False]
)

# ---------------------------------------------------------
# SAVE OUTPUTS
# ---------------------------------------------------------

classified_file = (
    "phase11C_trait_classified_associations.tsv"
)

summary_file = (
    "phase11C_external_validation_summary.tsv"
)

df.to_csv(
    classified_file,
    sep="\t",
    index=False
)

summary.to_csv(
    summary_file,
    sep="\t",
    index=False
)

# ---------------------------------------------------------
# DISPLAY SNP SUMMARY
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("SNP-LEVEL EXTERNAL VALIDATION SUMMARY")
print("=" * 70)

print(
    summary.to_string(index=False)
)

# ---------------------------------------------------------
# FINAL QC
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("FINAL QC")
print("=" * 70)

print("Input association rows:", len(df))
print("Output association rows:", len(df))
print("Unique input SNPs:", df["Query_SNP"].nunique())
print("Unique summary SNPs:", summary["Query_SNP"].nunique())

print("\nSaved:")
print(classified_file)
print(summary_file)

print("\n" + "=" * 70)
print("PHASE 11C.2 COMPLETE")
print("=" * 70)
