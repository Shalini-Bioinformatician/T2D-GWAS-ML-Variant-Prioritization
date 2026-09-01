
import pandas as pd
import numpy as np

input_file = "ST4_feature_engineered_v3.tsv"
output_file = "ST4_feature_engineered_v4.tsv"

df = pd.read_csv(input_file, sep="\t")

print("=" * 70)
print("FEATURE ENGINEERING - STEP 4")
print("Minor Allele Frequency (MAF)")
print("=" * 70)

print("\nOriginal shape:", df.shape)

# --------------------------------------------------
# 1. Calculate MAF
# --------------------------------------------------

df["MAF"] = np.minimum(
    df["Risk_Allele_Frequency"],
    1 - df["Risk_Allele_Frequency"]
)

# --------------------------------------------------
# 2. Basic validation
# --------------------------------------------------

print("\nMAF missing values:")
print(df["MAF"].isna().sum())

print("\nMAF infinite values:")
print(np.isinf(df["MAF"]).sum())

print("\nMAF outside valid range [0, 0.5]:")
invalid_maf = (
    (df["MAF"] < 0) |
    (df["MAF"] > 0.5)
).sum()

print(invalid_maf)

# --------------------------------------------------
# 3. Check relationship between RAF and MAF
# --------------------------------------------------

expected_maf = np.minimum(
    df["Risk_Allele_Frequency"],
    1 - df["Risk_Allele_Frequency"]
)

difference = np.abs(
    df["MAF"] - expected_maf
)

print("\nMaximum RAF → MAF calculation difference:")
print(difference.max())

# --------------------------------------------------
# 4. MAF summary
# --------------------------------------------------

print("\nMAF summary:")
print(df["MAF"].describe())

# --------------------------------------------------
# 5. Frequency categories
#
# These are descriptive only.
# We are NOT filtering variants based on them.
# --------------------------------------------------

def classify_maf(maf):
    if maf < 0.01:
        return "Rare"
    elif maf < 0.05:
        return "Low_frequency"
    else:
        return "Common"

df["MAF_category"] = df["MAF"].apply(classify_maf)

print("\nMAF category distribution:")
print(
    df["MAF_category"]
    .value_counts()
    .reindex(
        ["Rare", "Low_frequency", "Common"],
        fill_value=0
    )
)

# --------------------------------------------------
# 6. Display examples
# --------------------------------------------------

print("\nFirst 10 frequency records:")

print(
    df[
        [
            "Index_SNV",
            "Risk_Allele",
            "Other_Allele",
            "Risk_Allele_Frequency",
            "MAF",
            "MAF_category"
        ]
    ]
    .head(10)
    .to_string(index=False)
)

# --------------------------------------------------
# 7. Display rare variants
# --------------------------------------------------

print("\nLowest-frequency variants:")

print(
    df[
        [
            "Index_SNV",
            "Locus",
            "Risk_Allele_Frequency",
            "MAF",
            "MAF_category"
        ]
    ]
    .sort_values("MAF")
    .head(10)
    .to_string(index=False)
)

# --------------------------------------------------
# 8. Save
# --------------------------------------------------

df.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("\n" + "=" * 70)
print("FEATURE ENGINEERING STEP 4 COMPLETE")
print("=" * 70)

print("\nOutput file:", output_file)
print("Rows:", len(df))
print("Columns:", len(df.columns))
