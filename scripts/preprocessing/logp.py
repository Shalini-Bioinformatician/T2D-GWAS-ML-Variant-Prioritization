
import pandas as pd
import numpy as np

input_file = "ST4_feature_engineered_v2.tsv"
output_file = "ST4_feature_engineered_v3.tsv"

df = pd.read_csv(input_file, sep="\t")

print("=" * 70)
print("FEATURE ENGINEERING - STEP 3")
print("P-value transformation")
print("=" * 70)

print("\nOriginal shape:", df.shape)

# --------------------------------------------------
# 1. Create zero-P indicator
# --------------------------------------------------

df["P_is_zero"] = (
    df["MR_MEGA_Association_P"] == 0
).astype(int)

# --------------------------------------------------
# 2. Calculate -log10(P)
#
# Only calculate for positive P-values.
# Zero P-values remain NaN because their exact
# numerical values are unavailable in the source.
# --------------------------------------------------

df["neg_log10_P"] = np.where(
    df["MR_MEGA_Association_P"] > 0,
    -np.log10(df["MR_MEGA_Association_P"]),
    np.nan
)

# --------------------------------------------------
# 3. QC
# --------------------------------------------------

print("\nP-value records:")
print("Total:", len(df))
print(
    "P-values > 0:",
    (df["MR_MEGA_Association_P"] > 0).sum()
)
print(
    "P-values == 0:",
    (df["MR_MEGA_Association_P"] == 0).sum()
)

print("\nP_is_zero distribution:")
print(df["P_is_zero"].value_counts().sort_index())

print("\nMissing neg_log10_P:")
print(df["neg_log10_P"].isna().sum())

print("\nneg_log10_P summary:")
print(df["neg_log10_P"].describe())

# --------------------------------------------------
# 4. Check for invalid values
# --------------------------------------------------

print("\nInfinite neg_log10_P values:")
print(np.isinf(df["neg_log10_P"]).sum())

print("\nNegative neg_log10_P values:")
print(
    (df["neg_log10_P"].dropna() < 0).sum()
)

# --------------------------------------------------
# 5. Display most significant non-zero P-values
# --------------------------------------------------

print("\nTop 10 non-zero P-value signals:")

top = df[
    [
        "Index_SNV",
        "Locus",
        "MR_MEGA_Association_P",
        "neg_log10_P",
        "P_is_zero"
    ]
].sort_values(
    "neg_log10_P",
    ascending=False
)

print(top.head(10).to_string(index=False))

# --------------------------------------------------
# 6. Display the five zero-P variants
# --------------------------------------------------

print("\nSource-reported zero P-value variants:")

zero_rows = df[df["P_is_zero"] == 1]

print(
    zero_rows[
        [
            "Index_SNV",
            "Locus",
            "MR_MEGA_Association_P",
            "neg_log10_P",
            "P_is_zero"
        ]
    ].to_string(index=False)
)

# --------------------------------------------------
# 7. Save
# --------------------------------------------------

df.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("\n" + "=" * 70)
print("FEATURE ENGINEERING STEP 3 COMPLETE")
print("=" * 70)

print("\nOutput file:", output_file)
print("Rows:", len(df))
print("Columns:", len(df.columns))

