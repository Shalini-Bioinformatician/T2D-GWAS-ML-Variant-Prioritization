
import pandas as pd
import numpy as np

input_file = "ST4_feature_engineered.tsv"
output_file = "ST4_feature_engineered_v2.tsv"

# --------------------------------------------------
# Load existing feature-engineered dataset
# --------------------------------------------------

df = pd.read_csv(input_file, sep="\t")

print("=" * 70)
print("FEATURE ENGINEERING - STEP 2")
print("Effect-size transformations")
print("=" * 70)

print("\nOriginal shape:", df.shape)

# --------------------------------------------------
# 1. Log odds ratio
# --------------------------------------------------

df["log_OR"] = np.log(df["OR"])

# --------------------------------------------------
# 2. Confidence interval width
# --------------------------------------------------

df["CI_width"] = df["CI_upper"] - df["CI_lower"]

# --------------------------------------------------
# 3. Standard error of log(OR)
#
# SE = [ln(CI_upper) - ln(CI_lower)] / (2 × 1.96)
# --------------------------------------------------

df["SE_log_OR"] = (
    np.log(df["CI_upper"]) -
    np.log(df["CI_lower"])
) / (2 * 1.96)

# --------------------------------------------------
# Validation
# --------------------------------------------------

print("\nNEW FEATURES CREATED:")
print("log_OR")
print("CI_width")
print("SE_log_OR")

print("\nMissing values:")
print(
    df[
        ["log_OR", "CI_width", "SE_log_OR"]
    ].isna().sum()
)

print("\nInfinite values:")
print(
    np.isinf(
        df[
            ["log_OR", "CI_width", "SE_log_OR"]
        ]
    ).sum()
)

print("\nPositive-value checks:")

print(
    "CI_width <= 0:",
    (df["CI_width"] <= 0).sum()
)

print(
    "SE_log_OR <= 0:",
    (df["SE_log_OR"] <= 0).sum()
)

print("\nLog OR summary:")
print(df["log_OR"].describe())

print("\nCI width summary:")
print(df["CI_width"].describe())

print("\nSE(log OR) summary:")
print(df["SE_log_OR"].describe())

# --------------------------------------------------
# Check mathematical consistency
# --------------------------------------------------

expected_log_or = np.log(df["OR"])

difference = np.abs(
    df["log_OR"] - expected_log_or
)

print("\nMaximum log_OR calculation difference:")
print(difference.max())

# --------------------------------------------------
# Display examples
# --------------------------------------------------

print("\nFirst 10 transformed records:")

print(
    df[
        [
            "Index_SNV",
            "OR",
            "CI_lower",
            "CI_upper",
            "log_OR",
            "CI_width",
            "SE_log_OR"
        ]
    ]
    .head(10)
    .to_string(index=False)
)

# --------------------------------------------------
# Save
# --------------------------------------------------

df.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("\n" + "=" * 70)
print("FEATURE ENGINEERING STEP 2 COMPLETE")
print("=" * 70)

print("\nOutput file:", output_file)
print("Rows:", len(df))
print("Columns:", len(df.columns))

