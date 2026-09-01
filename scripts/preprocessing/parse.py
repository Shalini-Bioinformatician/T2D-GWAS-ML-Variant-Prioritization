
import pandas as pd

# Load cleaned ST4 dataset
file = "ST4_clean.tsv"
df = pd.read_csv(file, sep="\t")

print("Original shape:", df.shape)

# Extract OR and 95% confidence interval values
extracted = df["OR_95CI"].str.extract(
    r'^\s*([0-9.]+)\s*\(\s*([0-9.]+)\s*-\s*([0-9.]+)\s*\)\s*$'
)

# Rename extracted columns
extracted.columns = ["OR", "CI_lower", "CI_upper"]

# Convert to numeric
for col in extracted.columns:
    extracted[col] = pd.to_numeric(extracted[col], errors="coerce")

# Add new columns to the dataset
df = pd.concat([df, extracted], axis=1)

# Save feature-engineered dataset
output = "ST4_feature_engineered.tsv"
df.to_csv(output, sep="\t", index=False)

# --------------------------------------------------
# QC CHECKS
# --------------------------------------------------

print("\nFEATURE ENGINEERING OUTPUT")
print("=" * 50)
print("Output file:", output)
print("Rows:", len(df))
print("Columns:", len(df.columns))

print("\nMissing values after extraction:")
print(df[["OR", "CI_lower", "CI_upper"]].isna().sum())

print("\nPositive-value checks:")
print("OR <= 0:", (df["OR"] <= 0).sum())
print("CI_lower <= 0:", (df["CI_lower"] <= 0).sum())
print("CI_upper <= 0:", (df["CI_upper"] <= 0).sum())

print("\nConfidence interval ordering checks:")
print("CI_lower >= OR:", (df["CI_lower"] >= df["OR"]).sum())
print("OR >= CI_upper:", (df["OR"] >= df["CI_upper"]).sum())
print("Valid CI ordering (CI_lower < OR < CI_upper):",
      ((df["CI_lower"] < df["OR"]) &
       (df["OR"] < df["CI_upper"])).sum())

print("\nFirst 10 parsed records:")
print(
    df[
        ["Index_SNV", "OR_95CI", "OR", "CI_lower", "CI_upper"]
    ].head(10).to_string(index=False)
)

print("\nOR summary:")
print(df["OR"].describe())
