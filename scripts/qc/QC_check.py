import pandas as pd

file = "ST4_clean.tsv"

df = pd.read_csv(file, sep="\t")

print("=" * 70)
print("ST4 BASIC QC")
print("=" * 70)

# 1. Dataset dimensions
print("\n1. DATASET SIZE")
print(f"Rows:    {len(df):,}")
print(f"Columns: {len(df.columns)}")

# 2. Column names
print("\n2. COLUMNS")
for i, col in enumerate(df.columns, 1):
    print(f"{i:2}. {col}")

# 3. Data types
print("\n3. DATA TYPES")
print(df.dtypes)

# 4. Missing values
print("\n4. MISSING VALUES")
missing = df.isna().sum()
missing_pct = (missing / len(df) * 100).round(2)

missing_table = pd.DataFrame({
    "Missing_Count": missing,
    "Missing_Percent": missing_pct
})

print(missing_table)

# 5. Duplicate rows
print("\n5. DUPLICATE ROWS")
print("Duplicate complete rows:", df.duplicated().sum())

# 6. Duplicate Index SNVs
print("\n6. DUPLICATE INDEX SNVs")
print("Duplicate Index_SNVs:", df["Index_SNV"].duplicated().sum())

# 7. Duplicate genomic positions
print("\n7. DUPLICATE GENOMIC POSITIONS")
dup_pos = df.duplicated(
    subset=["Chromosome", "Position_b37"]
).sum()

print("Duplicate chromosome-position combinations:", dup_pos)

# 8. Chromosomes
print("\n8. CHROMOSOME DISTRIBUTION")
print(df["Chromosome"].value_counts().sort_index())

# 9. Invalid chromosome values
print("\n9. INVALID CHROMOSOMES")
valid_chr = set(range(1, 23))

invalid_chr = df.loc[
    ~df["Chromosome"].isin(valid_chr),
    "Chromosome"
].unique()

print("Invalid chromosomes:", invalid_chr)

# 10. Position checks
print("\n10. POSITION CHECK")
print("Minimum position:", df["Position_b37"].min())
print("Maximum position:", df["Position_b37"].max())
print(
    "Non-positive positions:",
    (df["Position_b37"] <= 0).sum()
)

# 11. Alleles
print("\n11. ALLELE CHECK")

valid_bases = {"A", "C", "G", "T"}

risk_invalid = ~df["Risk_Allele"].isin(valid_bases)
other_invalid = ~df["Other_Allele"].isin(valid_bases)

print("Invalid Risk_Alleles:", risk_invalid.sum())
print("Invalid Other_Alleles:", other_invalid.sum())

print("\nRisk allele distribution:")
print(df["Risk_Allele"].value_counts())

print("\nOther allele distribution:")
print(df["Other_Allele"].value_counts())

# 12. Same allele on both sides
print("\n12. SAME RISK/OTHER ALLELE")
same_allele = (
    df["Risk_Allele"] == df["Other_Allele"]
).sum()

print("Risk allele == Other allele:", same_allele)

# 13. RAF
print("\n13. RISK ALLELE FREQUENCY")
print(df["Risk_Allele_Frequency"].describe())

invalid_raf = (
    (df["Risk_Allele_Frequency"] < 0) |
    (df["Risk_Allele_Frequency"] > 1)
).sum()

print("RAF outside [0,1]:", invalid_raf)

print("\n" + "=" * 70)
print("BASIC QC COMPLETE")
print("=" * 70)
