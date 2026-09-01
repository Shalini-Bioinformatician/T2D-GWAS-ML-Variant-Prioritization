import pandas as pd

file = "Suzuki_et_al_2024_Supplementary_tables.xlsx"

# Read ST4 without assigning a header
raw = pd.read_excel(file, sheet_name="ST4", header=None)

# Data starts after the title and two header rows
data = raw.iloc[4:].copy()

# Assign clean column names
data.columns = [
    "Locus",
    "Chromosome",
    "Interval_b37",
    "Previously_reported",
    "Index_SNV",
    "Position_b37",
    "Risk_Allele",
    "Other_Allele",
    "Risk_Allele_Frequency",
    "MR_MEGA_Association_P",
    "OR_95CI",
    "Effective_Sample_Size"
]

# Keep rows that contain an actual index SNV
data = data.dropna(subset=["Index_SNV"]).copy()

# Forward-fill columns affected by merged cells in Excel
fill_columns = [
    "Locus",
    "Chromosome",
    "Interval_b37",
    "Previously_reported"
]

data[fill_columns] = data[fill_columns].ffill()

# Convert chromosome and position to numeric
data["Chromosome"] = pd.to_numeric(
    data["Chromosome"], errors="coerce"
)

data["Position_b37"] = pd.to_numeric(
    data["Position_b37"], errors="coerce"
)

# Remove genuinely invalid rows
data = data.dropna(
    subset=["Chromosome", "Position_b37", "Index_SNV"]
).copy()

# Convert to integers
data["Chromosome"] = data["Chromosome"].astype(int)
data["Position_b37"] = data["Position_b37"].astype(int)

# Remove duplicate index SNVs, if any
data = data.drop_duplicates(
    subset=["Index_SNV", "Chromosome", "Position_b37"]
).copy()

# Save corrected dataset
output = "ST4_clean.tsv"

data.to_csv(output, sep="\t", index=False)

print("Corrected clean file created:", output)
print("Number of index signals:", len(data))
print("Number of unique loci:", data["Locus"].nunique())
print("Number of columns:", len(data.columns))

print("\nFirst 10 rows:")
print(data.head(10).to_string(index=False))

print("\nExample rows from a locus with multiple signals:")
multi = data.groupby("Locus").filter(lambda x: len(x) > 1)
print(multi.head(10).to_string(index=False))

print("\nLast 5 rows:")
print(data.tail(5).to_string(index=False))
