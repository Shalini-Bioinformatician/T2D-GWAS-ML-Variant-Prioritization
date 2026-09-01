
import pandas as pd

input_file = "ST4_feature_engineered_v4.tsv"
output_file = "ST4_annotation_input.tsv"

df = pd.read_csv(input_file, sep="\t")

annotation = df[
    [
        "Index_SNV",
        "Chromosome",
        "Position_b37",
        "Risk_Allele",
        "Other_Allele"
    ]
].copy()

annotation.columns = [
    "SNP",
    "Chromosome",
    "Position_b37",
    "Risk_Allele",
    "Other_Allele"
]

# Check that we have the expected number of variants
print("Original feature-engineered rows:", len(df))
print("Annotation input rows:", len(annotation))
print("Unique SNPs:", annotation["SNP"].nunique())

# Check for missing values
print("\nMissing values:")
print(annotation.isnull().sum())

# Check chromosome range
print("\nChromosomes:")
print(sorted(annotation["Chromosome"].unique()))

# Save annotation input
annotation.to_csv(output_file, sep="\t", index=False)

print("\nAnnotation input created:")
print(output_file)

print("\nFirst 10 records:")
print(annotation.head(10).to_string(index=False))
