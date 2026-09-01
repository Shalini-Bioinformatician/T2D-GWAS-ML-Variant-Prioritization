import pandas as pd

# ---------------------------------------------------------
# FILES
# ---------------------------------------------------------

original_file = "ST4_vep_GRCh37_annotations.tsv"
recovered_file_1 = "ST4_missing_vep_annotations.tsv"
recovered_file_2 = "ST4_remaining_vep_annotations.tsv"

output_file = "ST4_combined_vep_annotations.tsv"

# ---------------------------------------------------------
# READ FILES
# ---------------------------------------------------------

print("=" * 70)
print("COMBINING VEP GRCh37 ANNOTATIONS")
print("=" * 70)

files = [
    original_file,
    recovered_file_1,
    recovered_file_2
]

dfs = []

for file in files:

    print(f"\nReading: {file}")

    df = pd.read_csv(
        file,
        sep="\t",
        dtype=str
    )

    print("  Rows:", len(df))
    print("  Unique SNPs:", df["SNP"].nunique())

    dfs.append(df)

# ---------------------------------------------------------
# COMBINE
# ---------------------------------------------------------

combined = pd.concat(
    dfs,
    ignore_index=True
)

print("\n" + "=" * 70)
print("COMBINED DATASET")
print("=" * 70)

print("Total annotation rows:", len(combined))
print("Unique SNPs:", combined["SNP"].nunique())

# ---------------------------------------------------------
# CHECK DUPLICATE SNP-LEVEL COVERAGE
# ---------------------------------------------------------

snp_counts = combined.groupby("SNP").size()

print("\nAnnotation rows per SNP:")
print(snp_counts.describe())

# ---------------------------------------------------------
# COMPARE AGAINST ORIGINAL 1289 SNPs
# ---------------------------------------------------------

input_file = "ST4_annotation_input.tsv"

input_df = pd.read_csv(
    input_file,
    sep="\t",
    dtype=str
)

input_snps = set(input_df["SNP"])

annotated_snps = set(combined["SNP"])

missing = input_snps - annotated_snps
extra = annotated_snps - input_snps

print("\n" + "=" * 70)
print("COVERAGE CHECK")
print("=" * 70)

print("Original SNPs:", len(input_snps))
print("Annotated SNPs:", len(annotated_snps))
print("Missing SNPs:", len(missing))
print("Extra SNPs:", len(extra))

coverage = (
    len(annotated_snps) /
    len(input_snps)
) * 100

print(f"Annotation coverage: {coverage:.2f}%")

if missing:
    print("\nWARNING: Missing SNPs:")
    for snp in sorted(missing):
        print(snp)

if extra:
    print("\nWARNING: Extra SNPs:")
    for snp in sorted(extra):
        print(snp)

# ---------------------------------------------------------
# REMOVE EXACT DUPLICATE ANNOTATION ROWS
# ---------------------------------------------------------

before = len(combined)

combined = combined.drop_duplicates()

after = len(combined)

print("\nExact duplicate rows removed:", before - after)
print("Final annotation rows:", after)
print("Final unique SNPs:", combined["SNP"].nunique())

# ---------------------------------------------------------
# SAVE
# ---------------------------------------------------------

combined.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("\n" + "=" * 70)
print("COMBINATION COMPLETE")
print("=" * 70)

print("Output:")
print(output_file)
