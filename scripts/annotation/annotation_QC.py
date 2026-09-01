
import pandas as pd

file = "ST4_vep_GRCh37_annotations.tsv"

df = pd.read_csv(file, sep="\t")

print("=" * 70)
print("ST4 VEP GRCh37 ANNOTATION QC")
print("=" * 70)

print("\n1. DATASET SIZE")
print("Annotation rows:", len(df))
print("Columns:", len(df.columns))
print("Unique SNPs:", df["SNP"].nunique())

print("\n2. COLUMNS")
for i, col in enumerate(df.columns, 1):
    print(f"{i}. {col}")

print("\n3. MISSING VALUES")
missing = pd.DataFrame({
    "Missing_Count": df.isna().sum(),
    "Missing_Percent": df.isna().mean() * 100
})

print(missing.to_string())

print("\n4. ASSEMBLY")
print(df["Assembly"].value_counts(dropna=False))

print("\n5. ANNOTATION CHROMOSOME")
print(df["Annotation_Chromosome"].value_counts(dropna=False).sort_index())

print("\n6. MOST SEVERE CONSEQUENCE")
print(
    df[
        ["SNP", "Most_severe_consequence"]
    ]
    .drop_duplicates()
    ["Most_severe_consequence"]
    .value_counts(dropna=False)
)

print("\n7. CONSEQUENCE")
print(
    df["Consequence"]
    .value_counts(dropna=False)
    .head(30)
)

print("\n8. IMPACT")
print(
    df["Impact"]
    .value_counts(dropna=False)
)

print("\n9. GENE ANNOTATION")
unique_snps = df["SNP"].nunique()

gene_annotated = (
    df.loc[
        df["Gene_symbol"].notna(),
        "SNP"
    ].nunique()
)

print("Unique SNPs:", unique_snps)
print("SNPs with gene annotation:", gene_annotated)
print(
    "Gene annotation rate:",
    round(gene_annotated / unique_snps * 100, 2),
    "%"
)

print("\n10. INTERGENIC VARIANTS")

intergenic = (
    df.loc[
        df["Most_severe_consequence"] == "intergenic_variant",
        "SNP"
    ].nunique()
)

print("Intergenic SNPs:", intergenic)
print(
    "Intergenic percentage:",
    round(intergenic / unique_snps * 100, 2),
    "%"
)

print("\n11. TRANSCRIPT MULTIPLICITY")

transcript_counts = (
    df.groupby("SNP")
      .size()
)

print(transcript_counts.describe())

print(
    "\nSNPs with >1 transcript annotation:",
    (transcript_counts > 1).sum()
)

print(
    "Maximum transcript annotations for one SNP:",
    transcript_counts.max()
)

print("\n12. EXAMPLE ANNOTATIONS")

example_snps = [
    "rs12138597",
    "rs7766070",
    "rs10811661",
    "rs7903146",
    "rs1558902"
]

print(
    df[
        df["SNP"].isin(example_snps)
    ][
        [
            "SNP",
            "Most_severe_consequence",
            "Gene_symbol",
            "Consequence",
            "Impact"
        ]
    ].to_string(index=False)
)

print("\n" + "=" * 70)
print("ANNOTATION QC COMPLETE")
print("=" * 70)
