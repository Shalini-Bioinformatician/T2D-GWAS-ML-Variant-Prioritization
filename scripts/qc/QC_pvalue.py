
import pandas as pd

file = "ST4_clean.tsv"

df = pd.read_csv(file, sep="\t")

zero_p = df[df["MR_MEGA_Association_P"] == 0].copy()

print("=" * 70)
print("ST4 QC-3: ZERO P-VALUE INVESTIGATION")
print("=" * 70)

print("\nNumber of zero P-values:", len(zero_p))

print("\nAffected records:")
print(
    zero_p[
        [
            "Locus",
            "Chromosome",
            "Index_SNV",
            "Position_b37",
            "Risk_Allele",
            "Other_Allele",
            "MR_MEGA_Association_P",
            "OR_95CI",
            "Effective_Sample_Size"
        ]
    ].to_string(index=False)
)

print("\n" + "=" * 70)
