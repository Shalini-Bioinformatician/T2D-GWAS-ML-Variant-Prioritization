import pandas as pd

print("=" * 70)
print("PHASE 11B.2C — FALSE-NEGATIVE BIOLOGICAL VALIDATION SET")
print("=" * 70)

input_file = "T2D_phase11B_candidate_panel.tsv"
output_file = "T2D_phase11B_false_negative_validation.tsv"

df = pd.read_csv(
    input_file,
    sep="\t"
)

print("\nCandidate panel shape:", df.shape)

# ---------------------------------------------------------
# SELECT FALSE NEGATIVES
# ---------------------------------------------------------

false_negative_validation = df[
    df["Candidate_Category"] == "False_Negative"
].copy()

# ---------------------------------------------------------
# SORT BY MODEL PROBABILITY
# ---------------------------------------------------------

false_negative_validation = (
    false_negative_validation
    .sort_values(
        "Protein_altering_probability",
        ascending=False
    )
    .reset_index(drop=True)
)

# ---------------------------------------------------------
# DISPLAY
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("FALSE-NEGATIVE CANDIDATES")
print("=" * 70)

print(
    false_negative_validation[
        [
            "Index_SNV",
            "Primary_gene",
            "Most_severe_consequence",
            "Protein_altering_probability",
            "Has_protein_altering",
            "Predicted_label",
            "OR",
            "neg_log10_P",
            "MAF",
            "Candidate_Category"
        ]
    ].to_string(index=False)
)

# ---------------------------------------------------------
# QC
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("QC")
print("=" * 70)

print(
    "Number of candidates:",
    len(false_negative_validation)
)

print(
    "Unique SNPs:",
    false_negative_validation["Index_SNV"].nunique()
)

print(
    "Actual protein-altering labels:",
    false_negative_validation[
        "Has_protein_altering"
    ].value_counts().to_dict()
)

print(
    "Predicted labels:",
    false_negative_validation[
        "Predicted_label"
    ].value_counts().to_dict()
)

# ---------------------------------------------------------
# SAVE
# ---------------------------------------------------------

false_negative_validation.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("\nOutput:")
print(output_file)

print("\n" + "=" * 70)
print("PHASE 11B.2C COMPLETE")
print("=" * 70)
