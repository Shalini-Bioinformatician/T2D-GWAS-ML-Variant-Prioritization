import pandas as pd

print("=" * 70)
print("PHASE 11B.2B — FALSE-POSITIVE BIOLOGICAL VALIDATION SET")
print("=" * 70)

# ---------------------------------------------------------
# LOAD CANDIDATE PANEL
# ---------------------------------------------------------

input_file = "T2D_phase11B_candidate_panel.tsv"
output_file = "T2D_phase11B_false_positive_validation.tsv"

df = pd.read_csv(
    input_file,
    sep="\t"
)

print("\nCandidate panel shape:", df.shape)

# ---------------------------------------------------------
# SELECT HIGH-CONFIDENCE FALSE POSITIVES
# ---------------------------------------------------------

false_positive_validation = df[
    df["Candidate_Category"] == "High_Confidence_False_Positive"
].copy()

# ---------------------------------------------------------
# SORT BY MODEL PROBABILITY
# ---------------------------------------------------------

false_positive_validation = (
    false_positive_validation
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
print("HIGH-CONFIDENCE FALSE-POSITIVE CANDIDATES")
print("=" * 70)

print(
    false_positive_validation[
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
    len(false_positive_validation)
)

print(
    "Unique SNPs:",
    false_positive_validation["Index_SNV"].nunique()
)

print(
    "Actual protein-altering labels:",
    false_positive_validation[
        "Has_protein_altering"
    ].value_counts().to_dict()
)

print(
    "Predicted labels:",
    false_positive_validation[
        "Predicted_label"
    ].value_counts().to_dict()
)

# ---------------------------------------------------------
# SAVE
# ---------------------------------------------------------

false_positive_validation.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("\nOutput:")
print(output_file)

print("\n" + "=" * 70)
print("PHASE 11B.2B COMPLETE")
print("=" * 70)
