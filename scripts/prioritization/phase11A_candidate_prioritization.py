import pandas as pd

from sklearn.model_selection import train_test_split

print("=" * 70)
print("ML PHASE 11A — BIOLOGICAL CANDIDATE PRIORITIZATION")
print("=" * 70)

# ---------------------------------------------------------
# LOAD COMPLETE DATASET
# ---------------------------------------------------------

input_file = "T2D_ML_protein_altering_dataset.tsv"
prediction_file = "T2D_final_safe_model_predictions.tsv"

df = pd.read_csv(
    input_file,
    sep="\t"
)

predictions = pd.read_csv(
    prediction_file,
    sep="\t"
)

print("\nComplete dataset shape:", df.shape)
print("Prediction file shape:", predictions.shape)


# ---------------------------------------------------------
# RECREATE EXACT PHASE 3 / PHASE 9 SPLIT
# ---------------------------------------------------------

target = "Has_protein_altering"

X = df.drop(
    columns=[target]
)

y = df[target]

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.20,
    stratify=y,
    random_state=42
)

print("\nRecreated test set size:", len(X_test))
print("Prediction rows:", len(predictions))


# ---------------------------------------------------------
# QC CHECK
# ---------------------------------------------------------

if len(X_test) != len(predictions):
    raise ValueError(
        "ERROR: Test set size does not match prediction file size."
    )

print("\nTest set and prediction file row counts match.")


# ---------------------------------------------------------
# GET ORIGINAL BIOLOGICAL METADATA
# ---------------------------------------------------------

metadata_columns = [
    "Index_SNV",
    "Locus",
    "Chromosome",
    "Position_b37",
    "Risk_Allele",
    "Other_Allele",
    "Risk_Allele_Frequency",
    "OR",
    "MR_MEGA_Association_P",
    "neg_log10_P",
    "MAF",
    "MAF_category",
    "Primary_gene",
    "Primary_biotype",
    "Most_severe_consequence",
    "Has_protein_altering"
]

metadata = df.loc[
    X_test.index,
    metadata_columns
].copy()


# ---------------------------------------------------------
# ATTACH MODEL PREDICTIONS
# ---------------------------------------------------------

metadata["Predicted_label"] = (
    predictions["Predicted_label"].values
)

metadata["Protein_altering_probability"] = (
    predictions[
        "Protein_altering_probability"
    ].values
)


# ---------------------------------------------------------
# SORT BY PREDICTED PROBABILITY
# ---------------------------------------------------------

metadata = metadata.sort_values(
    "Protein_altering_probability",
    ascending=False
).reset_index(drop=True)


# ---------------------------------------------------------
# OVERALL PRIORITIZATION SUMMARY
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("TOP 20 PRIORITIZED TEST-SET SNPs")
print("=" * 70)

print(
    metadata[
        [
            "Index_SNV",
            "Primary_gene",
            "OR",
            "neg_log10_P",
            "MAF",
            "Most_severe_consequence",
            "Has_protein_altering",
            "Predicted_label",
            "Protein_altering_probability"
        ]
    ]
    .head(20)
    .to_string(index=False)
)


# ---------------------------------------------------------
# TRUE POSITIVES
# ---------------------------------------------------------

true_positives = metadata[
    (metadata["Has_protein_altering"] == 1) &
    (metadata["Predicted_label"] == 1)
].copy()

print("\n" + "=" * 70)
print("TRUE POSITIVES")
print("=" * 70)

print("Number of true positives:", len(true_positives))

print(
    true_positives[
        [
            "Index_SNV",
            "Primary_gene",
            "OR",
            "neg_log10_P",
            "Most_severe_consequence",
            "Protein_altering_probability"
        ]
    ]
    .to_string(index=False)
)


# ---------------------------------------------------------
# FALSE NEGATIVES
# ---------------------------------------------------------

false_negatives = metadata[
    (metadata["Has_protein_altering"] == 1) &
    (metadata["Predicted_label"] == 0)
].copy()

print("\n" + "=" * 70)
print("FALSE NEGATIVES")
print("=" * 70)

print("Number of false negatives:", len(false_negatives))

print(
    false_negatives[
        [
            "Index_SNV",
            "Primary_gene",
            "OR",
            "neg_log10_P",
            "Most_severe_consequence",
            "Protein_altering_probability"
        ]
    ]
    .to_string(index=False)
)


# ---------------------------------------------------------
# FALSE POSITIVES
# ---------------------------------------------------------

false_positives = metadata[
    (metadata["Has_protein_altering"] == 0) &
    (metadata["Predicted_label"] == 1)
].copy()

print("\n" + "=" * 70)
print("TOP FALSE POSITIVE CANDIDATES")
print("=" * 70)

print("Number of false positives:", len(false_positives))

print(
    false_positives[
        [
            "Index_SNV",
            "Primary_gene",
            "OR",
            "neg_log10_P",
            "MAF",
            "Most_severe_consequence",
            "Protein_altering_probability"
        ]
    ]
    .sort_values(
        "Protein_altering_probability",
        ascending=False
    )
    .head(20)
    .to_string(index=False)
)


# ---------------------------------------------------------
# SAVE OUTPUTS
# ---------------------------------------------------------

metadata.to_csv(
    "T2D_phase11A_prioritized_candidates.tsv",
    sep="\t",
    index=False
)

true_positives.to_csv(
    "T2D_phase11A_true_positives.tsv",
    sep="\t",
    index=False
)

false_negatives.to_csv(
    "T2D_phase11A_false_negatives.tsv",
    sep="\t",
    index=False
)

false_positives.to_csv(
    "T2D_phase11A_false_positive_candidates.tsv",
    sep="\t",
    index=False
)


print("\n" + "=" * 70)
print("SAVED FILES")
print("=" * 70)

print("T2D_phase11A_prioritized_candidates.tsv")
print("T2D_phase11A_true_positives.tsv")
print("T2D_phase11A_false_negatives.tsv")
print("T2D_phase11A_false_positive_candidates.tsv")

print("\n" + "=" * 70)
print("ML PHASE 11A COMPLETE")
print("=" * 70)
