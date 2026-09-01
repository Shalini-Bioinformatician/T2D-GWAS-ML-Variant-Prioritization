
import pandas as pd

input_file = "T2D_phase11A_prioritized_candidates.tsv"
output_file = "T2D_phase11B_candidate_panel.tsv"

df = pd.read_csv(input_file, sep="\t")

# Top 10 candidates by model probability
top10 = df.sort_values(
    "Protein_altering_probability",
    ascending=False
).head(10).copy()

# All true positives
tp = df[df["Has_protein_altering"] == 1].copy()

# All false negatives
fn = df[
    (df["Has_protein_altering"] == 1) &
    (df["Predicted_label"] == 0)
].copy()

# Top 7 false positives
fp = df[
    (df["Has_protein_altering"] == 0) &
    (df["Predicted_label"] == 1)
].sort_values(
    "Protein_altering_probability",
    ascending=False
).head(7).copy()

# Combine
candidate = pd.concat(
    [top10, tp, fn, fp],
    ignore_index=True
)

# Remove duplicate SNPs
candidate = candidate.drop_duplicates(
    subset=["Index_SNV"]
)

# Add candidate category
candidate["Candidate_Category"] = "Other"

candidate.loc[
    candidate["Index_SNV"].isin(tp["Index_SNV"]),
    "Candidate_Category"
] = "True_Positive"

candidate.loc[
    candidate["Index_SNV"].isin(fn["Index_SNV"]),
    "Candidate_Category"
] = "False_Negative"

candidate.loc[
    candidate["Index_SNV"].isin(fp["Index_SNV"]),
    "Candidate_Category"
] = "High_Confidence_False_Positive"

# Sort by probability
candidate = candidate.sort_values(
    "Protein_altering_probability",
    ascending=False
).reset_index(drop=True)

candidate.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("=" * 70)
print("PHASE 11B.1 — CANDIDATE VALIDATION PANEL")
print("=" * 70)

print("\nTotal candidates:", len(candidate))

print("\nCandidate categories:")
print(candidate["Candidate_Category"].value_counts())

print("\nCandidates:")
print(
    candidate[
        [
            "Index_SNV",
            "Primary_gene",
            "Most_severe_consequence",
            "Protein_altering_probability",
            "Has_protein_altering",
            "Predicted_label",
            "Candidate_Category"
        ]
    ].to_string(index=False)
)

print("\nOutput:")
print(output_file)

print("\n" + "=" * 70)
print("PHASE 11B.1 COMPLETE")
print("=" * 70)

