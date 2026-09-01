# =============================================================================
# PHASE 13 — FIGURE 4
# MODEL VS BIOLOGICAL EVIDENCE DISCORDANCE
# =============================================================================

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt


OUTPUT_DIR = "../phase13_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)


print("=" * 80)
print("PHASE 13 — FIGURE 4: MODEL VS BIOLOGICAL EVIDENCE DISCORDANCE")
print("=" * 80)


# =============================================================================
# LOCATE DISCORDANCE FILE
# =============================================================================

candidates = [
    "../T2D_phase11D_model_biology_discordance.tsv",
    "T2D_phase11D_model_biology_discordance.tsv"
]

discordance_file = None

for file in candidates:
    if os.path.exists(file):
        discordance_file = file
        break


if discordance_file is None:

    matches = glob.glob(
        "../**/T2D_phase11D_model_biology_discordance.tsv",
        recursive=True
    )

    if matches:
        discordance_file = matches[0]


if discordance_file is None:
    raise FileNotFoundError(
        "Could not locate T2D_phase11D_model_biology_discordance.tsv"
    )


print(f"\nDiscordance file found: {discordance_file}")


# =============================================================================
# LOAD DATA
# =============================================================================

df = pd.read_csv(
    discordance_file,
    sep="\t"
)

print(f"Rows: {df.shape[0]}")
print(f"Unique SNPs: {df['Index_SNV'].nunique()}")

print("\nDiscordance type counts:")

discordance_counts = (
    df["Discordance_Type"]
    .value_counts()
)

print(discordance_counts)


# =============================================================================
# FIGURE 4A — DISCORDANCE TYPE COUNTS
# =============================================================================

plt.figure(figsize=(8, 6))

plt.bar(
    discordance_counts.index,
    discordance_counts.values
)

plt.xlabel("Discordance Type")
plt.ylabel("Number of Candidates")
plt.title("Model–Biology Discordance Among Final T2D Candidates")

plt.tight_layout()

output1 = os.path.join(
    OUTPUT_DIR,
    "Figure_4_Model_Biology_Discordance.png"
)

plt.savefig(
    output1,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(f"\nSaved: {output1}")


# =============================================================================
# FIGURE 4B — CANDIDATE-LEVEL DISCORDANCE
# =============================================================================

plot_df = df.copy()

plot_df["Candidate_Label"] = (
    plot_df["Index_SNV"].astype(str)
    + "\n"
    + plot_df["Primary_gene"].astype(str)
)

plot_df = plot_df.sort_values(
    "Integrated_Evidence_Score",
    ascending=False
)

plt.figure(figsize=(12, 8))

plt.barh(
    plot_df["Candidate_Label"],
    plot_df["Integrated_Evidence_Score"]
)

plt.xlabel("Integrated Evidence Score")
plt.ylabel("Discordant Candidate")
plt.title(
    "Integrated Evidence Scores of Model–Biology Discordant Candidates"
)

plt.gca().invert_yaxis()

plt.tight_layout()

output2 = os.path.join(
    OUTPUT_DIR,
    "Figure_4B_Discordant_Candidate_Evidence.png"
)

plt.savefig(
    output2,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print(f"Saved: {output2}")


# =============================================================================
# SAVE DISCORDANCE SUMMARY
# =============================================================================

summary = (
    df.groupby("Discordance_Type")
    .size()
    .reset_index(name="Candidate_Count")
)

summary_file = os.path.join(
    OUTPUT_DIR,
    "T2D_Model_Biology_Discordance_Summary.tsv"
)

summary.to_csv(
    summary_file,
    sep="\t",
    index=False
)

print(f"Saved: {summary_file}")


# =============================================================================
# COMPLETE
# =============================================================================

print("\n" + "=" * 80)
print("FIGURE 4 COMPLETE")
print("=" * 80)

print(f"\nTotal discordant candidates: {len(df)}")

print("\nFinal discordance summary:")
print(summary.to_string(index=False))
