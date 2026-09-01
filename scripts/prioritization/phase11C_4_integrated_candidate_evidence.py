import pandas as pd
import numpy as np

print("=" * 70)
print("PHASE 11C.4 — INTEGRATED CANDIDATE EVIDENCE RANKING")
print("=" * 70)


# =========================================================
# FILES
# =========================================================

ml_file = "../../T2D_phase11B_comparative_candidate_analysis.tsv"

gwas_summary_file = "phase11C_external_validation_summary.tsv"

gwas_ranking_file = "phase11C_external_GWAS_evidence_ranking.tsv"


# =========================================================
# LOAD DATA
# =========================================================

print("\nLoading ML/biological validation data...")
ml = pd.read_csv(
    ml_file,
    sep="\t"
)

print("Rows:", len(ml))
print("Unique SNPs:", ml["Index_SNV"].nunique())


print("\nLoading GWAS Catalog summary...")
gwas = pd.read_csv(
    gwas_summary_file,
    sep="\t"
)

print("Rows:", len(gwas))
print("Unique SNPs:", gwas["Query_SNP"].nunique())


print("\nLoading GWAS evidence ranking...")
gwas_rank = pd.read_csv(
    gwas_ranking_file,
    sep="\t"
)

print("Rows:", len(gwas_rank))
print("Unique SNPs:", gwas_rank["Query_SNP"].nunique())


# =========================================================
# BASIC QC
# =========================================================

print("\n" + "=" * 70)
print("INPUT QC")
print("=" * 70)

if ml["Index_SNV"].nunique() != len(ml):
    raise ValueError("Duplicate SNPs found in ML validation data.")

if gwas["Query_SNP"].nunique() != len(gwas):
    raise ValueError("Duplicate SNPs found in GWAS summary.")

if gwas_rank["Query_SNP"].nunique() != len(gwas_rank):
    raise ValueError("Duplicate SNPs found in GWAS ranking.")


ml_snps = set(ml["Index_SNV"])
gwas_snps = set(gwas["Query_SNP"])

if ml_snps != gwas_snps:
    print("\nWARNING: SNP sets do not match.")

    print("\nSNPs in ML but missing from GWAS:")
    print(sorted(ml_snps - gwas_snps))

    print("\nSNPs in GWAS but missing from ML:")
    print(sorted(gwas_snps - ml_snps))

    raise ValueError("SNP sets do not match.")

print("\nQC PASS")
print("ML SNPs:", len(ml_snps))
print("GWAS SNPs:", len(gwas_snps))
print("Matching SNPs:", len(ml_snps))


# =========================================================
# SELECT GWAS EVIDENCE COLUMNS
# =========================================================

gwas_columns = [
    "Query_SNP",
    "Total_GWAS_Associations",
    "Direct_T2D_Associations",
    "Glycemic_Associations",
    "Metabolic_Associations",
    "Lipid_Cardiovascular_Associations",
    "Has_Direct_T2D_Evidence",
    "Has_Glycemic_Evidence",
    "Has_Metabolic_Evidence",
    "Has_Lipid_Cardiovascular_Evidence",
    "Has_T2D_Related_Evidence"
]

gwas_rank_columns = [
    "Query_SNP",
    "External_GWAS_Evidence_Score",
    "External_Evidence_Tier",
    "External_GWAS_Rank"
]


# =========================================================
# MERGE
# =========================================================

print("\n" + "=" * 70)
print("MERGING ML + EXTERNAL GWAS EVIDENCE")
print("=" * 70)

integrated = ml.merge(
    gwas[gwas_columns],
    left_on="Index_SNV",
    right_on="Query_SNP",
    how="inner",
    validate="one_to_one"
)

integrated = integrated.merge(
    gwas_rank[gwas_rank_columns],
    left_on="Index_SNV",
    right_on="Query_SNP",
    how="left",
    validate="one_to_one",
    suffixes=("", "_ranking")
)

# Remove duplicate SNP column
integrated = integrated.drop(
    columns=["Query_SNP", "Query_SNP_ranking"],
    errors="ignore"
)

print("\nIntegrated rows:", len(integrated))
print("Integrated unique SNPs:", integrated["Index_SNV"].nunique())


if len(integrated) != 17:
    raise ValueError(
        "ERROR: Integrated dataset does not contain all 17 candidates."
    )


# =========================================================
# ML EVIDENCE COMPONENT
# =========================================================

# Model probability is already bounded between 0 and 1.

integrated["ML_Evidence_Score"] = (
    integrated["Protein_altering_probability"]
)


# =========================================================
# FUNCTIONAL EVIDENCE COMPONENT
# =========================================================

# Functional consequence is converted into a transparent
# biological evidence score.

functional_scores = {
    "stop_gained": 1.00,
    "missense_variant": 0.90,
    "splice_donor_region_variant": 0.90,
    "splice_region_variant": 0.80,
    "splice_polypyrimidine_tract_variant": 0.80,
    "synonymous_variant": 0.20,
    "3_prime_UTR_variant": 0.20,
    "5_prime_UTR_variant": 0.20,
    "upstream_gene_variant": 0.15,
    "downstream_gene_variant": 0.15,
    "intron_variant": 0.10,
    "intergenic_variant": 0.05,
    "non_coding_transcript_exon_variant": 0.10,
    "regulatory_region_variant": 0.30,
    "mature_miRNA_variant": 0.30
}

integrated["Functional_Evidence_Score"] = (
    integrated["Most_severe_consequence"]
    .map(functional_scores)
    .fillna(0.05)
)


# =========================================================
# GWAS SIGNIFICANCE COMPONENT
# =========================================================

# Convert -log10(P) into a bounded 0-1 evidence score.
#
# 7.3 corresponds approximately to genome-wide significance.
# Values above 7.3 are capped so extremely small P-values
# do not dominate the integrated score.

integrated["GWAS_Significance_Score"] = (
    integrated["neg_log10_P"]
    .clip(lower=0, upper=50)
    / 50
)


# =========================================================
# EXTERNAL T2D EVIDENCE COMPONENT
# =========================================================

# Direct T2D evidence is given the strongest external weight.
# log1p prevents very highly studied SNPs from dominating.

integrated["Direct_T2D_Evidence_Score"] = (
    np.log1p(
        integrated["Direct_T2D_Associations"]
    ) /
    np.log1p(
        integrated["Direct_T2D_Associations"].max()
    )
)


# Related metabolic evidence
integrated["Related_Metabolic_Evidence_Score"] = (
    (
        0.5 * np.log1p(
            integrated["Glycemic_Associations"]
        )
        +
        0.3 * np.log1p(
            integrated["Metabolic_Associations"]
        )
        +
        0.2 * np.log1p(
            integrated["Lipid_Cardiovascular_Associations"]
        )
    )
)


max_related = (
    integrated["Related_Metabolic_Evidence_Score"].max()
)

if max_related > 0:
    integrated["Related_Metabolic_Evidence_Score"] = (
        integrated["Related_Metabolic_Evidence_Score"]
        / max_related
    )


# =========================================================
# INTEGRATED EVIDENCE SCORE
# =========================================================

integrated["Integrated_Evidence_Score"] = (
    0.25 * integrated["ML_Evidence_Score"]
    +
    0.25 * integrated["Functional_Evidence_Score"]
    +
    0.15 * integrated["GWAS_Significance_Score"]
    +
    0.25 * integrated["Direct_T2D_Evidence_Score"]
    +
    0.10 * integrated["Related_Metabolic_Evidence_Score"]
)


# =========================================================
# FINAL EVIDENCE CATEGORY
# =========================================================

def classify_integrated_evidence(row):

    ml = row["ML_Evidence_Score"]
    functional = row["Functional_Evidence_Score"]
    direct = row["Direct_T2D_Associations"]
    score = row["Integrated_Evidence_Score"]

    # Strong biological + external support
    if (
        score >= 0.70
        and functional >= 0.80
        and direct >= 2
        and ml >= 0.70
    ):
        return "Tier_1_High_Priority"

    # Strong external T2D evidence with reasonable
    # biological/model support
    elif (
        score >= 0.55
        and direct >= 2
    ):
        return "Tier_2_Strong_Candidate"

    # Evidence supported by either model/functional
    # evidence or direct T2D evidence
    elif (
        direct >= 1
        or functional >= 0.80
        or ml >= 0.80
    ):
        return "Tier_3_Promising_Candidate"

    else:
        return "Tier_4_Lower_Priority"


integrated["Integrated_Evidence_Tier"] = (
    integrated.apply(
        classify_integrated_evidence,
        axis=1
    )
)


# =========================================================
# FINAL RANK
# =========================================================

integrated = integrated.sort_values(
    "Integrated_Evidence_Score",
    ascending=False
).reset_index(drop=True)

integrated["Integrated_Evidence_Rank"] = (
    integrated["Integrated_Evidence_Score"]
    .rank(
        ascending=False,
        method="min"
    )
    .astype(int)
)


# =========================================================
# DISPLAY FINAL RANKING
# =========================================================

print("\n" + "=" * 70)
print("INTEGRATED CANDIDATE RANKING")
print("=" * 70)

display_columns = [
    "Integrated_Evidence_Rank",
    "Index_SNV",
    "Primary_gene",
    "Validation_Group",
    "Most_severe_consequence",
    "Protein_altering_probability",
    "Has_protein_altering",
    "Predicted_label",
    "OR",
    "neg_log10_P",
    "MAF",
    "Direct_T2D_Associations",
    "Glycemic_Associations",
    "Metabolic_Associations",
    "Lipid_Cardiovascular_Associations",
    "External_GWAS_Evidence_Score",
    "External_Evidence_Tier",
    "Integrated_Evidence_Score",
    "Integrated_Evidence_Tier"
]

print(
    integrated[
        display_columns
    ].to_string(index=False)
)


# =========================================================
# TIER COUNTS
# =========================================================

print("\n" + "=" * 70)
print("INTEGRATED EVIDENCE TIER COUNTS")
print("=" * 70)

print(
    integrated[
        "Integrated_Evidence_Tier"
    ]
    .value_counts()
    .sort_index()
)


# =========================================================
# MODEL / EXTERNAL DISCORDANCE
# =========================================================

print("\n" + "=" * 70)
print("MODEL / EXTERNAL EVIDENCE DISCORDANCE")
print("=" * 70)

discordant = integrated[
    (
        (
            integrated["Protein_altering_probability"] >= 0.75
        )
        &
        (
            integrated["Direct_T2D_Associations"] == 0
        )
    )
    |
    (
        (
            integrated["Protein_altering_probability"] < 0.50
        )
        &
        (
            integrated["Direct_T2D_Associations"] >= 2
        )
    )
].copy()

print(
    "Discordant candidates:",
    len(discordant)
)

if len(discordant) > 0:
    print(
        discordant[
            [
                "Index_SNV",
                "Primary_gene",
                "Most_severe_consequence",
                "Protein_altering_probability",
                "Direct_T2D_Associations",
                "External_Evidence_Tier",
                "Validation_Group",
                "Integrated_Evidence_Score",
                "Integrated_Evidence_Tier"
            ]
        ]
        .sort_values(
            "Integrated_Evidence_Score",
            ascending=False
        )
        .to_string(index=False)
    )


# =========================================================
# SAVE COMPLETE INTEGRATED TABLE
# =========================================================

output_file = (
    "phase11C_integrated_candidate_evidence.tsv"
)

integrated.to_csv(
    output_file,
    sep="\t",
    index=False
)


# =========================================================
# SAVE COMPACT RANKING
# =========================================================

ranking_file = (
    "phase11C_integrated_candidate_ranking.tsv"
)

integrated[
    [
        "Integrated_Evidence_Rank",
        "Index_SNV",
        "Primary_gene",
        "Validation_Group",
        "Most_severe_consequence",
        "Protein_altering_probability",
        "OR",
        "neg_log10_P",
        "MAF",
        "Direct_T2D_Associations",
        "Glycemic_Associations",
        "Metabolic_Associations",
        "Lipid_Cardiovascular_Associations",
        "External_GWAS_Evidence_Score",
        "External_Evidence_Tier",
        "Integrated_Evidence_Score",
        "Integrated_Evidence_Tier"
    ]
].to_csv(
    ranking_file,
    sep="\t",
    index=False
)


# =========================================================
# COMPLETE
# =========================================================

print("\n" + "=" * 70)
print("SAVED FILES")
print("=" * 70)

print(output_file)
print(ranking_file)

print("\n" + "=" * 70)
print("PHASE 11C.4 COMPLETE")
print("=" * 70)