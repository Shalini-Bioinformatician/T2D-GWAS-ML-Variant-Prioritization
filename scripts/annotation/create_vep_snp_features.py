import pandas as pd

input_file = "ST4_combined_vep_annotations.tsv"
output_file = "ST4_vep_snp_features.tsv"

print("=" * 70)
print("VEP SNP-LEVEL FEATURE ENGINEERING")
print("=" * 70)

# ---------------------------------------------------------
# READ VEP ANNOTATIONS
# ---------------------------------------------------------

df = pd.read_csv(
    input_file,
    sep="\t",
    dtype=str
)

print("\nInput annotation rows:", len(df))
print("Unique SNPs:", df["SNP"].nunique())

# ---------------------------------------------------------
# CLEAN MISSING VALUES
# ---------------------------------------------------------

df = df.fillna("")

# ---------------------------------------------------------
# HELPER FUNCTION
# ---------------------------------------------------------

def normalize_consequences(values):
    """
    Combine VEP consequence values and normalize separators
    so individual consequence terms can be detected reliably.
    """
    text = " ".join(values.astype(str).tolist()).lower()

    # Normalize common separators
    text = text.replace(",", " ")
    text = text.replace("&", " ")
    text = text.replace(";", " ")

    return text

# ---------------------------------------------------------
# SNP-LEVEL FEATURE TABLE
# ---------------------------------------------------------

features = []

for snp, group in df.groupby("SNP"):

    consequences = " ".join(
        group["Consequence"].astype(str)
    ).lower()

    severe = " ".join(
        group["Most_severe_consequence"].astype(str)
    ).lower()

    impacts = set(
        group["Impact"]
        .astype(str)
        .str.upper()
    )

    genes = set(
        x for x in group["Gene_symbol"]
        if x != ""
    )

    transcripts = set(
        x for x in group["Transcript_ID"]
        if x != ""
    )

    biotypes = set(
        x for x in group["Biotype"]
        if x != ""
    )

    # -----------------------------------------------------
    # BASIC ANNOTATION
    # -----------------------------------------------------

    row = {
        "SNP": snp,

        "Gene_count":
            len(genes),

        "Transcript_count":
            len(transcripts),

        "Biotype_count":
            len(biotypes),

        # -------------------------------------------------
        # FUNCTIONAL CLASS FLAGS
        # -------------------------------------------------

        "Is_coding":
            int(
                any(
                    x in consequences
                    for x in [
                        "missense_variant",
                        "synonymous_variant",
                        "stop_gained",
                        "stop_lost",
                        "start_lost",
                        "frameshift_variant",
                        "inframe_insertion",
                        "inframe_deletion"
                    ]
                )
            ),

        "Is_missense":
            int("missense_variant" in consequences),

        "Is_synonymous":
            int("synonymous_variant" in consequences),

        "Is_stop_gained":
            int("stop_gained" in consequences),

        "Is_intronic":
            int("intron_variant" in consequences),

        "Is_intergenic":
            int("intergenic_variant" in severe),

        "Is_upstream":
            int("upstream_gene_variant" in consequences),

        "Is_downstream":
            int("downstream_gene_variant" in consequences),

        "Is_5UTR":
    	    int(
        	"5_prime_utr_variant" in consequences
        	or "5_prime_utr_variant" in severe
    	    ),

	"Is_3UTR":
   	    int(
        	"3_prime_utr_variant" in consequences
        	or "3_prime_utr_variant" in severe
   	    ),

        "Is_regulatory":
            int("regulatory_region_variant" in consequences),

        "Is_splice":
            int(
                any(
                    x in consequences
                    for x in [
                        "splice_donor",
                        "splice_acceptor",
                        "splice_region",
                        "splice_polypyrimidine"
                    ]
                )
            ),

        # -------------------------------------------------
        # IMPACT FEATURES
        # -------------------------------------------------

        "Has_HIGH_impact":
            int("HIGH" in impacts),

        "Has_MODERATE_impact":
            int("MODERATE" in impacts),

        "Has_LOW_impact":
            int("LOW" in impacts),

        "Has_MODIFIER_impact":
            int("MODIFIER" in impacts),

        # -------------------------------------------------
        # SEVERE FUNCTIONAL CONSEQUENCES
        # -------------------------------------------------

        "Has_protein_altering":
            int(
                any(
                    x in consequences
                    for x in [
                        "missense_variant",
                        "stop_gained",
                        "stop_lost",
                        "start_lost",
                        "frameshift_variant",
                        "inframe_insertion",
                        "inframe_deletion"
                    ]
                )
            ),

        # -------------------------------------------------
        # PRIMARY ANNOTATION
        # -------------------------------------------------

        "Most_severe_consequence":
            group["Most_severe_consequence"]
            .iloc[0],

        "Primary_gene":
            group["Gene_symbol"]
            .replace("", pd.NA)
            .dropna()
            .iloc[0]
            if group["Gene_symbol"]
               .replace("", pd.NA)
               .dropna()
               .shape[0] > 0
            else "",

        "Primary_biotype":
            group["Biotype"]
            .replace("", pd.NA)
            .dropna()
            .iloc[0]
            if group["Biotype"]
               .replace("", pd.NA)
               .dropna()
               .shape[0] > 0
            else ""
    }

    features.append(row)

# ---------------------------------------------------------
# CREATE DATAFRAME
# ---------------------------------------------------------

feature_df = pd.DataFrame(features)

# ---------------------------------------------------------
# SORT
# ---------------------------------------------------------

feature_df = feature_df.sort_values(
    "SNP"
).reset_index(drop=True)

# ---------------------------------------------------------
# SAVE
# ---------------------------------------------------------

feature_df.to_csv(
    output_file,
    sep="\t",
    index=False
)

# ---------------------------------------------------------
# QC
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("SNP-LEVEL FEATURE TABLE")
print("=" * 70)

print("Rows:", len(feature_df))
print("Columns:", len(feature_df.columns))
print("Unique SNPs:", feature_df["SNP"].nunique())

print("\nColumns:")
for col in feature_df.columns:
    print(" ", col)

print("\nMissing values:")
print(
    feature_df.isna().sum()
)

print("\nFunctional feature counts:")

feature_columns = [
    "Is_coding",
    "Is_missense",
    "Is_synonymous",
    "Is_stop_gained",
    "Is_intronic",
    "Is_intergenic",
    "Is_upstream",
    "Is_downstream",
    "Is_5UTR",
    "Is_3UTR",
    "Is_regulatory",
    "Is_splice",
    "Has_HIGH_impact",
    "Has_MODERATE_impact",
    "Has_LOW_impact",
    "Has_MODIFIER_impact",
    "Has_protein_altering"
]

for col in feature_columns:
    print(
        f"{col:30s}: "
        f"{feature_df[col].sum()}"
    )

print("\nMost severe consequence:")
print(
    feature_df[
        "Most_severe_consequence"
    ].value_counts()
)

print("\n" + "=" * 70)
print("FEATURE ENGINEERING COMPLETE")
print("=" * 70)

print("\nOutput file:")
print(output_file)
