import pandas as pd

from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

print("=" * 70)
print("ML PHASE 8 — LEAKAGE / ABLATION ANALYSIS")
print("=" * 70)

# ---------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------

df = pd.read_csv(
    "T2D_ML_protein_altering_dataset.tsv",
    sep="\t"
)

target = "Has_protein_altering"

y = df[target]

print("\nDataset shape:", df.shape)
print("\nTarget distribution:")
print(y.value_counts())

# ---------------------------------------------------------
# GWAS FEATURES
# ---------------------------------------------------------

gwas_features = [
    "Risk_Allele_Frequency",
    "MR_MEGA_Association_P",
    "Effective_Sample_Size",
    "OR",
    "CI_lower",
    "CI_upper",
    "log_OR",
    "CI_width",
    "SE_log_OR",
    "P_is_zero",
    "neg_log10_P",
    "MAF"
]

# ---------------------------------------------------------
# FULL VEP FEATURES
# ---------------------------------------------------------

full_vep_features = [
    "Gene_count",
    "Transcript_count",
    "Biotype_count",
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
    "Has_MODIFIER_impact"
]

# ---------------------------------------------------------
# SAFER VEP FEATURES
#
# Impact categories removed because they are too closely
# related to the definition of Has_protein_altering
# ---------------------------------------------------------

safe_vep_features = [
    "Gene_count",
    "Transcript_count",
    "Biotype_count",
    "Is_intronic",
    "Is_intergenic",
    "Is_upstream",
    "Is_downstream",
    "Is_5UTR",
    "Is_3UTR",
    "Is_regulatory",
    "Is_splice"
]

# ---------------------------------------------------------
# FEATURE SETS
# ---------------------------------------------------------

feature_sets = {
    "GWAS only":
        gwas_features,

    "GWAS + Full VEP":
        gwas_features + full_vep_features,

    "GWAS + Safe VEP":
        gwas_features + safe_vep_features
}

# ---------------------------------------------------------
# CROSS-VALIDATION SETUP
# ---------------------------------------------------------

cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=42
)

model = Pipeline([
    ("imputer", SimpleImputer(strategy="median")),
    ("scaler", StandardScaler()),
    (
        "classifier",
        LogisticRegression(
            class_weight="balanced",
            max_iter=2000,
            random_state=42
        )
    )
])

scoring = {
    "ROC_AUC": "roc_auc",
    "PR_AUC": "average_precision"
}

# ---------------------------------------------------------
# RUN ABLATION ANALYSIS
# ---------------------------------------------------------

results_list = []

for name, features in feature_sets.items():

    print("\n" + "-" * 70)
    print(name)
    print("-" * 70)

    X = df[features].copy()

    scores = cross_validate(
        model,
        X,
        y,
        cv=cv,
        scoring=scoring
    )

    roc_scores = scores["test_ROC_AUC"]
    pr_scores = scores["test_PR_AUC"]

    print("\nFeatures:", len(features))

    print("\nROC-AUC:")
    for i, score in enumerate(roc_scores, 1):
        print(f"Fold {i}: {score:.4f}")

    print(
        f"Mean: {roc_scores.mean():.4f} "
        f"± {roc_scores.std():.4f}"
    )

    print("\nPR-AUC:")
    for i, score in enumerate(pr_scores, 1):
        print(f"Fold {i}: {score:.4f}")

    print(
        f"Mean: {pr_scores.mean():.4f} "
        f"± {pr_scores.std():.4f}"
    )

    results_list.append({
        "Model": name,
        "Number_of_features": len(features),
        "Mean_ROC_AUC": roc_scores.mean(),
        "Std_ROC_AUC": roc_scores.std(),
        "Mean_PR_AUC": pr_scores.mean(),
        "Std_PR_AUC": pr_scores.std()
    })

# ---------------------------------------------------------
# SUMMARY
# ---------------------------------------------------------

results_df = pd.DataFrame(results_list)

print("\n" + "=" * 70)
print("ABLATION ANALYSIS SUMMARY")
print("=" * 70)

print(
    results_df.to_string(
        index=False
    )
)

results_df.to_csv(
    "T2D_ablation_analysis.tsv",
    sep="\t",
    index=False
)

print("\nSaved:")
print("T2D_ablation_analysis.tsv")

print("\n" + "=" * 70)
print("ML PHASE 8 COMPLETE")
print("=" * 70)
