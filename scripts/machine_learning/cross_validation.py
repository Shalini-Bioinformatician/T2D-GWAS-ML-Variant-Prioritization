import pandas as pd

from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

print("=" * 70)
print("ML PHASE 7 — STRATIFIED 5-FOLD CROSS-VALIDATION")
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

# ---------------------------------------------------------
# FEATURE SETS
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

vep_features = [
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

combined_features = gwas_features + vep_features

# ---------------------------------------------------------
# CROSS-VALIDATION
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
    "roc_auc": "roc_auc",
    "pr_auc": "average_precision"
}

# ---------------------------------------------------------
# FUNCTION TO EVALUATE A FEATURE SET
# ---------------------------------------------------------

def evaluate_feature_set(name, features):

    X = df[features].copy()

    results = cross_validate(
        model,
        X,
        y,
        cv=cv,
        scoring=scoring,
        return_train_score=False
    )

    roc_scores = results["test_roc_auc"]
    pr_scores = results["test_pr_auc"]

    print("\n" + "-" * 70)
    print(name)
    print("-" * 70)

    print("Number of features:", len(features))

    print("\nROC-AUC scores:")
    for i, score in enumerate(roc_scores, 1):
        print(f"Fold {i}: {score:.4f}")

    print(f"Mean ROC-AUC: {roc_scores.mean():.4f}")
    print(f"Std ROC-AUC:  {roc_scores.std():.4f}")

    print("\nPR-AUC scores:")
    for i, score in enumerate(pr_scores, 1):
        print(f"Fold {i}: {score:.4f}")

    print(f"Mean PR-AUC: {pr_scores.mean():.4f}")
    print(f"Std PR-AUC:  {pr_scores.std():.4f}")

    return {
        "Model": name,
        "Features": len(features),
        "Mean_ROC_AUC": roc_scores.mean(),
        "Std_ROC_AUC": roc_scores.std(),
        "Mean_PR_AUC": pr_scores.mean(),
        "Std_PR_AUC": pr_scores.std()
    }


# ---------------------------------------------------------
# MODEL A — GWAS ONLY
# ---------------------------------------------------------

gwas_results = evaluate_feature_set(
    "MODEL A — GWAS FEATURES ONLY",
    gwas_features
)

# ---------------------------------------------------------
# MODEL B — GWAS + VEP
# ---------------------------------------------------------

combined_results = evaluate_feature_set(
    "MODEL B — GWAS + VEP FEATURES",
    combined_features
)

# ---------------------------------------------------------
# FINAL COMPARISON
# ---------------------------------------------------------

summary = pd.DataFrame([
    gwas_results,
    combined_results
])

print("\n" + "=" * 70)
print("CROSS-VALIDATION SUMMARY")
print("=" * 70)

print(
    summary.to_string(
        index=False
    )
)

summary.to_csv(
    "T2D_cross_validation_results.tsv",
    sep="\t",
    index=False
)

print("\nSaved:")
print("T2D_cross_validation_results.tsv")

print("\n" + "=" * 70)
print("ML PHASE 7 COMPLETE")
print("=" * 70)
