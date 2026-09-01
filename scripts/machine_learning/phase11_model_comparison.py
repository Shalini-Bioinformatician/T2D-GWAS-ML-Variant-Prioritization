import pandas as pd

from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier

print("=" * 70)
print("ML PHASE 11 — LEAKAGE-SAFE MODEL COMPARISON")
print("=" * 70)

# ---------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------

df = pd.read_csv(
    "T2D_GWAS_ML_features.tsv",
    sep="\t"
)

print("\nDataset shape:")
print(df.shape)

# ---------------------------------------------------------
# TARGET
# ---------------------------------------------------------

target = "Has_protein_altering"

y = df[target]

# ---------------------------------------------------------
# LEAKAGE-SAFE FEATURES
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
    "Is_splice",
    "Has_LOW_impact",
    "Has_MODIFIER_impact"
]

features = gwas_features + safe_vep_features

X = df[features]

print("\nNumber of features:")
print(len(features))

print("\nTarget distribution:")
print(y.value_counts())

# ---------------------------------------------------------
# CROSS-VALIDATION
# ---------------------------------------------------------

cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=42
)

scoring = {
    "roc_auc": "roc_auc",
    "average_precision": "average_precision",
    "precision": "precision",
    "recall": "recall",
    "f1": "f1"
}

# ---------------------------------------------------------
# MODELS
# ---------------------------------------------------------

models = {

    "Logistic Regression":
    Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        (
            "model",
            LogisticRegression(
                class_weight="balanced",
                max_iter=5000,
                random_state=42
            )
        )
    ]),

    "Random Forest":
    Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        (
            "model",
            RandomForestClassifier(
                n_estimators=500,
                class_weight="balanced",
                random_state=42,
                n_jobs=-1
            )
        )
    ]),

    "HistGradientBoosting":
    Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        (
            "model",
            HistGradientBoostingClassifier(
                max_iter=300,
                learning_rate=0.05,
                max_leaf_nodes=15,
                random_state=42
            )
        )
    ])
}

# ---------------------------------------------------------
# RUN CROSS-VALIDATION
# ---------------------------------------------------------

results = []

for model_name, model in models.items():

    print("\n" + "-" * 70)
    print(model_name)
    print("-" * 70)

    scores = cross_validate(
        model,
        X,
        y,
        cv=cv,
        scoring=scoring,
        return_train_score=False,
        n_jobs=-1
    )

    roc_auc_mean = scores["test_roc_auc"].mean()
    roc_auc_std = scores["test_roc_auc"].std()

    pr_auc_mean = scores["test_average_precision"].mean()
    pr_auc_std = scores["test_average_precision"].std()

    precision_mean = scores["test_precision"].mean()
    recall_mean = scores["test_recall"].mean()
    f1_mean = scores["test_f1"].mean()

    print(f"\nROC-AUC: {roc_auc_mean:.4f} ± {roc_auc_std:.4f}")
    print(f"PR-AUC: {pr_auc_mean:.4f} ± {pr_auc_std:.4f}")
    print(f"Precision: {precision_mean:.4f}")
    print(f"Recall: {recall_mean:.4f}")
    print(f"F1-score: {f1_mean:.4f}")

    results.append({
        "Model": model_name,
        "Features": len(features),
        "Mean_ROC_AUC": roc_auc_mean,
        "Std_ROC_AUC": roc_auc_std,
        "Mean_PR_AUC": pr_auc_mean,
        "Std_PR_AUC": pr_auc_std,
        "Mean_Precision": precision_mean,
        "Mean_Recall": recall_mean,
        "Mean_F1": f1_mean
    })

# ---------------------------------------------------------
# SUMMARY
# ---------------------------------------------------------

results_df = pd.DataFrame(results)

results_df = results_df.sort_values(
    by="Mean_PR_AUC",
    ascending=False
)

print("\n" + "=" * 70)
print("MODEL COMPARISON SUMMARY")
print("=" * 70)

print(
    results_df.to_string(
        index=False
    )
)

# ---------------------------------------------------------
# SAVE RESULTS
# ---------------------------------------------------------

results_df.to_csv(
    "T2D_phase11_model_comparison.tsv",
    sep="\t",
    index=False
)

print("\nSaved:")
print("T2D_phase11_model_comparison.tsv")

print("\n" + "=" * 70)
print("ML PHASE 11 COMPLETE")
print("=" * 70)
