import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    confusion_matrix,
    classification_report,
    roc_auc_score,
    average_precision_score,
    precision_recall_curve,
    roc_curve
)

print("=" * 70)
print("ML PHASE 9 — FINAL LEAKAGE-SAFE MODEL")
print("=" * 70)

# ---------------------------------------------------------
# LOAD DATASET
# ---------------------------------------------------------

input_file = "T2D_ML_protein_altering_dataset.tsv"

df = pd.read_csv(
    input_file,
    sep="\t"
)

print("\nDataset shape:", df.shape)

# ---------------------------------------------------------
# TARGET
# ---------------------------------------------------------

target = "Has_protein_altering"

y = df[target]

# ---------------------------------------------------------
# DEFINE LEAKAGE-SAFE FEATURES
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
    "Has_MODIFIER_impact",
    "Primary_biotype"
]

# Keep only numeric leakage-safe features
safe_vep_features = [
    feature for feature in safe_vep_features
    if feature in df.columns
    and pd.api.types.is_numeric_dtype(df[feature])
]

feature_columns = gwas_features + safe_vep_features

print("\nNumber of GWAS features:", len(gwas_features))
print("Number of safe VEP features:", len(safe_vep_features))
print("Total features:", len(feature_columns))

print("\nFeatures used:")
for i, feature in enumerate(feature_columns, start=1):
    print(f"{i:2d}. {feature}")

X = df[feature_columns]

# ---------------------------------------------------------
# RECREATE PHASE 3 STRATIFIED SPLIT
# ---------------------------------------------------------

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.20,
    stratify=y,
    random_state=42
)

print("\n" + "=" * 70)
print("TRAIN / TEST SPLIT")
print("=" * 70)

print("\nTraining samples:", len(X_train))
print("Test samples:", len(X_test))

print("\nTraining positives:", y_train.sum())
print("Test positives:", y_test.sum())

# ---------------------------------------------------------
# BUILD LEAKAGE-SAFE PIPELINE
# ---------------------------------------------------------

model = Pipeline([
    (
        "imputer",
        SimpleImputer(strategy="median")
    ),
    (
        "scaler",
        StandardScaler()
    ),
    (
        "classifier",
        LogisticRegression(
            class_weight="balanced",
            max_iter=5000,
            random_state=42
        )
    )
])

# ---------------------------------------------------------
# TRAIN MODEL
# ---------------------------------------------------------

print("\nTraining final leakage-safe Logistic Regression...")

model.fit(
    X_train,
    y_train
)

print("Training complete.")

# ---------------------------------------------------------
# PREDICTIONS
# ---------------------------------------------------------

y_pred = model.predict(X_test)

y_prob = model.predict_proba(X_test)[:, 1]

# ---------------------------------------------------------
# EVALUATION
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("CONFUSION MATRIX")
print("=" * 70)

print(
    confusion_matrix(
        y_test,
        y_pred
    )
)

print("\n" + "=" * 70)
print("CLASSIFICATION REPORT")
print("=" * 70)

print(
    classification_report(
        y_test,
        y_pred,
        digits=4
    )
)

roc_auc = roc_auc_score(
    y_test,
    y_prob
)

pr_auc = average_precision_score(
    y_test,
    y_prob
)

print("\n" + "=" * 70)
print("FINAL PERFORMANCE")
print("=" * 70)

print(f"\nROC-AUC: {roc_auc:.4f}")
print(f"PR-AUC: {pr_auc:.4f}")

# ---------------------------------------------------------
# SAVE TEST PREDICTIONS
# ---------------------------------------------------------

results = X_test.copy()

results["True_label"] = y_test.values
results["Predicted_label"] = y_pred
results["Protein_altering_probability"] = y_prob

results.to_csv(
    "T2D_final_safe_model_predictions.tsv",
    sep="\t",
    index=False
)

# ---------------------------------------------------------
# SAVE MODEL COEFFICIENTS
# ---------------------------------------------------------

coefficients = pd.DataFrame({
    "Feature": feature_columns,
    "Coefficient": model.named_steps[
        "classifier"
    ].coef_[0]
})

coefficients["Absolute_coefficient"] = (
    coefficients["Coefficient"].abs()
)

coefficients = coefficients.sort_values(
    "Absolute_coefficient",
    ascending=False
)

coefficients.to_csv(
    "T2D_final_safe_model_coefficients.tsv",
    sep="\t",
    index=False
)

print("\nSaved predictions:")
print("T2D_final_safe_model_predictions.tsv")

print("\nSaved feature coefficients:")
print("T2D_final_safe_model_coefficients.tsv")

print("\n" + "=" * 70)
print("ML PHASE 9 COMPLETE")
print("=" * 70)
