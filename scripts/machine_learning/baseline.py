import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    confusion_matrix,
    classification_report,
    roc_auc_score,
    average_precision_score
)

# =========================================================
# ML PHASE 4 — BASELINE LOGISTIC REGRESSION
# =========================================================

print("=" * 70)
print("ML PHASE 4 — BASELINE LOGISTIC REGRESSION")
print("=" * 70)

# ---------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------

input_file = "T2D_ML_protein_altering_dataset.tsv"

df = pd.read_csv(
    input_file,
    sep="\t"
)

target = "Has_protein_altering"

features = [
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
    "MAF",
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

X = df[features].copy()
y = df[target].copy()

# ---------------------------------------------------------
# TRAIN / TEST SPLIT
# ---------------------------------------------------------

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.20,
    random_state=42,
    stratify=y
)

print("\nTraining samples:", len(X_train))
print("Test samples:", len(X_test))

print("\nTraining positives:", y_train.sum())
print("Test positives:", y_test.sum())

# ---------------------------------------------------------
# MODEL PIPELINE
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
            max_iter=2000,
            random_state=42
        )
    )
])

# ---------------------------------------------------------
# TRAIN
# ---------------------------------------------------------

print("\nTraining Logistic Regression...")

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
# CONFUSION MATRIX
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("CONFUSION MATRIX")
print("=" * 70)

cm = confusion_matrix(
    y_test,
    y_pred
)

print(cm)

# ---------------------------------------------------------
# CLASSIFICATION REPORT
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("CLASSIFICATION REPORT")
print("=" * 70)

print(
    classification_report(
        y_test,
        y_pred,
        digits=4,
        zero_division=0
    )
)

# ---------------------------------------------------------
# ROC-AUC
# ---------------------------------------------------------

roc_auc = roc_auc_score(
    y_test,
    y_prob
)

print("ROC-AUC:", round(roc_auc, 4))

# ---------------------------------------------------------
# PR-AUC
# ---------------------------------------------------------

pr_auc = average_precision_score(
    y_test,
    y_prob
)

print("PR-AUC:", round(pr_auc, 4))

# ---------------------------------------------------------
# BASELINE COMPARISON
# ---------------------------------------------------------

positive_rate = y_test.mean()

print("\n" + "=" * 70)
print("BASELINE COMPARISON")
print("=" * 70)

print(
    "Positive-class prevalence:",
    round(positive_rate, 4)
)

print(
    "Random PR-AUC baseline:",
    round(positive_rate, 4)
)

print("\n" + "=" * 70)
print("ML PHASE 4 COMPLETE")
print("=" * 70)
