import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    confusion_matrix,
    classification_report,
    roc_auc_score,
    average_precision_score
)

print("=" * 70)
print("ML PHASE 5 — RANDOM FOREST")
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

# ---------------------------------------------------------
# HANDLE MISSING VALUES
# ---------------------------------------------------------

imputer = SimpleImputer(strategy="median")

X_train = pd.DataFrame(
    imputer.fit_transform(X_train),
    columns=features,
    index=X_train.index
)

X_test = pd.DataFrame(
    imputer.transform(X_test),
    columns=features,
    index=X_test.index
)

# ---------------------------------------------------------
# RANDOM FOREST
# ---------------------------------------------------------

model = RandomForestClassifier(
    n_estimators=500,
    class_weight="balanced",
    random_state=42,
    n_jobs=-1,
    max_features="sqrt"
)

print("\nTraining Random Forest...")

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

print(
    confusion_matrix(
        y_test,
        y_pred
    )
)

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

# ---------------------------------------------------------
# PR-AUC
# ---------------------------------------------------------

pr_auc = average_precision_score(
    y_test,
    y_prob
)

print("ROC-AUC:", round(roc_auc, 4))
print("PR-AUC:", round(pr_auc, 4))

# ---------------------------------------------------------
# FEATURE IMPORTANCE
# ---------------------------------------------------------

importance = pd.DataFrame({
    "Feature": features,
    "Importance": model.feature_importances_
})

importance = importance.sort_values(
    "Importance",
    ascending=False
)

print("\n" + "=" * 70)
print("TOP 15 FEATURE IMPORTANCES")
print("=" * 70)

print(
    importance.head(15).to_string(
        index=False
    )
)

# ---------------------------------------------------------
# SAVE FEATURE IMPORTANCE
# ---------------------------------------------------------

importance.to_csv(
    "T2D_random_forest_feature_importance.tsv",
    sep="\t",
    index=False
)

print("\nFeature importance saved:")
print("T2D_random_forest_feature_importance.tsv")

print("\n" + "=" * 70)
print("ML PHASE 5 COMPLETE")
print("=" * 70)
