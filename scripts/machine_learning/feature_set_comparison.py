import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score
)

print("=" * 70)
print("ML PHASE 6 — GWAS-ONLY vs GWAS + VEP")
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

y = df[target]

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
# VEP FEATURES
# ---------------------------------------------------------

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

# ---------------------------------------------------------
# COMBINED FEATURES
# ---------------------------------------------------------

combined_features = gwas_features + vep_features

print("\nGWAS features:", len(gwas_features))
print("VEP features:", len(vep_features))
print("Combined features:", len(combined_features))

# ---------------------------------------------------------
# SAME STRATIFIED SPLIT FOR BOTH MODELS
# ---------------------------------------------------------

train_idx, test_idx = train_test_split(
    df.index,
    test_size=0.20,
    random_state=42,
    stratify=y
)

y_train = y.loc[train_idx]
y_test = y.loc[test_idx]

# ---------------------------------------------------------
# FUNCTION TO TRAIN + EVALUATE
# ---------------------------------------------------------

def evaluate_model(name, feature_list):

    X = df[feature_list].copy()

    X_train = X.loc[train_idx]
    X_test = X.loc[test_idx]

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

    model.fit(
        X_train,
        y_train
    )

    probability = model.predict_proba(
        X_test
    )[:, 1]

    roc_auc = roc_auc_score(
        y_test,
        probability
    )

    pr_auc = average_precision_score(
        y_test,
        probability
    )

    print("\n" + "-" * 70)
    print(name)
    print("-" * 70)

    print("Features:", len(feature_list))
    print("ROC-AUC:", round(roc_auc, 4))
    print("PR-AUC:", round(pr_auc, 4))

    return roc_auc, pr_auc


# ---------------------------------------------------------
# MODEL A
# ---------------------------------------------------------

gwas_roc, gwas_pr = evaluate_model(
    "MODEL A — GWAS FEATURES ONLY",
    gwas_features
)

# ---------------------------------------------------------
# MODEL B
# ---------------------------------------------------------

combined_roc, combined_pr = evaluate_model(
    "MODEL B — GWAS + VEP FEATURES",
    combined_features
)

# ---------------------------------------------------------
# IMPROVEMENT
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("PERFORMANCE COMPARISON")
print("=" * 70)

print(
    "ROC-AUC improvement:",
    round(combined_roc - gwas_roc, 4)
)

print(
    "PR-AUC improvement:",
    round(combined_pr - gwas_pr, 4)
)

print("\n" + "=" * 70)
print("ML PHASE 6 COMPLETE")
print("=" * 70)
