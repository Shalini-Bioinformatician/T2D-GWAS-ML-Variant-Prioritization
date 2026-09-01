import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shap

from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

warnings.filterwarnings("ignore")

print("=" * 80)
print("PHASE 10B — SHAP MODEL INTERPRETABILITY")
print("=" * 80)

# ============================================================
# PATHS
# ============================================================

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

INPUT_FILE = os.path.join(
    BASE_DIR,
    "T2D_ML_protein_altering_dataset.tsv"
)

PREDICTION_FILE = os.path.join(
    BASE_DIR,
    "T2D_final_safe_model_predictions.tsv"
)

COEFFICIENT_FILE = os.path.join(
    BASE_DIR,
    "T2D_final_safe_model_coefficients.tsv"
)

OUTPUT_DIR = os.path.join(
    BASE_DIR,
    "phase10B_SHAP_results"
)

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("\nInput dataset:")
print(INPUT_FILE)

print("\nOutput directory:")
print(OUTPUT_DIR)

# ============================================================
# LOAD DATA
# ============================================================

df = pd.read_csv(
    INPUT_FILE,
    sep="\t"
)

saved_predictions = pd.read_csv(
    PREDICTION_FILE,
    sep="\t"
)

saved_coefficients = pd.read_csv(
    COEFFICIENT_FILE,
    sep="\t"
)

print("\nDataset shape:", df.shape)
print("Saved prediction shape:", saved_predictions.shape)
print("Saved coefficient shape:", saved_coefficients.shape)

# ============================================================
# TARGET
# ============================================================

target = "Has_protein_altering"

y = df[target]

# ============================================================
# EXACT PHASE 9 FEATURES
# ============================================================

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

safe_vep_features = [
    feature
    for feature in safe_vep_features
    if feature in df.columns
    and pd.api.types.is_numeric_dtype(df[feature])
]

feature_columns = gwas_features + safe_vep_features

print("\n" + "=" * 80)
print("FEATURE RECONSTRUCTION")
print("=" * 80)

print("\nGWAS features:", len(gwas_features))
print("Numeric VEP features:", len(safe_vep_features))
print("Total features:", len(feature_columns))

for i, feature in enumerate(feature_columns, start=1):
    print(f"{i:2d}. {feature}")

# ============================================================
# VERIFY FEATURE COUNT
# ============================================================

if len(feature_columns) != 25:
    raise RuntimeError(
        f"Expected 25 Phase 9 features, found {len(feature_columns)}"
    )

if len(saved_coefficients) != 25:
    raise RuntimeError(
        f"Expected 25 saved coefficients, found {len(saved_coefficients)}"
    )

# ============================================================
# RECREATE EXACT PHASE 9 SPLIT
# ============================================================

X = df[feature_columns]

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.20,
    stratify=y,
    random_state=42
)

print("\n" + "=" * 80)
print("EXACT PHASE 9 TRAIN / TEST RECONSTRUCTION")
print("=" * 80)

print("\nTraining samples:", len(X_train))
print("Test samples:", len(X_test))

print("Training positives:", int(y_train.sum()))
print("Test positives:", int(y_test.sum()))

# ============================================================
# REBUILD EXACT PHASE 9 MODEL
# ============================================================

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

print("\nTraining reconstructed Phase 9 model...")

model.fit(
    X_train,
    y_train
)

print("Model training complete.")

# ============================================================
# REPRODUCE PREDICTIONS
# ============================================================

reconstructed_prob = model.predict_proba(
    X_test
)[:, 1]

reconstructed_pred = model.predict(
    X_test
)

# ============================================================
# VERIFY SAVED PREDICTIONS
# ============================================================

print("\n" + "=" * 80)
print("MODEL REPRODUCTION CHECK — PREDICTIONS")
print("=" * 80)

if len(saved_predictions) != len(X_test):
    raise RuntimeError(
        "Saved prediction file size does not match reconstructed test set."
    )

saved_prob = saved_predictions[
    "Protein_altering_probability"
].to_numpy()

saved_pred = saved_predictions[
    "Predicted_label"
].to_numpy()

prob_difference = np.max(
    np.abs(
        reconstructed_prob - saved_prob
    )
)

prediction_matches = np.array_equal(
    reconstructed_pred,
    saved_pred
)

print("\nMaximum probability difference:",
      f"{prob_difference:.12g}")

print("Predicted labels identical:",
      prediction_matches)

if not prediction_matches:
    raise RuntimeError(
        "Prediction reproduction FAILED."
    )

if prob_difference > 1e-10:
    raise RuntimeError(
        "Probability reproduction FAILED."
    )

print("\nPASS — Phase 9 predictions reproduced exactly.")

# ============================================================
# VERIFY COEFFICIENTS
# ============================================================

print("\n" + "=" * 80)
print("MODEL REPRODUCTION CHECK — COEFFICIENTS")
print("=" * 80)

reconstructed_coefficients = pd.DataFrame({
    "Feature": feature_columns,
    "Coefficient_Reconstructed":
        model.named_steps["classifier"].coef_[0]
})

coefficient_check = saved_coefficients[
    ["Feature", "Coefficient"]
].merge(
    reconstructed_coefficients,
    on="Feature",
    how="outer"
)

if coefficient_check.isna().any().any():
    raise RuntimeError(
        "Coefficient feature mismatch detected."
    )

coefficient_check[
    "Absolute_difference"
] = (
    coefficient_check["Coefficient"]
    - coefficient_check["Coefficient_Reconstructed"]
).abs()

max_coefficient_difference = (
    coefficient_check["Absolute_difference"].max()
)

print("\nMaximum coefficient difference:",
      f"{max_coefficient_difference:.12g}")

if max_coefficient_difference > 1e-10:
    raise RuntimeError(
        "Coefficient reproduction FAILED."
    )

print("\nPASS — Phase 9 coefficients reproduced exactly.")

coefficient_check.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "Phase10B_Model_Reproduction_Coefficient_Check.tsv"
    ),
    sep="\t",
    index=False
)

# ============================================================
# TRANSFORM DATA FOR SHAP
# ============================================================

print("\n" + "=" * 80)
print("PREPARING SHAP INPUT")
print("=" * 80)

imputer = model.named_steps["imputer"]
scaler = model.named_steps["scaler"]

X_train_imputed = imputer.transform(X_train)
X_test_imputed = imputer.transform(X_test)

X_train_scaled = scaler.transform(
    X_train_imputed
)

X_test_scaled = scaler.transform(
    X_test_imputed
)

X_train_scaled = pd.DataFrame(
    X_train_scaled,
    columns=feature_columns,
    index=X_train.index
)

X_test_scaled = pd.DataFrame(
    X_test_scaled,
    columns=feature_columns,
    index=X_test.index
)

# ============================================================
# SHAP EXPLAINER
# ============================================================

print("\nCreating SHAP LinearExplainer...")

classifier = model.named_steps["classifier"]

explainer = shap.LinearExplainer(
    classifier,
    X_train_scaled
)

shap_values = explainer(
    X_test_scaled
)

print("SHAP calculation complete.")

# ============================================================
# HANDLE SHAP OUTPUT
# ============================================================

values = shap_values.values

if values.ndim == 3:
    values = values[:, :, 1]

print("\nSHAP matrix shape:", values.shape)

if values.shape != X_test_scaled.shape:
    raise RuntimeError(
        "SHAP matrix dimensions do not match test feature matrix."
    )

# ============================================================
# GLOBAL SHAP IMPORTANCE
# ============================================================

mean_abs_shap = np.abs(values).mean(axis=0)
mean_shap = values.mean(axis=0)

shap_importance = pd.DataFrame({
    "Feature": feature_columns,
    "Mean_Absolute_SHAP": mean_abs_shap,
    "Mean_SHAP": mean_shap
})

shap_importance = shap_importance.sort_values(
    "Mean_Absolute_SHAP",
    ascending=False
).reset_index(drop=True)

shap_importance.insert(
    0,
    "SHAP_Rank",
    np.arange(1, len(shap_importance) + 1)
)

shap_importance.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "Phase10B_SHAP_Feature_Importance.tsv"
    ),
    sep="\t",
    index=False
)

print("\nTop SHAP features:")

print(
    shap_importance.head(15).to_string(
        index=False
    )
)

# ============================================================
# SHAP VALUES TABLE
# ============================================================

shap_values_df = pd.DataFrame(
    values,
    columns=feature_columns,
    index=X_test.index
)

shap_values_df.insert(
    0,
    "Index_SNV",
    df.loc[
        X_test.index,
        "Index_SNV"
    ].values
)

shap_values_df.insert(
    1,
    "True_label",
    y_test.values
)

shap_values_df.insert(
    2,
    "Predicted_probability",
    reconstructed_prob
)

shap_values_df.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "Phase10B_SHAP_Test_Set_Values.tsv"
    ),
    sep="\t",
    index=False
)

# ============================================================
# GLOBAL SHAP BAR PLOT
# ============================================================

print("\nCreating SHAP global importance plot...")

plt.figure(
    figsize=(10, 8)
)

shap.summary_plot(
    values,
    X_test_scaled,
    plot_type="bar",
    show=False
)

plt.title(
    "SHAP Global Feature Importance — Phase 9 Logistic Regression"
)

plt.tight_layout()

plt.savefig(
    os.path.join(
        OUTPUT_DIR,
        "Figure_10B_1_SHAP_Global_Importance.png"
    ),
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# ============================================================
# SHAP BEESWARM
# ============================================================

print("Creating SHAP beeswarm plot...")

plt.figure(
    figsize=(10, 9)
)

shap.summary_plot(
    values,
    X_test_scaled,
    show=False
)

plt.title(
    "SHAP Feature Contributions — Phase 9 Logistic Regression"
)

plt.tight_layout()

plt.savefig(
    os.path.join(
        OUTPUT_DIR,
        "Figure_10B_2_SHAP_Beeswarm.png"
    ),
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# ============================================================
# CANDIDATE-LEVEL SHAP
# ============================================================

print("\nCreating candidate-level SHAP explanations...")

candidate_columns = [
    "Index_SNV",
    "Primary_gene"
]

available_candidate_columns = [
    c for c in candidate_columns
    if c in df.columns
]

candidate_info = df.loc[
    X_test.index,
    available_candidate_columns
].copy()

candidate_info["Predicted_probability"] = reconstructed_prob

candidate_info["True_label"] = y_test.values

candidate_shap_rows = []

for i, idx in enumerate(X_test.index):

    row_values = values[i]

    top_indices = np.argsort(
        np.abs(row_values)
    )[::-1][:10]

    for rank, feature_index in enumerate(
        top_indices,
        start=1
    ):

        candidate_shap_rows.append({
            "Index_SNV":
                df.loc[idx, "Index_SNV"],
            "Primary_gene":
                df.loc[idx, "Primary_gene"]
                if "Primary_gene" in df.columns
                else "",
            "True_label":
                y_test.loc[idx],
            "Predicted_probability":
                reconstructed_prob[i],
            "SHAP_Rank":
                rank,
            "Feature":
                feature_columns[feature_index],
            "SHAP_value":
                row_values[feature_index],
            "Absolute_SHAP":
                abs(row_values[feature_index])
        })

candidate_shap = pd.DataFrame(
    candidate_shap_rows
)

candidate_shap.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "Phase10B_Candidate_SHAP_Explanations.tsv"
    ),
    sep="\t",
    index=False
)

# ============================================================
# REPRODUCIBILITY SUMMARY
# ============================================================

summary = pd.DataFrame({
    "Metric": [
        "Dataset_rows",
        "Dataset_columns",
        "Number_of_features",
        "Training_samples",
        "Test_samples",
        "SHAP_samples",
        "Prediction_reproduction",
        "Coefficient_reproduction",
        "Maximum_probability_difference",
        "Maximum_coefficient_difference",
        "SHAP_version"
    ],
    "Value": [
        len(df),
        len(df.columns),
        len(feature_columns),
        len(X_train),
        len(X_test),
        len(values),
        "PASS",
        "PASS",
        prob_difference,
        max_coefficient_difference,
        shap.__version__
    ]
})

summary.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "Phase10B_SHAP_QC_Summary.tsv"
    ),
    sep="\t",
    index=False
)

# ============================================================
# COMPLETION
# ============================================================

print("\n" + "=" * 80)
print("PHASE 10B COMPLETE")
print("=" * 80)

print("\nOutputs saved to:")
print(OUTPUT_DIR)

print("\nKey outputs:")

for filename in [
    "Phase10B_Model_Reproduction_Coefficient_Check.tsv",
    "Phase10B_SHAP_Feature_Importance.tsv",
    "Phase10B_SHAP_Test_Set_Values.tsv",
    "Phase10B_Candidate_SHAP_Explanations.tsv",
    "Phase10B_SHAP_QC_Summary.tsv",
    "Figure_10B_1_SHAP_Global_Importance.png",
    "Figure_10B_2_SHAP_Beeswarm.png"
]:
    print(" -", filename)

print("\nModel reproduction:")
print("Prediction check: PASS")
print("Coefficient check: PASS")

print("\nSHAP interpretation generated successfully.")
