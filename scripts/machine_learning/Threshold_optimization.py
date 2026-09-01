import pandas as pd
import numpy as np

from sklearn.metrics import (
    precision_score,
    recall_score,
    f1_score,
    confusion_matrix
)

print("=" * 70)
print("ML PHASE 10B — THRESHOLD OPTIMIZATION")
print("=" * 70)

# ---------------------------------------------------------
# LOAD PREDICTIONS
# ---------------------------------------------------------

df = pd.read_csv(
    "T2D_final_safe_model_predictions.tsv",
    sep="\t"
)

print("\nPrediction file shape:")
print(df.shape)

# ---------------------------------------------------------
# DEFINE TRUE LABELS AND PROBABILITIES
# ---------------------------------------------------------

y_true = df["True_label"]

y_prob = df["Protein_altering_probability"]

print("\nPositive SNPs in test set:")
print(y_true.sum())

# ---------------------------------------------------------
# TEST DIFFERENT THRESHOLDS
# ---------------------------------------------------------

thresholds = [
    0.05,
    0.10,
    0.15,
    0.20,
    0.25,
    0.30,
    0.35,
    0.40,
    0.45,
    0.50,
    0.55,
    0.60,
    0.65,
    0.70,
    0.75,
    0.80,
    0.85,
    0.90
]

results = []

print("\n" + "-" * 70)
print("THRESHOLD PERFORMANCE")
print("-" * 70)

for threshold in thresholds:

    y_pred = (y_prob >= threshold).astype(int)

    precision = precision_score(
        y_true,
        y_pred,
        zero_division=0
    )

    recall = recall_score(
        y_true,
        y_pred,
        zero_division=0
    )

    f1 = f1_score(
        y_true,
        y_pred,
        zero_division=0
    )

    tn, fp, fn, tp = confusion_matrix(
        y_true,
        y_pred
    ).ravel()

    results.append({
        "Threshold": threshold,
        "Precision": precision,
        "Recall": recall,
        "F1_score": f1,
        "True_positive": tp,
        "False_positive": fp,
        "True_negative": tn,
        "False_negative": fn
    })

results_df = pd.DataFrame(results)

print(
    results_df.to_string(
        index=False
    )
)

# ---------------------------------------------------------
# BEST THRESHOLD BY F1 SCORE
# ---------------------------------------------------------

best_row = results_df.loc[
    results_df["F1_score"].idxmax()
]

print("\n" + "=" * 70)
print("BEST THRESHOLD — BASED ON F1 SCORE")
print("=" * 70)

print(best_row)

# ---------------------------------------------------------
# SAVE RESULTS
# ---------------------------------------------------------

results_df.to_csv(
    "T2D_threshold_optimization.tsv",
    sep="\t",
    index=False
)

print("\nSaved:")
print("T2D_threshold_optimization.tsv")

print("\n" + "=" * 70)
print("ML PHASE 10B COMPLETE")
print("=" * 70)