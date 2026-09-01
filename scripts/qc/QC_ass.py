import pandas as pd
import numpy as np

file = "ST4_clean.tsv"

df = pd.read_csv(file, sep="\t")

print("=" * 70)
print("ST4 QC-2: ASSOCIATION STATISTICS")
print("=" * 70)

# --------------------------------------------------
# 1. P-VALUE SUMMARY
# --------------------------------------------------

p = df["MR_MEGA_Association_P"]

print("\n1. P-VALUE SUMMARY")
print(p.describe())

print("\nNumber of P-values equal to exactly 0:")
print((p == 0).sum())

print("\nNumber of P-values > 0:")
print((p > 0).sum())

print("\nSmallest non-zero P-value:")
nonzero = p[p > 0]

if len(nonzero) > 0:
    print(nonzero.min())
else:
    print("No non-zero P-values found.")

print("\nLargest P-value:")
print(p.max())

# --------------------------------------------------
# 2. UNIQUE P-VALUE EXAMPLES
# --------------------------------------------------

print("\n2. FIRST 20 UNIQUE P-VALUE VALUES")
print(p.drop_duplicates().head(20).to_string(index=False))

# --------------------------------------------------
# 3. P-VALUE DISTRIBUTION
# --------------------------------------------------

print("\n3. P-VALUE THRESHOLDS")

for threshold in [0.05, 0.01, 1e-3, 1e-5, 5e-8]:
    count = (p < threshold).sum()
    print(f"P < {threshold:g}: {count}")

# --------------------------------------------------
# 4. -LOG10(P)
# --------------------------------------------------

print("\n4. -LOG10(P) CHECK")

if (p > 0).all():
    neglogp = -np.log10(p)

    print(neglogp.describe())

    print("\nTop 10 most significant signals:")
    temp = df[["Index_SNV", "Locus", "Chromosome",
               "Position_b37", "MR_MEGA_Association_P"]].copy()

    temp["-log10P"] = neglogp

    print(
        temp.sort_values("-log10P", ascending=False)
        .head(10)
        .to_string(index=False)
    )

else:
    print(
        "WARNING: Zero P-values detected, so -log10(P) "
        "cannot be calculated directly."
    )

# --------------------------------------------------
# 5. ODDS RATIO TEXT
# --------------------------------------------------

print("\n5. OR_95CI EXAMPLES")

print(df["OR_95CI"].head(10).to_string(index=False))

# --------------------------------------------------
# 6. CHECK OR FORMAT
# --------------------------------------------------

import re

pattern = re.compile(
    r"^\s*([0-9.]+)\s*\(([0-9.]+)-([0-9.]+)\)\s*$"
)

valid_or = df["OR_95CI"].astype(str).apply(
    lambda x: bool(pattern.match(x))
)

print("\n6. OR FORMAT CHECK")
print("Valid OR_95CI format:", valid_or.sum())
print("Invalid OR_95CI format:", (~valid_or).sum())

# --------------------------------------------------
# 7. SAMPLE SIZE
# --------------------------------------------------

print("\n7. EFFECTIVE SAMPLE SIZE")

print(df["Effective_Sample_Size"].describe())

print(
    "Non-positive sample sizes:",
    (df["Effective_Sample_Size"] <= 0).sum()
)

# --------------------------------------------------

print("\n" + "=" * 70)
print("QC-2 COMPLETE")
print("=" * 70)
