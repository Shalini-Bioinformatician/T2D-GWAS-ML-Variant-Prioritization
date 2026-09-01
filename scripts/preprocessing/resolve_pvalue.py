
import pandas as pd

file = "Suzuki_et_al_2024_Supplementary_tables.xlsx"

target_snps = {
    "rs7766070",
    "rs10811661",
    "rs7903146",
    "rs2237897",
    "rs1558902"
}

raw = pd.read_excel(
    file,
    sheet_name="ST4",
    header=None
)

print("=" * 80)
print("ORIGINAL ST4 P-VALUE INSPECTION")
print("=" * 80)

for snp in target_snps:

    matches = []

    for idx, row in raw.iterrows():
        if row.astype(str).str.contains(
            snp,
            case=False,
            na=False,
            regex=False
        ).any():
            matches.append(idx)

    print("\n" + "-" * 80)
    print("SNP:", snp)
    print("Matching Excel rows:", matches)

    for idx in matches:

        start = max(0, idx - 1)
        end = min(len(raw), idx + 2)

        print("\nRows around match:")

        print(
            raw.iloc[start:end].to_string(
                index=True,
                header=False
            )
        )

print("\n" + "=" * 80)
print("INSPECTION COMPLETE")
print("=" * 80)
