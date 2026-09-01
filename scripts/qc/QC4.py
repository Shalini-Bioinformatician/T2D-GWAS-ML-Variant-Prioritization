
import pandas as pd

file = "Suzuki_et_al_2024_Supplementary_tables.xlsx"

# Read ST4 as raw values
raw = pd.read_excel(
    file,
    sheet_name="ST4",
    header=None,
    dtype=object
)

# Search for the five SNPs
target_snvs = [
    "rs7766070",
    "rs10811661",
    "rs7903146",
    "rs2237897",
    "rs1558902"
]

for snv in target_snvs:
    matches = raw[
        raw.apply(
            lambda row: row.astype(str).str.strip().eq(snv).any(),
            axis=1
        )
    ]

    print("\n" + "="*80)
    print("SNV:", snv)
    print("="*80)

    if len(matches) == 0:
        print("NOT FOUND")
    else:
        print(matches.to_string(index=True, header=False))