import pandas as pd
import requests
import json

input_file = "ST4_annotation_input.tsv"

df = pd.read_csv(input_file, sep="\t")

test_df = df.head(10)

variants = []

for _, row in test_df.iterrows():

    variant = (
        f"{int(row['Chromosome'])} "
        f"{int(row['Position_b37'])} "
        f"{row['SNP']} "
        f"{row['Risk_Allele']} "
        f"{row['Other_Allele']} . . ."
    )

    variants.append(variant)

print("=" * 70)
print("GRCh37 VEP REGION ANNOTATION TEST")
print("=" * 70)

print("\nVariants being sent:\n")

for v in variants:
    print(v)

url = "https://grch37.rest.ensembl.org/vep/homo_sapiens/region"

headers = {
    "Content-Type": "application/json",
    "Accept": "application/json"
}

payload = {
    "variants": variants
}

response = requests.post(
    url,
    headers=headers,
    json=payload,
    timeout=180
)

print("\nHTTP status:", response.status_code)

if response.status_code != 200:
    print("\nERROR:")
    print(response.text[:3000])
    raise SystemExit(1)

results = response.json()

print("Number of returned records:", len(results))

print("\n" + "=" * 70)

for result in results:

    print("\nInput:")
    print(result.get("input"))

    print(
        "Most severe consequence:",
        result.get("most_severe_consequence")
    )

    print(
        "Transcript consequences:",
        len(result.get("transcript_consequences", []))
    )

    mappings = result.get("mappings", [])

    print(
        "Mappings returned:",
        len(mappings)
    )

    if mappings:

        mapping = mappings[0]

        print(
            "Assembly:",
            mapping.get("assembly_name")
        )

        print(
            "Chromosome:",
            mapping.get("seq_region_name")
        )

        print(
            "Position:",
            mapping.get("start")
        )

        print(
            "Allele:",
            mapping.get("allele_string")
        )

print("\n" + "=" * 70)
print("REGION ANNOTATION TEST COMPLETE")
print("=" * 70)
