import requests
import json
import time

test_variants = [
    "rs10088883",
    "rs10096860",
    "rs10113902",
    "rs10116566",
    "rs10121760",
    "rs10147577",
    "rs10152999",
    "rs1022736",
    "rs1029016",
    "rs1044531"
]

url = "https://grch37.rest.ensembl.org/vep/human/id"

headers = {
    "Content-Type": "application/json",
    "Accept": "application/json"
}

payload = {
    "ids": test_variants
}

print("=" * 70)
print("TESTING MISSING SNPs DIRECTLY AGAINST ENSEMBL GRCh37")
print("=" * 70)

print("\nVariants being tested:")
print(", ".join(test_variants))

response = requests.post(
    url,
    headers=headers,
    data=json.dumps(payload),
    timeout=120
)

print("\nHTTP status:", response.status_code)

if response.status_code != 200:
    print("\nERROR:")
    print(response.text[:3000])
    raise SystemExit(1)

results = response.json()

print("Number of returned records:", len(results))

returned = set()

for result in results:

    snp = result.get("input")
    returned.add(snp)

    print("\n" + "-" * 70)
    print("Input SNP:", snp)

    print(
        "Most severe consequence:",
        result.get("most_severe_consequence")
    )

    mappings = result.get("mappings", [])

    if mappings:
        m = mappings[0]

        print("Assembly:", m.get("assembly_name"))
        print("Chromosome:", m.get("seq_region_name"))
        print("Position:", m.get("start"))
        print("Allele:", m.get("allele_string"))
    else:
        print("No mapping information returned")

    transcripts = result.get("transcript_consequences", [])

    print("Transcript consequences:", len(transcripts))

    for tc in transcripts[:3]:

        print(
            "  Gene:",
            tc.get("gene"),
            "| Symbol:",
            tc.get("gene_symbol"),
            "| Consequence:",
            ",".join(tc.get("consequence_terms", [])),
            "| Impact:",
            tc.get("impact")
        )

print("\n" + "=" * 70)
print("RETURN SUMMARY")
print("=" * 70)

print("Requested:", len(test_variants))
print("Returned:", len(returned))
print("Not returned:", len(test_variants) - len(returned))

not_returned = set(test_variants) - returned

if not_returned:
    print("\nSNPs NOT returned by Ensembl:")
    for snp in sorted(not_returned):
        print(snp)

print("\n" + "=" * 70)
print("TEST COMPLETE")
print("=" * 70)

