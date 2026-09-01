import requests
import json
import time

test_variants = [
    "rs12138597",
    "rs7766070",
    "rs10811661",
    "rs7903146",
    "rs1558902"
]

url = "https://grch37.rest.ensembl.org/vep/human/id"

headers = {
    "Content-Type": "application/json",
    "Accept": "application/json"
}

payload = {
    "ids": test_variants
}

print("Sending 5-variant annotation test to Ensembl GRCh37...")
print("Variants:", ", ".join(test_variants))

response = requests.post(
    url,
    headers=headers,
    data=json.dumps(payload),
    timeout=120
)

print("\nHTTP status:", response.status_code)

if response.status_code != 200:
    print("\nERROR:")
    print(response.text[:2000])
    raise SystemExit(1)

results = response.json()

print("\nNumber of returned records:", len(results))

for result in results:
    print("\n" + "=" * 70)
    print("Input variant:", result.get("input"))
    print("Most severe consequence:",
          result.get("most_severe_consequence"))

    mappings = result.get("mappings", [])

    if mappings:
        m = mappings[0]

        print("Assembly:", m.get("assembly_name"))
        print("Chromosome:", m.get("seq_region_name"))
        print("Start:", m.get("start"))
        print("Allele:", m.get("allele_string"))

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
print("5-VARIANT ANNOTATION TEST COMPLETE")
