
import pandas as pd
import requests
import time
import json
import os

input_file = "ST4_annotation_input.tsv"
output_file = "ST4_vep_GRCh37_annotations.tsv"

df = pd.read_csv(input_file, sep="\t")

variants = df["SNP"].tolist()

url = "https://grch37.rest.ensembl.org/vep/human/id"

headers = {
    "Content-Type": "application/json",
    "Accept": "application/json"
}

batch_size = 100

results = []

total = len(variants)

print("=" * 70)
print("FULL GRCh37 VARIANT ANNOTATION")
print("=" * 70)
print("Total variants:", total)
print("Batch size:", batch_size)
print("Number of batches:", (total + batch_size - 1) // batch_size)
print()

for start in range(0, total, batch_size):

    batch = variants[start:start + batch_size]

    print(
        f"Processing variants {start + 1}-{start + len(batch)} "
        f"of {total}..."
    )

    payload = {
        "ids": batch
    }

    success = False

    for attempt in range(3):

        try:

            response = requests.post(
                url,
                headers=headers,
                data=json.dumps(payload),
                timeout=180
            )

            if response.status_code == 200:

                batch_results = response.json()
                results.extend(batch_results)

                print(
                    f"  SUCCESS: {len(batch_results)} annotations returned"
                )

                success = True
                break

            else:

                print(
                    f"  Attempt {attempt + 1}: "
                    f"HTTP {response.status_code}"
                )

                time.sleep(5)

        except Exception as e:

            print(
                f"  Attempt {attempt + 1} failed: {e}"
            )

            time.sleep(5)

    if not success:

        print(
            f"  FAILED batch starting at variant {start + 1}"
        )

    time.sleep(1)

print()
print("=" * 70)
print("ANNOTATION REQUESTS COMPLETE")
print("=" * 70)

print("Input variants:", total)
print("Returned records:", len(results))

# ------------------------------------------------------------------
# Extract useful annotation fields
# ------------------------------------------------------------------

annotation_rows = []

for result in results:

    snp = result.get("input")

    most_severe = result.get(
        "most_severe_consequence"
    )

    mappings = result.get("mappings", [])

    if mappings:

        mapping = mappings[0]

        assembly = mapping.get("assembly_name")
        chromosome = mapping.get("seq_region_name")
        position = mapping.get("start")
        allele_string = mapping.get("allele_string")

    else:

        assembly = None
        chromosome = None
        position = None
        allele_string = None

    transcripts = result.get(
        "transcript_consequences",
        []
    )

    if transcripts:

        for tc in transcripts:

            annotation_rows.append({

                "SNP": snp,

                "Assembly": assembly,

                "Annotation_Chromosome": chromosome,

                "Annotation_Position": position,

                "Allele_string": allele_string,

                "Most_severe_consequence":
                    most_severe,

                "Gene_symbol":
                    tc.get("gene_symbol"),

                "Gene_ID":
                    tc.get("gene"),

                "Transcript_ID":
                    tc.get("transcript_id"),

                "Consequence":
                    ",".join(
                        tc.get(
                            "consequence_terms",
                            []
                        )
                    ),

                "Impact":
                    tc.get("impact"),

                "Biotype":
                    tc.get("biotype"),

                "Exon":
                    tc.get("exon"),

                "Protein_position":
                    tc.get("protein_start"),

                "Amino_acids":
                    tc.get("amino_acids")

            })

    else:

        annotation_rows.append({

            "SNP": snp,

            "Assembly": assembly,

            "Annotation_Chromosome": chromosome,

            "Annotation_Position": position,

            "Allele_string": allele_string,

            "Most_severe_consequence":
                most_severe,

            "Gene_symbol": None,

            "Gene_ID": None,

            "Transcript_ID": None,

            "Consequence": most_severe,

            "Impact": None,

            "Biotype": None,

            "Exon": None,

            "Protein_position": None,

            "Amino_acids": None

        })

annotation_df = pd.DataFrame(annotation_rows)

annotation_df.to_csv(
    output_file,
    sep="\t",
    index=False
)

print()
print("Output file:")
print(output_file)

print()
print("Annotation rows:", len(annotation_df))
print("Unique SNPs annotated:",
      annotation_df["SNP"].nunique())

print()
print("Most severe consequence distribution:")

print(
    annotation_df[
        [
            "SNP",
            "Most_severe_consequence"
        ]
    ]
    .drop_duplicates()
    ["Most_severe_consequence"]
    .value_counts()
)

print()
print("Annotation complete.")