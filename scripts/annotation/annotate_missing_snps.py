import pandas as pd
import requests
import time
import json

input_file = "ST4_missing_annotation_snps.tsv"
output_file = "ST4_missing_vep_annotations.tsv"

df = pd.read_csv(input_file, sep="\t", dtype=str)

variants = df["SNP"].tolist()

url = "https://grch37.rest.ensembl.org/vep/human/id"

headers = {
    "Content-Type": "application/json",
    "Accept": "application/json"
}

batch_size = 100

results = []

print("=" * 70)
print("ANNOTATING MISSING SNPs")
print("=" * 70)

print("Missing SNPs:", len(variants))
print("Batch size:", batch_size)
print()

for start in range(0, len(variants), batch_size):

    batch = variants[start:start + batch_size]

    print(
        f"Processing {start + 1}-{start + len(batch)} "
        f"of {len(variants)}..."
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
                json=payload,
                timeout=180
            )

            print(
                "  HTTP status:",
                response.status_code
            )

            if response.status_code == 200:

                batch_results = response.json()

                results.extend(batch_results)

                print(
                    "  Returned:",
                    len(batch_results)
                )

                success = True
                break

            else:

                print(
                    "  Attempt",
                    attempt + 1,
                    "failed"
                )

                print(
                    response.text[:500]
                )

                time.sleep(5)

        except Exception as e:

            print(
                "  Attempt",
                attempt + 1,
                "error:",
                e
            )

            time.sleep(5)

    if not success:

        print(
            "  WARNING: batch failed"
        )

    time.sleep(1)


print()
print("=" * 70)
print("REQUESTS COMPLETE")
print("=" * 70)

print(
    "Requested SNPs:",
    len(variants)
)

print(
    "Returned records:",
    len(results)
)


# --------------------------------------------------
# Extract annotations
# --------------------------------------------------

annotation_rows = []

for result in results:

    snp = result.get("input")

    most_severe = result.get(
        "most_severe_consequence"
    )

    transcripts = result.get(
        "transcript_consequences",
        []
    )

    if transcripts:

        for tc in transcripts:

            annotation_rows.append({

                "SNP": snp,

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

            "Most_severe_consequence":
                most_severe,

            "Gene_symbol": None,

            "Gene_ID": None,

            "Transcript_ID": None,

            "Consequence":
                most_severe,

            "Impact": None,

            "Biotype": None,

            "Exon": None,

            "Protein_position": None,

            "Amino_acids": None

        })


annotation_df = pd.DataFrame(
    annotation_rows
)

annotation_df.to_csv(
    output_file,
    sep="\t",
    index=False
)


print()
print("=" * 70)
print("MISSING SNP ANNOTATION COMPLETE")
print("=" * 70)

print(
    "Annotation rows:",
    len(annotation_df)
)

print(
    "Unique SNPs annotated:",
    annotation_df["SNP"].nunique()
)

print()
print("Output:")
print(output_file)
