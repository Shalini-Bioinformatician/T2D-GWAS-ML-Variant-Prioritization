import pandas as pd
import requests
import time
import json

input_file = "ST4_missing_annotation_snps.tsv"
successful_file = "ST4_missing_vep_annotations.tsv"

df = pd.read_csv(
    input_file,
    sep="\t",
    dtype=str
)

successful = pd.read_csv(
    successful_file,
    sep="\t",
    dtype=str
)

requested_snps = set(df["SNP"])

successful_snps = set(successful["SNP"])

remaining_snps = sorted(
    requested_snps - successful_snps
)

print("=" * 70)
print("RETRYING UNANNOTATED SNPs")
print("=" * 70)

print("Original missing SNPs:", len(requested_snps))
print("Already annotated:", len(successful_snps))
print("Remaining SNPs:", len(remaining_snps))

url = "https://grch37.rest.ensembl.org/vep/human/id"

headers = {
    "Content-Type": "application/json",
    "Accept": "application/json"
}


def query_batch(batch):

    payload = {
        "ids": batch
    }

    try:

        response = requests.post(
            url,
            headers=headers,
            json=payload,
            timeout=180
        )

        print(
            f"Batch size {len(batch)} → "
            f"HTTP {response.status_code}"
        )

        if response.status_code == 200:

            return response.json()

        return None

    except Exception as e:

        print("Request error:", e)

        return None


# --------------------------------------------------
# Recursive batch splitting
# --------------------------------------------------

def annotate_recursive(batch):

    if not batch:
        return []

    print()
    print(
        "Trying batch:",
        len(batch),
        "SNPs"
    )

    results = query_batch(batch)

    if results is not None:

        print(
            "SUCCESS:",
            len(results),
            "records returned"
        )

        return results

    # If one SNP fails, isolate it
    if len(batch) == 1:

        print(
            "FAILED INDIVIDUAL SNP:",
            batch[0]
        )

        return []

    midpoint = len(batch) // 2

    first_half = batch[:midpoint]
    second_half = batch[midpoint:]

    print(
        "Splitting into:",
        len(first_half),
        "+",
        len(second_half)
    )

    time.sleep(2)

    results = []

    results.extend(
        annotate_recursive(first_half)
    )

    time.sleep(2)

    results.extend(
        annotate_recursive(second_half)
    )

    return results


# --------------------------------------------------
# Run
# --------------------------------------------------

results = annotate_recursive(
    remaining_snps
)


print()
print("=" * 70)
print("RETRY COMPLETE")
print("=" * 70)

print(
    "Remaining SNPs before retry:",
    len(remaining_snps)
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


output_file = (
    "ST4_remaining_vep_annotations.tsv"
)

annotation_df = pd.DataFrame(
    annotation_rows
)

annotation_df.to_csv(
    output_file,
    sep="\t",
    index=False
)

print()
print("Output file:")
print(output_file)

print(
    "Annotation rows:",
    len(annotation_df)
)

print(
    "Unique SNPs recovered:",
    annotation_df["SNP"].nunique()
)
