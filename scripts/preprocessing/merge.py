
import csv

gw = "ST4_feature_engineered_v4.tsv"
vep = "ST4_vep_snp_features.tsv"

def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

gw_df = read_tsv(gw)
vep_df = read_tsv(vep)

gw_snps = {r["Index_SNV"] for r in gw_df}
vep_snps = {r["SNP"] for r in vep_df}

print("=" * 70)
print("PRE-MERGE QC")
print("=" * 70)

print("GWAS rows:", len(gw_df))
print("VEP rows:", len(vep_df))

print("GWAS unique SNPs:", len(gw_snps))
print("VEP unique SNPs:", len(vep_snps))

print("Duplicate GWAS SNPs:", len(gw_df) - len(gw_snps))
print("Duplicate VEP SNPs:", len(vep_df) - len(vep_snps))

print("SNPs in GWAS but missing from VEP:", len(gw_snps - vep_snps))
print("SNPs in VEP but missing from GWAS:", len(vep_snps - gw_snps))

print("Common SNPs:", len(gw_snps & vep_snps))

print("\nExpected:")
print("Common SNPs should be: 1289")

print("\n" + "=" * 70)
print("UTR FEATURE CHECK")
print("=" * 70)

utr5 = sum(r["Is_5UTR"] == "1" for r in vep_df)
utr3 = sum(r["Is_3UTR"] == "1" for r in vep_df)

print("SNPs with Is_5UTR = 1:", utr5)
print("SNPs with Is_3UTR = 1:", utr3)

print("\n" + "=" * 70)
print("PRE-MERGE QC COMPLETE")
print("=" * 70)
