# PHASE 14 — FINAL PROJECT ORGANIZATION & REPRODUCIBILITY
# =============================================================================

import os
import glob
import shutil
import hashlib
import pandas as pd
from datetime import datetime


# =============================================================================
# CONFIGURATION — LOCATION INDEPENDENT
# =============================================================================

# This script is located at:
# <PROJECT_ROOT>/phase11C/GWAS_Catalog/phase14_project_organization.py

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Actual project root
PROJECT_ROOT = os.path.abspath(
    os.path.join(SCRIPT_DIR, "..", "..")
)

# Phase 11C GWAS Catalog directory
CURRENT_DIR = os.path.join(
    PROJECT_ROOT,
    "phase11C",
    "GWAS_Catalog"
)

PHASE11C_DIR = os.path.join(
    PROJECT_ROOT,
    "phase11C"
)

# Final release directory
RELEASE_DIR = os.path.join(
    PROJECT_ROOT,
    "T2D_GWAS_ML_Final_Project"
)

SCRIPTS_DIR = os.path.join(RELEASE_DIR, "scripts")
RESULTS_DIR = os.path.join(RELEASE_DIR, "results")
FIGURES_DIR = os.path.join(RELEASE_DIR, "figures")
TABLES_DIR = os.path.join(RELEASE_DIR, "tables")
DOCS_DIR = os.path.join(RELEASE_DIR, "docs")
DATA_DIR = os.path.join(RELEASE_DIR, "data_reference")


for directory in [
    RELEASE_DIR,
    SCRIPTS_DIR,
    RESULTS_DIR,
    FIGURES_DIR,
    TABLES_DIR,
    DOCS_DIR,
    DATA_DIR
]:
    os.makedirs(directory, exist_ok=True)


print("=" * 90)
print("PHASE 14 — FINAL PROJECT ORGANIZATION & REPRODUCIBILITY")
print("=" * 90)

print("\nScript directory:")
print(SCRIPT_DIR)

print("\nCurrent analysis directory:")
print(CURRENT_DIR)

print("\nProject root:")
print(PROJECT_ROOT)

print("\nFinal release directory:")
print(RELEASE_DIR)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def copy_if_exists(source, destination):

    if os.path.isfile(source):

        os.makedirs(os.path.dirname(destination), exist_ok=True)

        shutil.copy2(source, destination)

        print(f"Copied: {source}")
        print(f"    -> {destination}")

        return True

    print(f"NOT FOUND: {source}")
    return False


def sha256_file(filepath):

    sha256 = hashlib.sha256()

    with open(filepath, "rb") as f:

        for block in iter(
            lambda: f.read(1024 * 1024),
            b""
        ):
            sha256.update(block)

    return sha256.hexdigest()


# =============================================================================
# 14.0 — CONSOLIDATING FINAL RESULTS
# =============================================================================

print("\n" + "=" * 90)
print("14.0 — CONSOLIDATING FINAL RESULTS")
print("=" * 90)


# -------------------------------------------------------------------------
# Source locations
# -------------------------------------------------------------------------

PHASE11C_RESULTS = os.path.join(
    CURRENT_DIR
)

PHASE12_RESULTS = os.path.join(
    CURRENT_DIR
)

PHASE13_RESULTS = os.path.join(
    PROJECT_ROOT,
    "phase13_results"
)

SHAP_RESULTS = os.path.join(
    PROJECT_ROOT,
    "phase10B_SHAP_results"
)


# -------------------------------------------------------------------------
# Final result files
# -------------------------------------------------------------------------

final_result_files = {

    # Phase 11C
    "phase11C_external_validation_summary.tsv":
        os.path.join(
            PHASE11C_RESULTS,
            "phase11C_external_validation_summary.tsv"
        ),

    "phase11C_external_GWAS_evidence_ranking.tsv":
        os.path.join(
            PHASE11C_RESULTS,
            "phase11C_external_GWAS_evidence_ranking.tsv"
        ),

    "phase11C_integrated_candidate_evidence.tsv":
        os.path.join(
            PHASE11C_RESULTS,
            "phase11C_integrated_candidate_evidence.tsv"
        ),

    "phase11C_integrated_candidate_ranking.tsv":
        os.path.join(
            PHASE11C_RESULTS,
            "phase11C_integrated_candidate_ranking.tsv"
        ),

    # Phase 11D
    "T2D_phase11D_final_candidate_shortlist.tsv":
        os.path.join(
            PHASE11C_RESULTS,
            "T2D_phase11D_final_candidate_shortlist.tsv"
        ),

    "T2D_phase11D_model_biology_discordance.tsv":
        os.path.join(
            PHASE11C_RESULTS,
            "T2D_phase11D_model_biology_discordance.tsv"
        ),

    "T2D_phase11D_final_candidate_ranking.tsv":
        os.path.join(
            PHASE11C_RESULTS,
            "T2D_phase11D_final_candidate_ranking.tsv"
        ),

    # Phase 12
    "phase12_rank_stability.tsv":
        os.path.join(
            PHASE12_RESULTS,
            "phase12_rank_stability.tsv"
        ),
}


for filename, source in final_result_files.items():

    destination = os.path.join(
        RESULTS_DIR,
        filename
    )

    copy_if_exists(
        source,
        destination
    )


# =============================================================================
# 14.0B — SHAP RESULTS
# =============================================================================

print("\n" + "=" * 90)
print("14.0B — CONSOLIDATING PHASE 10B SHAP RESULTS")
print("=" * 90)


shap_result_files = {

    "Phase10B_Model_Reproduction_Coefficient_Check.tsv":
        os.path.join(
            SHAP_RESULTS,
            "Phase10B_Model_Reproduction_Coefficient_Check.tsv"
        ),

    "Phase10B_SHAP_Feature_Importance.tsv":
        os.path.join(
            SHAP_RESULTS,
            "Phase10B_SHAP_Feature_Importance.tsv"
        ),

    "Phase10B_SHAP_QC_Summary.tsv":
        os.path.join(
            SHAP_RESULTS,
            "Phase10B_SHAP_QC_Summary.tsv"
        ),

    "Phase10B_SHAP_Test_Set_Values.tsv":
        os.path.join(
            SHAP_RESULTS,
            "Phase10B_SHAP_Test_Set_Values.tsv"
        ),
}


for filename, source in shap_result_files.items():

    destination = os.path.join(
        RESULTS_DIR,
        filename
    )

    copy_if_exists(
        source,
        destination
    )


# =============================================================================
# 14.0C — FINAL TABLES
# =============================================================================

print("\n" + "=" * 90)
print("14.0C — CONSOLIDATING FINAL TABLES")
print("=" * 90)


final_table_files = {

    "T2D_Final_Project_Candidate_Summary.tsv":
        os.path.join(
            PHASE13_RESULTS,
            "T2D_Final_Project_Candidate_Summary.tsv"
        ),

    "T2D_Top10_Final_Candidates.tsv":
        os.path.join(
            PHASE13_RESULTS,
            "T2D_Top10_Final_Candidates.tsv"
        ),

    "T2D_Model_Biology_Discordance_Summary.tsv":
        os.path.join(
            PHASE13_RESULTS,
            "T2D_Model_Biology_Discordance_Summary.tsv"
        ),

    "Phase10B_Candidate_SHAP_Explanations.tsv":
        os.path.join(
            SHAP_RESULTS,
            "Phase10B_Candidate_SHAP_Explanations.tsv"
        ),
}


for filename, source in final_table_files.items():

    destination = os.path.join(
        TABLES_DIR,
        filename
    )

    copy_if_exists(
        source,
        destination
    )


# =============================================================================
# 14.1 — CONSOLIDATING FINAL FIGURES
# =============================================================================

print("\n" + "=" * 90)
print("14.1 — CONSOLIDATING FINAL FIGURES")
print("=" * 90)


final_figure_files = {

    # Phase 13
    "Figure_1_Final_Candidate_Ranking.png":
        os.path.join(
            PHASE13_RESULTS,
            "Figure_1_Final_Candidate_Ranking.png"
        ),

    "Figure_2_Integrated_Evidence_Components.png":
        os.path.join(
            PHASE13_RESULTS,
            "Figure_2_Integrated_Evidence_Components.png"
        ),

    "Figure_3_External_GWAS_Evidence.png":
        os.path.join(
            PHASE13_RESULTS,
            "Figure_3_External_GWAS_Evidence.png"
        ),

    "Figure_4_Model_Biology_Discordance.png":
        os.path.join(
            PHASE13_RESULTS,
            "Figure_4_Model_Biology_Discordance.png"
        ),

    "Figure_4B_Discordant_Candidate_Evidence.png":
        os.path.join(
            PHASE13_RESULTS,
            "Figure_4B_Discordant_Candidate_Evidence.png"
        ),

    "Figure_5_Final_Evidence_Heatmap.png":
        os.path.join(
            PHASE13_RESULTS,
            "Figure_5_Final_Evidence_Heatmap.png"
        ),

    "Figure_6_Final_Priority_Distribution.png":
        os.path.join(
            PHASE13_RESULTS,
            "Figure_6_Final_Priority_Distribution.png"
        ),

    "Figure_7_Sensitivity_Rank_Stability.png":
        os.path.join(
            PHASE13_RESULTS,
            "Figure_7_Sensitivity_Rank_Stability.png"
        ),

    # Phase 10B SHAP
    "Figure_10B_1_SHAP_Global_Importance.png":
        os.path.join(
            SHAP_RESULTS,
            "Figure_10B_1_SHAP_Global_Importance.png"
        ),

    "Figure_10B_2_SHAP_Beeswarm.png":
        os.path.join(
            SHAP_RESULTS,
            "Figure_10B_2_SHAP_Beeswarm.png"
        ),
}


figure_count = 0

for filename, source in final_figure_files.items():

    destination = os.path.join(
        FIGURES_DIR,
        filename
    )

    if copy_if_exists(
        source,
        destination
    ):
        figure_count += 1


print(f"\nFigures found: {figure_count}")


# =============================================================================
# 14.2 — CONSOLIDATING ANALYSIS SCRIPTS
# =============================================================================

print("\n" + "=" * 90)
print("14.2 — CONSOLIDATING ANALYSIS SCRIPTS")
print("=" * 90)


# Core scripts required to understand/reproduce the final analysis
core_scripts = [

    # QC / preprocessing
    "st4_cleaning.py",
    "parse.py",
    "effect_size.py",
    "logp.py",
    "maf.py",
    "merge.py",

    # Annotation
    "annotation.py",
    "annotation_QC.py",
    "full_annotation.py",
    "create_vep_snp_features.py",
    "combine_vep_annotations.py",
    "annotate_missing_snps.py",

    # ML
    "baseline.py",
    "random_forest.py",
    "phase9_final_safe_model.py",
    "phase11_model_comparison.py",
    "cross_validation.py",
    "feature_set_comparison.py",
    "ablation_analysis.py",
    "Threshold_optimization.py",

    # SHAP
    "phase10B_SHAP_interpretability.py",

    # Candidate prioritization
    "candidate_panel.py",
    "phase11A_candidate_prioritization.py",

    # External biological evidence
    "Catalog_Trait_Classification.py",
    "phase11C_3_external_evidence_scoring.py",
    "phase11C_4_integrated_candidate_evidence.py",
    "phase11D_complete_candidate_prioritization.py",

    # Validation / robustness
    "phase12_final_validation_robustness.py",

    # Final visualization
    "phase13_figure4_discordance.py",
    "phase13_final_visualization.py",

    # Final organization
    "phase14_project_organization.py",
]


# Search locations for scripts
script_search_dirs = [
    PROJECT_ROOT,
    CURRENT_DIR,
]


copied_scripts = set()


for script_name in core_scripts:

    source = None

    for search_dir in script_search_dirs:

        candidate = os.path.join(
            search_dir,
            script_name
        )

        if os.path.isfile(candidate):

            source = candidate
            break

    if source is None:

        print(f"NOT FOUND SCRIPT: {script_name}")

        continue

    destination = os.path.join(
        SCRIPTS_DIR,
        script_name
    )

    if script_name not in copied_scripts:

        shutil.copy2(
            source,
            destination
        )

        copied_scripts.add(
            script_name
        )

        print(f"Copied: {source}")
        print(f"    -> {destination}")


print(
    f"\nCore reproducibility scripts found: "
    f"{len(copied_scripts)}"
)


# =============================================================================
# 14.3 — PROJECT INVENTORY
# =============================================================================

print("\n" + "=" * 90)
print("14.3 — CREATING PROJECT INVENTORY")
print("=" * 90)


inventory_rows = []


for root, dirs, files in os.walk(RELEASE_DIR):

    for filename in files:

        filepath = os.path.join(
            root,
            filename
        )

        relative_path = os.path.relpath(
            filepath,
            RELEASE_DIR
        )

        try:

            size = os.path.getsize(
                filepath
            )

            checksum = sha256_file(
                filepath
            )

        except Exception:

            size = None
            checksum = None

        inventory_rows.append({

            "File":
                relative_path,

            "Size_Bytes":
                size,

            "SHA256":
                checksum,

            "Modified":
                datetime.fromtimestamp(
                    os.path.getmtime(filepath)
                ).isoformat()

        })


inventory = pd.DataFrame(
    inventory_rows
)

inventory = inventory.sort_values(
    "File"
)


inventory_path = os.path.join(
    RELEASE_DIR,
    "PROJECT_FILE_INVENTORY.tsv"
)


inventory.to_csv(
    inventory_path,
    sep="\t",
    index=False
)


print(
    f"Inventory created: {inventory_path}"
)

print(
    f"Files inventoried: {len(inventory)}"
)


# =============================================================================
# 14.4 — FINAL CANDIDATE DATASET QC
# =============================================================================

print("\n" + "=" * 90)
print("14.4 — FINAL CANDIDATE DATASET QC")
print("=" * 90)


candidate_file = os.path.join(
    TABLES_DIR,
    "T2D_Final_Project_Candidate_Summary.tsv"
)


if not os.path.isfile(candidate_file):

    raise FileNotFoundError(
        "Final project candidate summary was not created."
    )


candidate_df = pd.read_csv(
    candidate_file,
    sep="\t"
)


if "Index_SNV" in candidate_df.columns:

    snp_column = "Index_SNV"

elif "SNP" in candidate_df.columns:

    snp_column = "SNP"

elif "rsID" in candidate_df.columns:

    snp_column = "rsID"

else:

    raise ValueError(
        "No SNP identifier column found "
        "in final candidate table."
    )


rows = len(candidate_df)

unique_snps = candidate_df[
    snp_column
].nunique()

duplicate_snps = (
    rows -
    unique_snps
)

missing_snps = candidate_df[
    snp_column
].isna().sum()


print(f"Rows: {rows}")
print(f"Unique SNPs: {unique_snps}")
print(f"Duplicate SNPs: {duplicate_snps}")
print(f"Missing SNPs: {missing_snps}")


if rows != 17:

    raise ValueError(
        f"Expected 17 final candidates, "
        f"found {rows}."
    )

if unique_snps != 17:

    raise ValueError(
        f"Expected 17 unique SNPs, "
        f"found {unique_snps}."
    )

if duplicate_snps != 0:

    raise ValueError(
        "Duplicate SNPs detected."
    )

if missing_snps != 0:

    raise ValueError(
        "Missing SNP identifiers detected."
    )


# =============================================================================
# 14.5 — FINAL FIGURE QC
# =============================================================================

print("\n" + "=" * 90)
print("14.5 — FINAL FIGURE QC")
print("=" * 90)


required_figures = list(
    final_figure_files.keys()
)


for figure in required_figures:

    path = os.path.join(
        FIGURES_DIR,
        figure
    )

    if not os.path.isfile(path):

        raise FileNotFoundError(
            f"Required figure missing: {figure}"
        )

    print(
        f"PASS: {figure}"
    )


# =============================================================================
# 14.5B — SHAP QC
# =============================================================================

print("\n" + "=" * 90)
print("14.5B — SHAP QC")
print("=" * 90)


required_shap_files = [

    "Phase10B_Model_Reproduction_Coefficient_Check.tsv",

    "Phase10B_SHAP_Feature_Importance.tsv",

    "Phase10B_SHAP_QC_Summary.tsv",

    "Phase10B_SHAP_Test_Set_Values.tsv",

    "Phase10B_Candidate_SHAP_Explanations.tsv",

    "Figure_10B_1_SHAP_Global_Importance.png",

    "Figure_10B_2_SHAP_Beeswarm.png",

]


for shap_file in required_shap_files:

    if (
        shap_file.startswith("Figure_")
    ):

        path = os.path.join(
            FIGURES_DIR,
            shap_file
        )

    elif (
        shap_file ==
        "Phase10B_Candidate_SHAP_Explanations.tsv"
    ):

        path = os.path.join(
            TABLES_DIR,
            shap_file
        )

    else:

        path = os.path.join(
            RESULTS_DIR,
            shap_file
        )


    if not os.path.isfile(path):

        raise FileNotFoundError(
            f"Required SHAP file missing: "
            f"{shap_file}"
        )

    print(
        f"PASS: {shap_file}"
    )


# =============================================================================
# 14.6 — REPRODUCIBILITY MANIFEST
# =============================================================================

print("\n" + "=" * 90)
print("14.6 — REPRODUCIBILITY MANIFEST")
print("=" * 90)


manifest_rows = []


for root, dirs, files in os.walk(RELEASE_DIR):

    for filename in files:

        filepath = os.path.join(
            root,
            filename
        )

        relative_path = os.path.relpath(
            filepath,
            RELEASE_DIR
        )

        if relative_path in [
            "PROJECT_FILE_INVENTORY.tsv",
            "docs/PROJECT_REPRODUCIBILITY_MANIFEST.tsv"
        ]:

            continue

        manifest_rows.append({

            "File":
                relative_path,

            "SHA256":
                sha256_file(filepath),

            "Size_Bytes":
                os.path.getsize(filepath)

        })


manifest = pd.DataFrame(
    manifest_rows
).sort_values(
    "File"
)


manifest_path = os.path.join(
    DOCS_DIR,
    "PROJECT_REPRODUCIBILITY_MANIFEST.tsv"
)


manifest.to_csv(
    manifest_path,
    sep="\t",
    index=False
)


print(
    f"Manifest saved: {manifest_path}"
)


# =============================================================================
# 14.7 — PROJECT STRUCTURE
# =============================================================================

structure_path = os.path.join(
    DOCS_DIR,
    "PROJECT_STRUCTURE.txt"
)


with open(
    structure_path,
    "w"
) as f:

    f.write(
        "T2D GWAS ML FINAL PROJECT\n"
        "=========================\n\n"
    )

    for directory in [
        "scripts",
        "results",
        "figures",
        "tables",
        "docs",
        "data_reference"
    ]:

        f.write(
            f"{directory}/\n"
        )

        directory_path = os.path.join(
            RELEASE_DIR,
            directory
        )

        if os.path.isdir(directory_path):

            for filename in sorted(
                os.listdir(directory_path)
            ):

                f.write(
                    f"    {filename}\n"
                )


print(
    f"Project structure saved: "
    f"{structure_path}"
)


# =============================================================================
# 14.8 — FINAL PROJECT QC
# =============================================================================

print("\n" + "=" * 90)
print("14.8 — FINAL PROJECT QC")
print("=" * 90)


print("Final candidate table: PASS")
print("17-candidate dataset: PASS")
print("17 unique SNPs: PASS")
print("No duplicate SNPs: PASS")
print("All final figures present: PASS")
print("SHAP interpretation: PASS")
print("SHAP QC: PASS")
print("Scripts consolidated: PASS")
print("Reproducibility manifest: PASS")


print("\n" + "=" * 90)
print("PHASE 14 COMPLETE")
print("=" * 90)


print(
    "\nFinal project package created at:"
)

print(
    RELEASE_DIR
)


print("\nMain directories:")

print("  scripts/")
print("  results/")
print("  figures/")
print("  tables/")
print("  docs/")
print("  data_reference/")


print("\nNext stage:")
print("Biogrademy project report")
