# T2D-GWAS-ML-Variant-Prioritization
Multi-ancestry Type 2 Diabetes GWAS analysis with variant annotation, machine learning, SHAP-based interpretation, and candidate variant prioritization.
# Integrative Genomic Variant Prioritization for Type 2 Diabetes

## Overview

## Research Question

## Objectives

## Dataset

## Workflow

## Methods

## Machine Learning

## Explainable AI / SHAP

## Candidate Prioritization

## Key Results

## Repository Structure

## Reproducibility

## Limitations

## Citation
## 📁 Repository Structure

```text
T2D-GWAS-ML-Variant-Prioritization/
│
├── README.md
│
├── data_reference/
│   └── README.md
│
├── scripts/
│   ├── st4_cleaning.py
│   ├── parse.py
│   ├── effect_size.py
│   ├── logp.py
│   ├── maf.py
│   ├── annotation.py
│   ├── full_annotation.py
│   ├── phase9_final_safe_model.py
│   ├── phase10B_SHAP_interpretability.py
│   ├── phase11A_candidate_prioritization.py
│   ├── phase11C_3_external_evidence_scoring.py
│   ├── phase11C_4_integrated_candidate_evidence.py
│   ├── phase11D_complete_candidate_prioritization.py
│   ├── phase12_final_validation_robustness.py
│   ├── phase13_final_visualization.py
|   ├── phase13_figure4_discordance.py
│   └── ...
│
├── results/
│   ├── phase11C_external_validation_summary.tsv
│   ├── phase11C_external_GWAS_evidence_ranking.tsv
│   ├── phase11C_integrated_candidate_evidence.tsv
│   ├── phase11C_integrated_candidate_ranking.tsv
│   ├── T2D_phase11D_final_candidate_ranking.tsv
│   ├── T2D_phase11D_final_candidate_shortlist.tsv
|   ├── T2D_phase11D_model_biology_discordance.tsv
│   ├── phase12_rank_stability.tsv
│   ├── Phase10B_Model_Reproduction_Coefficient_Check.tsv
│   ├── Phase10B_SHAP_Feature_Importance.tsv
│   ├── Phase10B_SHAP_QC_Summary.tsv
│   ├── Phase10B_SHAP_Test_Set_Values.tsv
│   ├── phase14_figure_qc.tsv
│
├── tables/
│   └── Phase10B_Candidate_SHAP_Explanations.tsv
│
├── figures/
│   ├── Figure_1_Final_Candidate_Ranking.png
│   ├── Figure_2_Integrated_Evidence_Components.png
│   ├── Figure_3_External_GWAS_Evidence.png
│   ├── Figure_4_Model_Biology_Discordance.png
│   ├── Figure_4B_Discordant_Candidate_Evidence.png
│   ├── Figure_5_Final_Evidence_Heatmap.png
│   ├── Figure_6_Final_Priority_Distribution.png
│   ├── Figure_7_Sensitivity_Rank_Stability.png
│   ├── Figure_10B_1_SHAP_Global_Importance.png
│   └── Figure_10B_2_SHAP_Beeswarm.png
│
└── docs/
    ├── PROJECT_STRUCTURE.txt
    ├── PROJECT_REPRODUCIBILITY_MANIFEST.tsv
    └── ...
