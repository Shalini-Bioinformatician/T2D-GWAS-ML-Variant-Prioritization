# Integrative Genomic Variant Prioritization for Type 2 Diabetes

## Multi-Ancestry GWAS, Variant Annotation, Machine Learning, SHAP Interpretation, External GWAS Evidence, and Robust Candidate Prioritization

---

## Overview

Type 2 Diabetes (T2D) is a complex metabolic disease influenced by genetic, environmental, and lifestyle factors. Genome-wide association studies (GWAS) have identified thousands of variants associated with T2D and related metabolic traits, but identifying the variants most likely to have functional and biological relevance remains challenging.

This project develops an **integrative computational framework for prioritizing candidate genetic variants associated with Type 2 Diabetes**.

The framework combines:

- Multi-ancestry T2D GWAS evidence
- Variant-level quality control and preprocessing
- Effect-size transformation
- P-value transformation
- Minor allele frequency (MAF)
- Ensembl Variant Effect Predictor (VEP) annotation
- Protein-altering consequence prediction
- Machine learning
- Cross-validation
- Random Forest analysis
- Feature ablation analysis
- Threshold optimization
- Explainable AI using SHAP
- GWAS Catalog external evidence
- Cross-trait evidence integration
- Candidate-level biological interpretation
- Sensitivity analysis
- Rank stability analysis
- Final candidate prioritization

The objective is not simply to identify statistically significant variants, but to determine which variants show the **strongest combination of genetic association, functional annotation, machine-learning support, external GWAS evidence, and robustness**.

---

# Research Question

> **Can an integrative machine-learning and evidence-based framework prioritize genetically and biologically plausible Type 2 Diabetes candidate variants from multi-ancestry GWAS data?**

A secondary objective is to determine whether integrating independent GWAS evidence and functional annotation can identify candidates that may be missed by a machine-learning model based on a single evidence source.

---

# Objectives

### Primary Objective

Develop a reproducible computational pipeline for prioritizing candidate T2D variants using multiple complementary evidence sources.

### Specific Objectives

1. Process and quality-control T2D GWAS summary statistics.
2. Engineer informative variant-level features.
3. Annotate variants using Ensembl VEP.
4. Identify protein-altering variants and functional consequences.
5. Construct a machine-learning dataset.
6. Train and evaluate predictive models.
7. Compare alternative feature sets and models.
8. Optimize the classification threshold.
9. Apply SHAP for model interpretation.
10. Identify machine-learning-supported candidate variants.
11. Integrate external GWAS Catalog evidence.
12. Classify evidence across T2D-related traits.
13. Construct an integrated candidate evidence score.
14. Evaluate false-positive and false-negative candidates.
15. Perform sensitivity and rank-stability analysis.
16. Generate a final ranked candidate list.
17. Identify robust high-priority candidates.

---

# Dataset

## Primary GWAS Dataset

The primary dataset is derived from the multi-ancestry Type 2 Diabetes GWAS resource associated with:

**Suzuki et al., Nature, 2024**

The analysis uses the relevant supplementary GWAS information and index-signal information for downstream variant prioritization.

The processed project dataset contains **1,289 index signals across 611 unique loci** in the ST4-based analysis.

---

## External GWAS Evidence

External genetic evidence was incorporated using GWAS Catalog associations.

The external evidence was classified into categories including:

- Direct T2D associations
- Glycemic traits
- Metabolic traits
- Lipid/cardiovascular traits
- Other T2D-related evidence

This provides an independent evidence layer beyond the primary GWAS and machine-learning model.

---

# Analytical Workflow

```text
Primary T2D GWAS
       │
       ▼
Data Cleaning & QC
       │
       ▼
Effect Size Processing
       │
       ▼
P-value Transformation
       │
       ▼
MAF Calculation
       │
       ▼
Variant Feature Engineering
       │
       ▼
VEP Functional Annotation
       │
       ▼
Protein-Altering Feature Construction
       │
       ▼
Machine-Learning Dataset
       │
       ├──────────────► Baseline Model
       │
       ├──────────────► Cross-Validation
       │
       ├──────────────► Random Forest
       │
       ├──────────────► Feature Comparison
       │
       ├──────────────► Ablation Analysis
       │
       └──────────────► Threshold Optimization
       │
       ▼
Final Safe Model
       │
       ▼
SHAP Interpretation
       │
       ▼
Initial Candidate Prioritization
       │
       ▼
External GWAS Catalog Evidence
       │
       ▼
Cross-Trait Evidence Integration
       │
       ▼
Integrated Candidate Ranking
       │
       ▼
Biological Validation
       │
       ▼
Final Candidate Prioritization
       │
       ▼
Sensitivity & Rank Stability
       │
       ▼
Final Robust Candidate Set
