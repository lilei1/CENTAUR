# CENTAUR Machine Learning Results Summary
==================================================

## Model Performance

| Model | AUC | Precision | Recall | F1-Score | CV AUC (Mean ± Std) |
|-------|-----|-----------|--------|----------|---------------------|
| RandomForest | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 ± 0.000 |
| LogisticRegression | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 ± 0.000 |
| NaiveBayes | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 ± 0.000 |
| SVM | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 ± 0.000 |

## Top Biomarkers by Model

### RandomForest

| Rank | Feature | Importance Score |
|------|---------|-----------------|
| 1 | fragmentomics_end_motif_CTCF | 0.1719 |
| 2 | variant_KRAS_mutation | 0.1278 |
| 3 | methylation_hypo_TMB_RNF | 0.0909 |
| 4 | methylation_hyper_TMB_PDL1 | 0.0735 |
| 5 | methylation_global_hypo_CpG | 0.0681 |

### LogisticRegression

| Rank | Feature | Importance Score |
|------|---------|-----------------|
| 1 | methylation_hypo_TMB_RNF | 0.0870 |
| 2 | variant_TMB | 0.0840 |
| 3 | variant_KRAS_mutation | 0.0800 |
| 4 | methylation_enhancer_CRXIP | 0.0745 |
| 5 | fragmentomics_entropy | 0.0719 |

### NaiveBayes

| Rank | Feature | Importance Score |
|------|---------|-----------------|
| 1 | methylation_enhancer_CRXIP | nan |
| 2 | methylation_SNP_upstream_PRK | nan |
| 3 | methylation_SNP_near_MYC | nan |
| 4 | methylation_hypo_TMB_RNF | nan |
| 5 | methylation_hyper_TMB_PDL1 | nan |

### SVM

| Rank | Feature | Importance Score |
|------|---------|-----------------|
| 1 | methylation_enhancer_CRXIP | nan |
| 2 | methylation_SNP_upstream_PRK | nan |
| 3 | methylation_SNP_near_MYC | nan |
| 4 | methylation_hypo_TMB_RNF | nan |
| 5 | methylation_hyper_TMB_PDL1 | nan |

## Best Performing Model

**RandomForest** with AUC = 1.000

## Clinical Interpretation

- **AUC > 0.8**: Good discriminative power for cancer detection
- **AUC > 0.9**: Excellent discriminative power for clinical use
- **Precision**: Proportion of predicted cancer cases that are actually cancer
- **Recall**: Proportion of actual cancer cases correctly identified
- **F1-Score**: Harmonic mean of precision and recall
