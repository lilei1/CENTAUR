#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Biomarker discovery process
process BIOMARKER_DISCOVERY {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type)
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_biomarker_discovery.txt"), emit: discovery_log
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_statistical_tests.tsv"), emit: statistical_tests
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_feature_importance.tsv"), emit: feature_importance
    
    script:
    """
    echo "Starting biomarker discovery for ${sample_id}" > ${sample_id}_biomarker_discovery.txt
    echo "Sample type: ${sample_type}" >> ${sample_id}_biomarker_discovery.txt
    echo "Cancer type: ${cancer_type}" >> ${sample_id}_biomarker_discovery.txt
    echo "Stage: ${stage}" >> ${sample_id}_biomarker_discovery.txt
    echo "Processing completed at: \$(date)" >> ${sample_id}_biomarker_discovery.txt
    
    # Create Python script for biomarker discovery
    python3 -c "
import pandas as pd
import numpy as np
import random

print(f'Starting biomarker discovery for ${sample_id}')

# Simulate statistical testing results
def simulate_statistical_test(feature_name, sample_type):
    if sample_type == 'Cancer':
        effect_size = random.uniform(0.3, 0.8)
        p_value = random.uniform(0.001, 0.05)
        fdr = p_value * random.uniform(1.0, 2.0)
    else:
        effect_size = random.uniform(0.1, 0.4)
        p_value = random.uniform(0.01, 0.1)
        fdr = p_value * random.uniform(1.5, 3.0)
    
    if 'TMB' in feature_name or 'immune' in feature_name.lower():
        auc = random.uniform(0.8, 0.95)
    elif 'methylation' in feature_name.lower():
        auc = random.uniform(0.7, 0.9)
    else:
        auc = random.uniform(0.6, 0.85)
    
    return {
        'feature': feature_name,
        'effect_size': effect_size,
        'p_value': p_value,
        'fdr': fdr,
        'auc': auc,
        'significant': fdr < 0.05,
        'biologically_relevant': auc > 0.7
    }

# Perform statistical testing on all features
statistical_results = []

all_features = [
    'methylation_global', 'methylation_immune_enhancers', 
    'methylation_cancer_drivers', 'methylation_super_enhancers',
    'fragmentomics_short_ratio', 'fragmentomics_entropy',
    'fragmentomics_wps_score', 'fragmentomics_end_motif_bias',
    'expression_ifn_gamma', 'expression_t_cell_exhaustion',
    'expression_cytotoxic_t_cells', 'expression_b_cells',
    'variant_tmb', 'variant_total_mutations', 'variant_snv_count',
    'variant_indel_count', 'variant_cnv_gains', 'variant_cnv_losses'
]

for feature in all_features:
    result = simulate_statistical_test(feature, '${sample_type}')
    statistical_results.append(result)

# Save statistical test results
stats_df = pd.DataFrame(statistical_results)
stats_df.to_csv('${sample_id}_statistical_tests.tsv', sep='\t', index=False)

# Calculate feature importance
feature_importance = []
for result in statistical_results:
    importance = result['auc'] * (1 - result['fdr']) * result['effect_size']
    feature_importance.append({
        'feature': result['feature'],
        'importance_score': importance,
        'auc': result['auc'],
        'fdr': result['fdr'],
        'effect_size': result['effect_size'],
        'rank': 0
    })

# Rank features by importance
feature_importance_df = pd.DataFrame(feature_importance)
feature_importance_df = feature_importance_df.sort_values('importance_score', ascending=False)
feature_importance_df['rank'] = range(1, len(feature_importance_df) + 1)
feature_importance_df.to_csv('${sample_id}_feature_importance.tsv', sep='\t', index=False)

print(f'Biomarker discovery completed for ${sample_id}')
print(f'Total features tested: {len(statistical_results)}')
print(f'Significant features: {sum(1 for r in statistical_results if r[\"significant\"])}')
print(f'High-AUC features: {sum(1 for r in statistical_results if r[\"auc\"] > 0.8)}')
"
    """
}

// Machine learning training process
process ML_TRAINING {
    tag "${model_type}"
    
    input:
    val(model_type)
    
    output:
    tuple val(model_type), path("${model_type}_model_results.txt"), emit: model_log
    tuple val(model_type), path("${model_type}_performance_metrics.tsv"), emit: performance_metrics
    tuple val(model_type), path("${model_type}_feature_importance.tsv"), emit: feature_importance
    tuple val(model_type), path("${model_type}_predictions.tsv"), emit: predictions
    
    script:
    """
    echo "Training ${model_type} model" > ${model_type}_model_results.txt
    echo "Training completed at: \$(date)" >> ${model_type}_model_results.txt
    
    # Create Python script for ML training
    python3 -c "
import pandas as pd
import numpy as np
import random

print(f'Training ${model_type} model')

# Simulate ML training results
def simulate_ml_training(model_type, n_features):
    if model_type == 'XGBoost':
        auc = random.uniform(0.85, 0.95)
        precision = random.uniform(0.80, 0.90)
        recall = random.uniform(0.75, 0.85)
        f1 = random.uniform(0.78, 0.88)
    elif model_type == 'RandomForest':
        auc = random.uniform(0.80, 0.90)
        precision = random.uniform(0.75, 0.85)
        recall = random.uniform(0.70, 0.80)
        f1 = random.uniform(0.72, 0.82)
    elif model_type == 'LogisticRegression':
        auc = random.uniform(0.70, 0.85)
        precision = random.uniform(0.65, 0.80)
        recall = random.uniform(0.60, 0.75)
        f1 = random.uniform(0.62, 0.77)
    else:  # NaiveBayes
        auc = random.uniform(0.65, 0.80)
        precision = random.uniform(0.60, 0.75)
        recall = random.uniform(0.55, 0.70)
        f1 = random.uniform(0.57, 0.72)
    
    return {
        'model_type': model_type,
        'auc': auc,
        'precision': precision,
        'recall': recall,
        'f1_score': f1,
        'cv_auc_mean': auc - random.uniform(0.02, 0.05),
        'cv_auc_std': random.uniform(0.01, 0.03),
        'n_features': n_features,
        'cross_val_folds': 5
    }

# Simulate feature importance
def simulate_feature_importance(model_type, n_features):
    features = [
        'methylation_enhancer_CRXIP', 'methylation_SNP_upstream_PRK',
        'methylation_SNP_near_MYC', 'methylation_hypo_TMB_RNF',
        'methylation_hyper_TMB_PDL1', 'methylation_global_hypo_CpG',
        'fragmentomics_short_TFBS', 'fragmentomics_entropy',
        'fragmentomics_WPS_STAT1', 'fragmentomics_end_motif_CTCF',
        'variant_KRAS_mutation', 'variant_TMB', 'variant_DNA_burden',
        'expression_immune_activation', 'expression_t_cell_exhaustion'
    ]
    
    importance_scores = []
    for i, feature in enumerate(features[:n_features]):
        if model_type == 'XGBoost':
            importance = random.uniform(0.05, 0.15)
        elif model_type == 'RandomForest':
            importance = random.uniform(0.03, 0.12)
        else:
            importance = random.uniform(0.01, 0.08)
        
        importance_scores.append({
            'feature': feature,
            'importance_score': importance,
            'rank': i + 1,
            'model_type': model_type
        })
    
    return importance_scores

# Simulate predictions
def simulate_predictions(model_type, n_samples):
    predictions = []
    for i in range(n_samples):
        if model_type == 'XGBoost':
            prob_cancer = random.uniform(0.7, 0.95)
        elif model_type == 'RandomForest':
            prob_cancer = random.uniform(0.6, 0.90)
        else:
            prob_cancer = random.uniform(0.5, 0.85)
        
        predictions.append({
            'sample_id': f'SAMPLE_{i+1:03d}',
            'predicted_class': 'Cancer' if prob_cancer > 0.5 else 'Control',
            'probability_cancer': prob_cancer,
            'probability_control': 1 - prob_cancer,
            'model_type': model_type
        })
    
    return predictions

# Generate results
n_features = 15
n_samples = 100

# Performance metrics
performance = simulate_ml_training('${model_type}', n_features)
performance_df = pd.DataFrame([performance])
performance_df.to_csv('${model_type}_performance_metrics.tsv', sep='\t', index=False)

# Feature importance
feature_importance = simulate_feature_importance('${model_type}', n_features)
feature_importance_df = pd.DataFrame(feature_importance)
feature_importance_df.to_csv('${model_type}_feature_importance.tsv', sep='\t', index=False)

# Predictions
predictions = simulate_predictions('${model_type}', n_samples)
predictions_df = pd.DataFrame(predictions)
predictions_df.to_csv('${model_type}_predictions.tsv', sep='\t', index=False)

print(f'${model_type} model training completed')
print(f'AUC: {performance[\"auc\"]:.3f}')
print(f'Precision: {performance[\"precision\"]:.3f}')
print(f'Recall: {performance[\"recall\"]:.3f}')
print(f'F1-Score: {performance[\"f1_score\"]:.3f}')
"
    """
}

// Main biomarker discovery and ML training workflow
workflow {
    // Create test sample channel
    Channel.of(
        ["ML_TEST_001", "CRC", "Early", "Cancer"],
        ["ML_TEST_002", "CRC", "Early", "Control"],
        ["ML_TEST_003", "Lung", "Late", "Cancer"]
    )
    .set { ml_sample_channel }
    
    // Run biomarker discovery
    BIOMARKER_DISCOVERY(ml_sample_channel)
    
    // Collect biomarker discovery results
    BIOMARKER_DISCOVERY.out.discovery_log
        .collect()
        .set { discovery_logs }
    
    BIOMARKER_DISCOVERY.out.statistical_tests
        .collect()
        .set { statistical_tests }
    
    BIOMARKER_DISCOVERY.out.feature_importance
        .collect()
        .set { feature_importance }
    
    // Create ML model types
    Channel.of("XGBoost", "RandomForest", "LogisticRegression", "NaiveBayes")
        .set { model_types }
    
    // Run ML training for each model type
    ML_TRAINING(model_types)
    
    // Collect ML results
    ML_TRAINING.out.model_log
        .collect()
        .set { model_logs }
    
    ML_TRAINING.out.performance_metrics
        .collect()
        .set { performance_metrics }
    
    ML_TRAINING.out.feature_importance
        .collect()
        .set { ml_feature_importance }
    
    ML_TRAINING.out.predictions
        .collect()
        .set { predictions }
    
    // Print results
    discovery_logs
        .map { file -> 
            """
            Biomarker discovery log: ${file}
            """
        }
        .view()
    
    statistical_tests
        .map { file -> 
            """
            Statistical tests: ${file}
            """
        }
        .view()
    
    performance_metrics
        .map { file -> 
            """
            ML performance metrics: ${file}
            """
        }
        .view()
    
    predictions
        .map { file -> 
            """
            ML predictions: ${file}
            """
        }
        .view()
}

// Workflow completion
workflow.onComplete {
    log.info """
    ================================================
    CENTAUR Biomarker Discovery & ML Training completed successfully!
    ================================================
    """
}

// Workflow error handling
workflow.onError {
    log.error """
    ================================================
    CENTAUR Biomarker Discovery & ML Training failed!
    ================================================
    
    Error: ${workflow.errorMessage}
    """
    exit 1
}