#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Multi-omics integration test process
process MULTIOMICS_INTEGRATION_TEST {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type)
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_integration_test.txt"), emit: test_output
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_harmonized_features.tsv"), emit: harmonized_features
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_biomarker_features.tsv"), emit: biomarker_features
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_ml_features.tsv"), emit: ml_features
    
    script:
    """
    echo "Processing multi-omics integration for ${sample_id}" > ${sample_id}_integration_test.txt
    echo "Sample type: ${sample_type}" >> ${sample_id}_integration_test.txt
    echo "Cancer type: ${cancer_type}" >> ${sample_id}_integration_test.txt
    echo "Stage: ${stage}" >> ${sample_id}_integration_test.txt
    echo "Processing completed at: \$(date)" >> ${sample_id}_integration_test.txt
    
    # Create Python script for multi-omics integration
    cat > multiomics_integration.py << 'EOF'
import pandas as pd
import numpy as np
import random

print(f'Processing ${sample_id} - ${sample_type}')

# Simulate features from each omics type based on previous test results
# Methylation features (from EM-seq test)
methylation_features = {
    'methylation_global': random.uniform(0.3, 0.7),
    'methylation_immune_enhancers': random.uniform(0.2, 0.8),
    'methylation_cancer_drivers': random.uniform(0.1, 0.9),
    'methylation_super_enhancers': random.uniform(0.2, 0.7)
}

# Fragmentomics features (from methylation test)
fragmentomics_features = {
    'fragmentomics_short_ratio': random.uniform(0.1, 0.3),
    'fragmentomics_entropy': random.uniform(2.5, 3.5),
    'fragmentomics_wps_score': random.uniform(0.1, 0.9),
    'fragmentomics_end_motif_bias': random.uniform(0.05, 0.15)
}

# Expression features (from RNA-seq test)
if '${sample_type}' == 'Cancer':
    expression_features = {
        'expression_ifn_gamma': random.uniform(1000, 2000),  # High in cancer
        'expression_t_cell_exhaustion': random.uniform(700, 1200),
        'expression_cytotoxic_t_cells': random.uniform(500, 800),
        'expression_b_cells': random.uniform(300, 600),
        'expression_nk_cells': random.uniform(500, 700)
    }
else:
    expression_features = {
        'expression_ifn_gamma': random.uniform(600, 1000),  # Lower in control
        'expression_t_cell_exhaustion': random.uniform(700, 1000),
        'expression_cytotoxic_t_cells': random.uniform(700, 900),
        'expression_b_cells': random.uniform(600, 800),
        'expression_nk_cells': random.uniform(600, 800)
    }

# Variant features (from WES test)
if '${sample_type}' == 'Cancer':
    variant_features = {
        'variant_tmb': random.uniform(0.5, 1.5),  # Higher TMB in cancer
        'variant_total_mutations': random.randint(15, 35),
        'variant_snv_count': random.randint(10, 25),
        'variant_indel_count': random.randint(5, 15),
        'variant_cnv_gains': random.randint(1, 3),
        'variant_cnv_losses': random.randint(1, 3),
        'variant_key_mutations': random.randint(1, 4)
    }
else:
    variant_features = {
        'variant_tmb': random.uniform(0.1, 0.5),  # Lower TMB in control
        'variant_total_mutations': random.randint(2, 10),
        'variant_snv_count': random.randint(2, 8),
        'variant_indel_count': random.randint(0, 5),
        'variant_cnv_gains': random.randint(0, 2),
        'variant_cnv_losses': random.randint(0, 2),
        'variant_key_mutations': random.randint(0, 2)
    }

# Add key mutation features
key_mutations = ['TP53', 'KRAS', 'PIK3CA', 'APC', 'EGFR', 'BRAF', 'MYC', 'TERT']
for gene in key_mutations:
    variant_features[f'variant_{gene}_mutated'] = random.random() < 0.3

# Combine all features
all_features = {
    'sample_id': '${sample_id}',
    'cancer_type': '${cancer_type}',
    'stage': '${stage}',
    'sample_type': '${sample_type}',
    **methylation_features,
    **fragmentomics_features,
    **expression_features,
    **variant_features
}

# Create harmonized features DataFrame
harmonized_df = pd.DataFrame([all_features])
harmonized_df.to_csv('${sample_id}_harmonized_features.tsv', sep='\t', index=False)

# Extract biomarker features based on experimental design
biomarker_features = {
    'sample_id': '${sample_id}',
    'cancer_type': '${cancer_type}',
    'stage': '${stage}',
    'sample_type': '${sample_type}',
    
    # Methylation biomarkers (6 features)
    'methylation_enhancer_CRXIP': methylation_features['methylation_immune_enhancers'],
    'methylation_SNP_upstream_PRK': methylation_features['methylation_cancer_drivers'],
    'methylation_SNP_near_MYC': methylation_features['methylation_cancer_drivers'],
    'methylation_hypo_TMB_RNF': methylation_features['methylation_global'],
    'methylation_hyper_TMB_PDL1': methylation_features['methylation_immune_enhancers'],
    'methylation_global_hypo_CpG': methylation_features['methylation_global'],
    
    # Fragmentomics biomarkers (4 features)
    'fragmentomics_short_TFBS': fragmentomics_features['fragmentomics_short_ratio'],
    'fragmentomics_entropy': fragmentomics_features['fragmentomics_entropy'],
    'fragmentomics_WPS_STAT1': fragmentomics_features['fragmentomics_wps_score'],
    'fragmentomics_end_motif_CTCF': fragmentomics_features['fragmentomics_end_motif_bias'],
    
    # Somatic mutation biomarkers (3 features)
    'variant_KRAS_mutation': variant_features['variant_KRAS_mutated'],
    'variant_TMB': variant_features['variant_tmb'],
    'variant_DNA_burden': variant_features['variant_total_mutations'],
    
    # RNA expression biomarkers (2 features)
    'expression_immune_activation': expression_features['expression_ifn_gamma'],
    'expression_t_cell_exhaustion': expression_features['expression_t_cell_exhaustion']
}

# Create biomarker features DataFrame
biomarker_df = pd.DataFrame([biomarker_features])
biomarker_df.to_csv('${sample_id}_biomarker_features.tsv', sep='\t', index=False)

# Create ML-ready features (15 key biomarkers)
ml_features = {
    'sample_id': '${sample_id}',
    'cancer_type': '${cancer_type}',
    'stage': '${stage}',
    'sample_type': '${sample_type}',
    
    # 15 key biomarkers from experimental design
    'methylation_enhancer_CRXIP': biomarker_features['methylation_enhancer_CRXIP'],
    'methylation_SNP_upstream_PRK': biomarker_features['methylation_SNP_upstream_PRK'],
    'methylation_SNP_near_MYC': biomarker_features['methylation_SNP_near_MYC'],
    'methylation_hypo_TMB_RNF': biomarker_features['methylation_hypo_TMB_RNF'],
    'methylation_hyper_TMB_PDL1': biomarker_features['methylation_hyper_TMB_PDL1'],
    'methylation_global_hypo_CpG': biomarker_features['methylation_global_hypo_CpG'],
    
    'fragmentomics_short_TFBS': biomarker_features['fragmentomics_short_TFBS'],
    'fragmentomics_entropy': biomarker_features['fragmentomics_entropy'],
    'fragmentomics_WPS_STAT1': biomarker_features['fragmentomics_WPS_STAT1'],
    'fragmentomics_end_motif_CTCF': biomarker_features['fragmentomics_end_motif_CTCF'],
    
    'variant_KRAS_mutation': biomarker_features['variant_KRAS_mutation'],
    'variant_TMB': biomarker_features['variant_TMB'],
    'variant_DNA_burden': biomarker_features['variant_DNA_burden'],
    
    'expression_immune_activation': biomarker_features['expression_immune_activation'],
    'expression_t_cell_exhaustion': biomarker_features['expression_t_cell_exhaustion']
}

# Create ML features DataFrame
ml_df = pd.DataFrame([ml_features])
ml_df.to_csv('${sample_id}_ml_features.tsv', sep='\t', index=False)

print(f'Generated multi-omics integration for ${sample_id}')
print(f'Total harmonized features: {len(all_features)}')
print(f'Biomarker features: {len(biomarker_features)}')
print(f'ML-ready features: {len(ml_features)}')
EOF

    # Run Python script
    python3 multiomics_integration.py
    """
}

// Test workflow
workflow {
    // Create integrated sample channel
    Channel.of(
        ["INTEGRATION_TEST_001", "CRC", "Early", "Cancer"],
        ["INTEGRATION_TEST_002", "CRC", "Early", "Control"],
        ["INTEGRATION_TEST_003", "Lung", "Late", "Cancer"]
    )
    .set { integration_sample_channel }
    
    // Run multi-omics integration test
    MULTIOMICS_INTEGRATION_TEST(integration_sample_channel)
    
    // Collect results
    MULTIOMICS_INTEGRATION_TEST.out.test_output
        .collect()
        .set { test_results }
    
    MULTIOMICS_INTEGRATION_TEST.out.harmonized_features
        .collect()
        .set { harmonized_results }
    
    MULTIOMICS_INTEGRATION_TEST.out.biomarker_features
        .collect()
        .set { biomarker_results }
    
    MULTIOMICS_INTEGRATION_TEST.out.ml_features
        .collect()
        .set { ml_results }
    
    // Print results
    test_results
        .map { file -> 
            """
            Integration test file: ${file}
            """
        }
        .view()
    
    harmonized_results
        .map { file -> 
            """
            Harmonized features file: ${file}
            """
        }
        .view()
    
    biomarker_results
        .map { file -> 
            """
            Biomarker features file: ${file}
            """
        }
        .view()
    
    ml_results
        .map { file -> 
            """
            ML-ready features file: ${file}
            """
        }
        .view()
}

// Workflow completion
workflow.onComplete {
    log.info """
    ================================================
    CENTAUR Multi-Omics Integration Test completed successfully!
    ================================================
    """
}

// Workflow error handling
workflow.onError {
    log.error """
    ================================================
    CENTAUR Multi-Omics Integration Test failed!
    ================================================
    
    Error: ${workflow.errorMessage}
    """
    exit 1
}