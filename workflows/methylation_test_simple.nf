#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Simple methylation test process
process METHYLATION_TEST {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(cfdna_file)
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_methylation_test.txt"), emit: test_output
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_fragmentomics_test.tsv"), emit: fragmentomics
    
    script:
    """
    echo "Processing methylation data for ${sample_id}" > ${sample_id}_methylation_test.txt
    echo "Sample type: ${sample_type}" >> ${sample_id}_methylation_test.txt
    echo "Cancer type: ${cancer_type}" >> ${sample_id}_methylation_test.txt
    echo "Stage: ${stage}" >> ${sample_id}_methylation_test.txt
    echo "Input file: ${cfdna_file}" >> ${sample_id}_methylation_test.txt
    echo "Processing completed at: \$(date)" >> ${sample_id}_methylation_test.txt
    
    # Simulate fragmentomics analysis
    python3 -c "
    import pandas as pd
    import numpy as np
    import random
    
    # Simulate fragmentomics features
    features = {
        'sample_id': '${sample_id}',
        'total_fragments': random.randint(800, 1200),
        'mean_insert_size': random.uniform(160, 180),
        'median_insert_size': random.uniform(165, 175),
        'std_insert_size': random.uniform(20, 30),
        'short_fragments_ratio': random.uniform(0.1, 0.3),
        'long_fragments_ratio': random.uniform(0.05, 0.15),
        'fragment_entropy': random.uniform(2.5, 3.5),
        'CCCA_frequency': random.uniform(0.05, 0.15),
        'CCCG_frequency': random.uniform(0.03, 0.12),
        'CCCT_frequency': random.uniform(0.02, 0.10),
        'CCCC_frequency': random.uniform(0.01, 0.08)
    }
    
    df = pd.DataFrame([features])
    df.to_csv('${sample_id}_fragmentomics_test.tsv', sep='\t', index=False)
    print(f'Generated fragmentomics features for ${sample_id}')
    "
    """
}

// Test workflow
workflow {
    // Create test sample channel
    Channel.of(
        ["TEST_001", "CRC", "Early", "Cancer", file("${baseDir}/data/test_data/methylation/TEST_001/TEST_001_R1.fastq.gz")],
        ["TEST_002", "CRC", "Early", "Control", file("${baseDir}/data/test_data/methylation/TEST_002/TEST_002_R1.fastq.gz")],
        ["TEST_003", "Lung", "Late", "Cancer", file("${baseDir}/data/test_data/methylation/TEST_003/TEST_003_R1.fastq.gz")]
    )
    .set { sample_channel }
    
    // Run methylation test
    METHYLATION_TEST(sample_channel)
    
    // Collect results
    METHYLATION_TEST.out.test_output
        .collect()
        .set { test_results }
    
    METHYLATION_TEST.out.fragmentomics
        .collect()
        .set { fragmentomics_results }
    
    // Print results
    test_results
        .map { file -> 
            """
            Test file: ${file}
            Content:
            ${file.text}
            """
        }
        .view()
    
    fragmentomics_results
        .map { file -> 
            """
            Fragmentomics file: ${file}
            """
        }
        .view()
}

// Workflow completion
workflow.onComplete {
    log.info """
    ================================================
    CENTAUR Methylation Test completed successfully!
    ================================================
    """
}

// Workflow error handling
workflow.onError {
    log.error """
    ================================================
    CENTAUR Methylation Test failed!
    ================================================
    
    Error: ${workflow.errorMessage}
    """
    exit 1
}
