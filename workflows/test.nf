#!/usr/bin/env nextflow

/*
 * CENTAUR Multi-Omics Pipeline - Test Workflow
 * Simplified test workflow for validation
 */

nextflow.enable.dsl = 2

// Test parameters
params {
    test_mode = true
    test_samples = 3
    test_output_dir = "$baseDir/data/test_results"
}

// Simple test process
process TEST_PROCESS {
    tag "${sample_id}"
    
    input:
    val sample_id
    
    output:
    path("${sample_id}_test_output.txt"), emit: test_output
    
    script:
    """
    echo "Testing sample: ${sample_id}" > ${sample_id}_test_output.txt
    echo "Test completed successfully" >> ${sample_id}_test_output.txt
    echo "Timestamp: $(date)" >> ${sample_id}_test_output.txt
    """
}

// Test workflow
workflow {
    // Create test sample channel
    Channel.of("SAMPLE_001", "SAMPLE_002", "SAMPLE_003")
        .set { test_samples }
    
    // Run test process
    TEST_PROCESS(test_samples)
    
    // Collect results
    TEST_PROCESS.out.test_output
        .collect()
        .set { test_results }
    
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
}

// Workflow completion
workflow.onComplete {
    log.info """
    ================================================
    CENTAUR Test Workflow completed successfully!
    ================================================
    
    Test results are available in: ${params.test_output_dir}
    """
}

// Workflow error handling
workflow.onError {
    log.error """
    ================================================
    CENTAUR Test Workflow failed!
    ================================================
    
    Error: ${workflow.errorMessage}
    """
    exit 1
}
