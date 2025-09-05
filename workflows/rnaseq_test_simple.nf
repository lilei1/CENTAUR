#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Simple RNA-seq test process
process RNASEQ_TEST {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(rnaseq_file)
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_rnaseq_test.txt"), emit: test_output
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_gene_counts.txt"), emit: gene_counts
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_immune_signatures.tsv"), emit: immune_signatures
    
    script:
    """
    echo "Processing RNA-seq data for ${sample_id}" > ${sample_id}_rnaseq_test.txt
    echo "Sample type: ${sample_type}" >> ${sample_id}_rnaseq_test.txt
    echo "Cancer type: ${cancer_type}" >> ${sample_id}_rnaseq_test.txt
    echo "Stage: ${stage}" >> ${sample_id}_rnaseq_test.txt
    echo "Input file: ${rnaseq_file}" >> ${sample_id}_rnaseq_test.txt
    echo "Processing completed at: \$(date)" >> ${sample_id}_rnaseq_test.txt
    
    # Simulate gene expression analysis
    python3 -c "
    import pandas as pd
    import numpy as np
    import random
    
    # Define mock genes
    genes = [
        'STAT1', 'IRF1', 'CXCL10', 'CXCL9', 'IDO1',  # IFN-gamma response
        'PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT',  # T-cell exhaustion
        'GZMA', 'GZMB', 'PRF1', 'GNLY', 'NKG7',       # Cytotoxic T-cells
        'CD19', 'CD20', 'CD22', 'CD79A', 'CD79B',      # B-cells
        'KIR2DL1', 'KIR2DL2', 'KIR2DL3', 'KIR2DL4',   # NK-cells
        'GAPDH', 'ACTB', 'TUBB', 'RPL13A', 'RPS18'    # Housekeeping genes
    ]
    
    # Simulate gene expression counts based on sample type
    gene_counts = {}
    for gene in genes:
        if '${sample_type}' == 'Cancer':
            # Cancer samples: higher immune gene expression, lower housekeeping
            if gene in ['STAT1', 'IRF1', 'CXCL10', 'CXCL9', 'IDO1', 'PDCD1', 'CTLA4']:
                count = random.randint(500, 2000)  # High immune gene expression
            elif gene in ['GAPDH', 'ACTB', 'TUBB', 'RPL13A', 'RPS18']:
                count = random.randint(100, 500)  # Lower housekeeping
            else:
                count = random.randint(200, 800)  # Medium expression
        else:
            # Control samples: balanced expression
            count = random.randint(300, 1000)
        
        gene_counts[gene] = count
    
    # Create gene counts DataFrame
    counts_df = pd.DataFrame(list(gene_counts.items()), columns=['Gene', 'Count'])
    counts_df.to_csv('${sample_id}_gene_counts.txt', sep='\t', index=False)
    
    # Calculate immune signatures
    immune_signatures = {
        'sample_id': '${sample_id}',
        'IFN_gamma_response': np.mean([gene_counts[g] for g in ['STAT1', 'IRF1', 'CXCL10', 'CXCL9', 'IDO1']]),
        'T_cell_exhaustion': np.mean([gene_counts[g] for g in ['PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT']]),
        'Cytotoxic_T_cells': np.mean([gene_counts[g] for g in ['GZMA', 'GZMB', 'PRF1', 'GNLY', 'NKG7']]),
        'B_cells': np.mean([gene_counts[g] for g in ['CD19', 'CD20', 'CD22', 'CD79A', 'CD79B']]),
        'NK_cells': np.mean([gene_counts[g] for g in ['KIR2DL1', 'KIR2DL2', 'KIR2DL3', 'KIR2DL4']]),
        'Housekeeping': np.mean([gene_counts[g] for g in ['GAPDH', 'ACTB', 'TUBB', 'RPL13A', 'RPS18']])
    }
    
    # Create immune signatures DataFrame
    immune_df = pd.DataFrame([immune_signatures])
    immune_df.to_csv('${sample_id}_immune_signatures.tsv', sep='\t', index=False)
    
    print(f'Generated gene expression data for ${sample_id}')
    print(f'Total genes: {len(genes)}')
    print(f'Immune signatures calculated')
    "
    """
}

// Test workflow
workflow {
    // Create test sample channel
    Channel.of(
        ["RNA_TEST_001", "CRC", "Early", "Cancer", file("${baseDir}/data/test_data/rnaseq/RNA_TEST_001/RNA_TEST_001_R1.fastq.gz")],
        ["RNA_TEST_002", "CRC", "Early", "Control", file("${baseDir}/data/test_data/rnaseq/RNA_TEST_002/RNA_TEST_002_R1.fastq.gz")],
        ["RNA_TEST_003", "Lung", "Late", "Cancer", file("${baseDir}/data/test_data/rnaseq/RNA_TEST_003/RNA_TEST_003_R1.fastq.gz")]
    )
    .set { sample_channel }
    
    // Run RNA-seq test
    RNASEQ_TEST(sample_channel)
    
    // Collect results
    RNASEQ_TEST.out.test_output
        .collect()
        .set { test_results }
    
    RNASEQ_TEST.out.gene_counts
        .collect()
        .set { gene_counts_results }
    
    RNASEQ_TEST.out.immune_signatures
        .collect()
        .set { immune_signatures_results }
    
    // Print results
    test_results
        .map { file -> 
            """
            Test file: ${file}
            """
        }
        .view()
    
    gene_counts_results
        .map { file -> 
            """
            Gene counts file: ${file}
            """
        }
        .view()
    
    immune_signatures_results
        .map { file -> 
            """
            Immune signatures file: ${file}
            """
        }
        .view()
}

// Workflow completion
workflow.onComplete {
    log.info """
    ================================================
    CENTAUR RNA-seq Test completed successfully!
    ================================================
    """
}

// Workflow error handling
workflow.onError {
    log.error """
    ================================================
    CENTAUR RNA-seq Test failed!
    ================================================
    
    Error: ${workflow.errorMessage}
    """
    exit 1
}
