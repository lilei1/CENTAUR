#!/usr/bin/env nextflow

/*
 * CENTAUR Multi-Omics Pipeline
 * A Scalable Platform for Early Cancer Biomarker Discovery
 * 
 * This pipeline processes cfDNA methylation, RNA-seq, and WES data
 * to identify biomarkers for early cancer detection
 */

nextflow.enable.dsl = 2

// Import sub-workflows
include { METHYLATION_PROCESSING } from '../processes/methylation.nf'
include { RNASEQ_PROCESSING } from '../processes/rnaseq.nf'
include { WES_PROCESSING } from '../processes/wes.nf'
include { MULTIOMICS_INTEGRATION } from '../processes/integration.nf'

// Pipeline parameters
params {
    // Input data
    input_dir = "$baseDir/data/raw"
    sample_sheet = "$baseDir/data/sample_sheet.csv"
    
    // Reference files
    genome_fasta = "$baseDir/references/hg38.fa"
    genome_index = "$baseDir/references/hg38"
    annotation_gtf = "$baseDir/references/hg38.gtf"
    
    // Output directory
    output_dir = "$baseDir/data/results"
    
    // Processing parameters
    methylation_coverage = 20
    rnaseq_reads = 50
    wes_coverage = 30
    
    // Multi-omics integration
    dmr_fdr_threshold = 0.05
    methylation_diff_threshold = 20
    tmb_threshold = 10
    
    // Resource allocation
    max_cpus = 16
    max_memory = '64.GB'
    max_time = '24.h'
    
    // Help information
    help = false
}

// Print help information
if (params.help) {
    log.info """
    CENTAUR Multi-Omics Pipeline
    
    Usage:
        nextflow run main.nf --input_dir <path> --sample_sheet <path>
    
    Parameters:
        --input_dir              Input directory containing raw sequencing data
        --sample_sheet           CSV file with sample metadata
        --genome_fasta           Reference genome FASTA file
        --genome_index           Reference genome index directory
        --annotation_gtf        Gene annotation GTF file
        --output_dir             Output directory for results
        --methylation_coverage   Target coverage for methylation analysis (default: 20)
        --rnaseq_reads           Target read count for RNA-seq (default: 50M)
        --wes_coverage           Target coverage for WES (default: 30)
        --dmr_fdr_threshold      FDR threshold for DMR calling (default: 0.05)
        --methylation_diff_threshold  Methylation difference threshold (default: 20%)
        --tmb_threshold          TMB threshold for filtering (default: 10)
        --max_cpus               Maximum number of CPUs to use (default: 16)
        --max_memory             Maximum memory to use (default: 64.GB)
        --max_time               Maximum time per process (default: 24.h)
        --help                   Show this help message
    
    Example:
        nextflow run main.nf --input_dir data/raw --sample_sheet data/sample_sheet.csv
    """
    exit 0
}

// Validate input parameters
if (!params.input_dir || !params.sample_sheet) {
    log.error "Missing required parameters: --input_dir and --sample_sheet"
    exit 1
}

// Create output directory
workflow {
    // Parse sample sheet
    Channel.fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> 
            [
                row.sample_id,
                row.cancer_type,
                row.stage,
                row.sample_type,
                row.cfDNA_file,
                row.rnaseq_file,
                row.wes_tumor_file,
                row.wes_normal_file
            ]
        }
        .set { sample_channel }
    
    // Process methylation data (cfDNA EM-seq)
    METHYLATION_PROCESSING (
        sample_channel,
        params.genome_fasta,
        params.genome_index,
        params.methylation_coverage
    )
    
    // Process RNA-seq data (PBMC)
    RNASEQ_PROCESSING (
        sample_channel,
        params.genome_fasta,
        params.genome_index,
        params.annotation_gtf,
        params.rnaseq_reads
    )
    
    // Process WES data (tumor + normal)
    WES_PROCESSING (
        sample_channel,
        params.genome_fasta,
        params.genome_index,
        params.wes_coverage
    )
    
    // Integrate multi-omics data and generate ML-ready features
    MULTIOMICS_INTEGRATION (
        METHYLATION_PROCESSING.out.methylation_features,
        METHYLATION_PROCESSING.out.fragmentomics_features,
        RNASEQ_PROCESSING.out.expression_features,
        RNASEQ_PROCESSING.out.immune_signatures,
        WES_PROCESSING.out.variant_features,
        WES_PROCESSING.out.tmb_features,
        sample_channel,
        params.dmr_fdr_threshold,
        params.methylation_diff_threshold,
        params.tmb_threshold
    )
    
    // Generate final reports
    MULTIOMICS_INTEGRATION.out.final_features
        .map { features -> 
            """
            Sample: ${features.sample_id}
            Cancer Type: ${features.cancer_type}
            Stage: ${features.stage}
            Methylation Features: ${features.methylation_count}
            Fragmentomics Features: ${features.fragmentomics_count}
            Expression Features: ${features.expression_count}
            Variant Features: ${features.variant_count}
            Total Features: ${features.total_features}
            """
        }
        .collect()
        .set { final_report }
}

// Workflow completion
workflow.onComplete {
    log.info """
    ================================================
    CENTAUR Multi-Omics Pipeline completed successfully!
    ================================================
    
    Results are available in: ${params.output_dir}
    
    Pipeline Summary:
    - Methylation analysis: ${METHYLATION_PROCESSING.out.methylation_features.count()} samples
    - RNA-seq analysis: ${RNASEQ_PROCESSING.out.expression_features.count()} samples  
    - WES analysis: ${WES_PROCESSING.out.variant_features.count()} samples
    - Multi-omics integration: ${MULTIOMICS_INTEGRATION.out.final_features.count()} samples
    
    Next steps:
    1. Review QC reports in ${params.output_dir}/qc/
    2. Examine biomarker features in ${params.output_dir}/features/
    3. Run downstream ML analysis using the generated feature matrices
    """
}

// Workflow error handling
workflow.onError {
    log.error """
    ================================================
    CENTAUR Multi-Omics Pipeline failed!
    ================================================
    
    Error: ${workflow.errorMessage}
    
    Please check the logs and try again.
    """
    exit 1
}
