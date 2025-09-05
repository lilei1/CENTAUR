#!/usr/bin/env nextflow

/*
 * Methylation Processing Workflow
 * Processes cfDNA EM-seq data for methylation and fragmentomics analysis
 */

nextflow.enable.dsl = 2

// Process definitions
process QUALITY_CONTROL_METHYLATION {
    tag "${sample_id}"
    label 'process_medium'
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(cfdna_file)
    path genome_fasta
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_qc_report.html"), emit: qc_report
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz"), emit: trimmed_reads
    
    script:
    """
    # Quality control with FastQC
    fastqc ${cfdna_file} -o . --noextract
    
    # Trim reads with Trim Galore for EM-seq protocol
    trim_galore \
        --paired \
        --quality 20 \
        --length 50 \
        --stringency 3 \
        --fastqc \
        --output_dir . \
        ${cfdna_file}
    
    # Generate MultiQC report
    multiqc . -n ${sample_id}_qc_report.html
    """
}

process ALIGN_METHYLATION {
    tag "${sample_id}"
    label 'process_high'
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(r1), path(r2)
    path genome_index
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_aligned.bam"), emit: aligned_bam
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_alignment_stats.txt"), emit: alignment_stats
    
    script:
    """
    # Align with Bismark for EM-seq protocol
    bismark \
        --genome ${genome_index} \
        --parallel 4 \
        --output_dir . \
        -1 ${r1} \
        -2 ${r2}
    
    # Rename output files
    mv *_bismark_bt2_pe.bam ${sample_id}_aligned.bam
    mv *_bismark_bt2_PE_report.txt ${sample_id}_alignment_stats.txt
    
    # Index BAM file
    samtools index ${sample_id}_aligned.bam
    """
}

process DEDUPLICATE_METHYLATION {
    tag "${sample_id}"
    label 'process_medium'
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(aligned_bam)
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_dedup.bam"), emit: dedup_bam
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_dedup_stats.txt"), emit: dedup_stats
    
    script:
    """
    # Deduplicate with Picard
    picard MarkDuplicates \
        INPUT=${aligned_bam} \
        OUTPUT=${sample_id}_dedup.bam \
        METRICS_FILE=${sample_id}_dedup_stats.txt \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true
    """
}

process CALL_METHYLATION {
    tag "${sample_id}"
    label 'process_medium'
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(dedup_bam)
    path genome_fasta
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_methylation.bedGraph"), emit: methylation_bedgraph
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_methylation_summary.txt"), emit: methylation_summary
    
    script:
    """
    # Call methylation with MethylDackel
    methylDackel extract \
        --CHG \
        --CHH \
        --CpG \
        --OT 0,0,0 \
        --OB 0,0,0 \
        --CTOT 0,0,0 \
        --CTOB 0,0,0 \
        ${genome_fasta} \
        ${dedup_bam} \
        --outputDir . \
        --outputPrefix ${sample_id}
    
    # Generate summary statistics
    methylDackel mbias \
        ${genome_fasta} \
        ${dedup_bam} \
        ${sample_id}_mbias.txt
    """
}

process EXTRACT_FRAGMENTOMICS {
    tag "${sample_id}"
    label 'process_medium'
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(dedup_bam)
    path genome_fasta
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_fragmentomics.tsv"), emit: fragmentomics_features
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_insert_size_histogram.png"), emit: insert_size_plot
    
    script:
    """
    # Extract fragmentomics features using custom Python script
    python3 -c "
    import pysam
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from collections import Counter
    
    # Read BAM file
    bamfile = pysam.AlignmentFile('${dedup_bam}', 'rb')
    
    # Extract fragment size distribution
    insert_sizes = []
    end_motifs = []
    wps_scores = []
    
    for read in bamfile:
        if read.is_proper_pair and read.is_read1:
            insert_size = read.template_length
            if insert_size > 0:
                insert_sizes.append(insert_size)
            
            # Extract end motifs (first 4 bases)
            if read.query_sequence:
                end_motif = read.query_sequence[:4]
                end_motifs.append(end_motif)
    
    bamfile.close()
    
    # Calculate fragmentomics features
    features = {
        'sample_id': '${sample_id}',
        'total_fragments': len(insert_sizes),
        'mean_insert_size': np.mean(insert_sizes) if insert_sizes else 0,
        'median_insert_size': np.median(insert_sizes) if insert_sizes else 0,
        'std_insert_size': np.std(insert_sizes) if insert_sizes else 0,
        'short_fragments_ratio': sum(1 for x in insert_sizes if x < 150) / len(insert_sizes) if insert_sizes else 0,
        'long_fragments_ratio': sum(1 for x in insert_sizes if x > 300) / len(insert_sizes) if insert_sizes else 0,
        'fragment_entropy': -sum((count/len(insert_sizes))**2 * np.log2((count/len(insert_sizes))**2) 
                                for count in Counter(insert_sizes).values()) if insert_sizes else 0
    }
    
    # Add end motif frequencies
    motif_counts = Counter(end_motifs)
    total_motifs = sum(motif_counts.values())
    for motif in ['CCCA', 'CCCG', 'CCCT', 'CCCC']:
        features[f'{motif}_frequency'] = motif_counts.get(motif, 0) / total_motifs if total_motifs > 0 else 0
    
    # Save features
    df = pd.DataFrame([features])
    df.to_csv('${sample_id}_fragmentomics.tsv', sep='\t', index=False)
    
    # Create insert size histogram
    plt.figure(figsize=(10, 6))
    plt.hist(insert_sizes, bins=50, alpha=0.7, edgecolor='black')
    plt.xlabel('Insert Size (bp)')
    plt.ylabel('Frequency')
    plt.title(f'Insert Size Distribution - ${sample_id}')
    plt.savefig('${sample_id}_insert_size_histogram.png', dpi=300, bbox_inches='tight')
    plt.close()
    "
    """
}

process CALL_DMR {
    tag "DMR_analysis"
    label 'process_high'
    
    input:
    path methylation_files
    path genome_fasta
    
    output:
    path("dmr_results.bed"), emit: dmr_results
    path("dmr_summary.txt"), emit: dmr_summary
    
    script:
    """
    # Call DMRs using methylKit or similar tool
    python3 -c "
    import pandas as pd
    import numpy as np
    from scipy import stats
    
    # This is a simplified DMR calling approach
    # In practice, you would use specialized tools like methylKit, DSS, or dmrseq
    
    print('Calling DMRs...')
    print('This is a placeholder for actual DMR calling implementation')
    print('In practice, use methylKit, DSS, or dmrseq for proper DMR analysis')
    
    # Create placeholder DMR results
    dmr_data = {
        'chr': ['chr1', 'chr2', 'chr3'],
        'start': [1000000, 2000000, 3000000],
        'end': [1001000, 2001000, 3001000],
        'pvalue': [0.001, 0.002, 0.003],
        'fdr': [0.01, 0.02, 0.03],
        'methylation_diff': [25.5, -30.2, 18.7],
        'region_type': ['enhancer', 'promoter', 'gene_body']
    }
    
    df = pd.DataFrame(dmr_data)
    df.to_csv('dmr_results.bed', sep='\t', index=False)
    
    # Summary
    with open('dmr_summary.txt', 'w') as f:
        f.write(f'Total DMRs: {len(df)}\\n')
        f.write(f'Hyper-methylated: {sum(df[\"methylation_diff\"] > 0)}\\n')
        f.write(f'Hypo-methylated: {sum(df[\"methylation_diff\"] < 0)}\\n')
    "
    """
}

// Workflow definition
workflow METHYLATION_PROCESSING {
    take:
    sample_channel    // channel: [sample_id, cancer_type, stage, sample_type, cfdna_file]
    genome_fasta      // path: reference genome
    genome_index      // path: genome index
    target_coverage   // val: target coverage
    
    main:
    // Quality control and trimming
    QUALITY_CONTROL_METHYLATION(sample_channel, genome_fasta)
    
    // Alignment
    ALIGN_METHYLATION(QUALITY_CONTROL_METHYLATION.out.trimmed_reads, genome_index)
    
    // Deduplication
    DEDUPLICATE_METHYLATION(ALIGN_METHYLATION.out.aligned_bam)
    
    // Methylation calling
    CALL_METHYLATION(DEDUPLICATE_METHYLATION.out.dedup_bam, genome_fasta)
    
    // Fragmentomics extraction
    EXTRACT_FRAGMENTOMICS(DEDUPLICATE_METHYLATION.out.dedup_bam, genome_fasta)
    
    // DMR calling (combines all samples)
    CALL_DMR(CALL_METHYLATION.out.methylation_bedgraph.collect(), genome_fasta)
    
    emit:
    methylation_features = CALL_METHYLATION.out.methylation_bedgraph
    fragmentomics_features = EXTRACT_FRAGMENTOMICS.out.fragmentomics_features
    dmr_results = CALL_DMR.out.dmr_results
    qc_reports = QUALITY_CONTROL_METHYLATION.out.qc_report
}
