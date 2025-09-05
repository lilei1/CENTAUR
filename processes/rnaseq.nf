#!/usr/bin/env nextflow

/*
 * RNA-seq Processing Workflow
 * Processes bulk RNA-seq data from PBMC samples for expression analysis and immune profiling
 */

nextflow.enable.dsl = 2

// Process definitions
process QUALITY_CONTROL_RNA {
    tag "${sample_id}"
    label 'process_medium'
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(rnaseq_file)
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_qc_report.html"), emit: qc_report
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz"), emit: trimmed_reads
    
    script:
    """
    # Quality control with FastQC
    fastqc ${rnaseq_file} -o . --noextract
    
    # Trim reads with Trim Galore
    trim_galore \
        --paired \
        --quality 20 \
        --length 50 \
        --stringency 3 \
        --fastqc \
        --output_dir . \
        ${rnaseq_file}
    
    # Generate MultiQC report
    multiqc . -n ${sample_id}_qc_report.html
    """
}

process ALIGN_RNA {
    tag "${sample_id}"
    label 'process_high'
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(r1), path(r2)
    path genome_index
    path annotation_gtf
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_aligned.bam"), emit: aligned_bam
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_alignment_stats.txt"), emit: alignment_stats
    
    script:
    """
    # Align with STAR
    STAR \
        --genomeDir ${genome_index} \
        --readFilesIn ${r1} ${r2} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${sample_id}_ \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --sjdbGTFfile ${annotation_gtf} \
        --runThreadN 8 \
        --outSAMattributes Standard
    
    # Rename output files
    mv ${sample_id}_Aligned.sortedByCoord.out.bam ${sample_id}_aligned.bam
    mv ${sample_id}_Log.final.out ${sample_id}_alignment_stats.txt
    
    # Index BAM file
    samtools index ${sample_id}_aligned.bam
    """
}

process QUANTIFY_GENES {
    tag "${sample_id}"
    label 'process_medium'
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(aligned_bam)
    path annotation_gtf
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_gene_counts.txt"), emit: gene_counts
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_quantification_stats.txt"), emit: quant_stats
    
    script:
    """
    # Quantify gene expression with featureCounts
    featureCounts \
        -a ${annotation_gtf} \
        -o ${sample_id}_gene_counts.txt \
        -T 4 \
        -p \
        -B \
        -C \
        ${aligned_bam}
    
    # Generate quantification statistics
    echo "Gene quantification completed for ${sample_id}" > ${sample_id}_quantification_stats.txt
    echo "Total genes quantified: $(tail -n +3 ${sample_id}_gene_counts.txt | wc -l)" >> ${sample_id}_quantification_stats.txt
    """
}

process NORMALIZE_EXPRESSION {
    tag "normalization"
    label 'process_medium'
    
    input:
    path gene_counts_files
    path sample_metadata
    
    output:
    path("normalized_counts.tsv"), emit: normalized_counts
    path("normalization_report.txt"), emit: normalization_report
    
    script:
    """
    # Normalize expression data with DESeq2
    python3 -c "
    import pandas as pd
    import numpy as np
    from scipy import stats
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # Read gene count files
    count_files = '${gene_counts_files}'.split()
    sample_ids = []
    count_data = {}
    
    for file in count_files:
        sample_id = file.split('/')[-1].replace('_gene_counts.txt', '')
        sample_ids.append(sample_id)
        
        # Read counts (skip first 2 lines which are comments)
        df = pd.read_csv(file, sep='\t', skiprows=1)
        count_data[sample_id] = df.iloc[:, -1]  # Last column contains counts
    
    # Create count matrix
    gene_ids = df['Geneid'].values
    count_matrix = pd.DataFrame(count_data, index=gene_ids)
    
    # Basic normalization (CPM)
    cpm_matrix = count_matrix.div(count_matrix.sum(axis=0), axis=1) * 1e6
    
    # Log2 transformation (add pseudocount to avoid log(0))
    log2_cpm_matrix = np.log2(cpm_matrix + 1)
    
    # Save normalized counts
    log2_cpm_matrix.to_csv('normalized_counts.tsv', sep='\t')
    
    # Generate normalization report
    with open('normalization_report.txt', 'w') as f:
        f.write(f'Total samples: {len(sample_ids)}\\n')
        f.write(f'Total genes: {len(gene_ids)}\\n')
        f.write(f'Mean counts per sample: {count_matrix.sum(axis=0).mean():.0f}\\n')
        f.write(f'Median counts per sample: {count_matrix.sum(axis=0).median():.0f}\\n')
    
    # Create QC plots
    plt.figure(figsize=(12, 8))
    
    # Library size distribution
    plt.subplot(2, 2, 1)
    plt.hist(count_matrix.sum(axis=0), bins=20, alpha=0.7)
    plt.xlabel('Total Counts')
    plt.ylabel('Number of Samples')
    plt.title('Library Size Distribution')
    
    # Gene expression distribution
    plt.subplot(2, 2, 2)
    plt.hist(log2_cpm_matrix.values.flatten(), bins=50, alpha=0.7)
    plt.xlabel('Log2 CPM')
    plt.ylabel('Frequency')
    plt.title('Gene Expression Distribution')
    
    # Sample correlation heatmap
    plt.subplot(2, 2, 3)
    correlation_matrix = log2_cpm_matrix.corr()
    sns.heatmap(correlation_matrix, cmap='coolwarm', center=0)
    plt.title('Sample Correlation Matrix')
    
    # PCA plot
    plt.subplot(2, 2, 4)
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(log2_cpm_matrix.T)
    plt.scatter(pca_result[:, 0], pca_result[:, 1], alpha=0.7)
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
    plt.title('PCA Plot')
    
    plt.tight_layout()
    plt.savefig('normalization_qc_plots.png', dpi=300, bbox_inches='tight')
    plt.close()
    "
    """
}

process DIFFERENTIAL_EXPRESSION {
    tag "DE_analysis"
    label 'process_medium'
    
    input:
    path normalized_counts
    path sample_metadata
    
    output:
    path("differential_expression_results.tsv"), emit: de_results
    path("de_summary.txt"), emit: de_summary
    
    script:
    """
    # Perform differential expression analysis
    python3 -c "
    import pandas as pd
    import numpy as np
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    
    # Read normalized counts
    expr_data = pd.read_csv('${normalized_counts}', sep='\t', index_col=0)
    
    # Read sample metadata (simplified for this example)
    # In practice, you would read the actual sample sheet
    sample_ids = expr_data.columns.tolist()
    
    # Create mock metadata for demonstration
    # In practice, this would come from your sample sheet
    metadata = {
        'sample_id': sample_ids,
        'cancer_type': ['CRC'] * (len(sample_ids) // 3) + ['Lung'] * (len(sample_ids) // 3) + ['Breast'] * (len(sample_ids) - 2 * (len(sample_ids) // 3)),
        'stage': ['Early'] * (len(sample_ids) // 2) + ['Late'] * (len(sample_ids) - len(sample_ids) // 2),
        'sample_type': ['Cancer'] * (len(sample_ids) // 2) + ['Control'] * (len(sample_ids) - len(sample_ids) // 2)
    }
    
    metadata_df = pd.DataFrame(metadata)
    
    # Perform differential expression analysis
    de_results = []
    
    for gene in expr_data.index:
        gene_expr = expr_data.loc[gene]
        
        # Compare cancer vs control
        cancer_expr = gene_expr[metadata_df['sample_type'] == 'Cancer']
        control_expr = gene_expr[metadata_df['sample_type'] == 'Control']
        
        if len(cancer_expr) > 1 and len(control_expr) > 1:
            # t-test
            t_stat, p_value = stats.ttest_ind(cancer_expr, control_expr)
            
            # Calculate fold change
            mean_cancer = cancer_expr.mean()
            mean_control = control_expr.mean()
            fold_change = mean_cancer / mean_control if mean_control > 0 else np.inf
            
            de_results.append({
                'gene_id': gene,
                'mean_cancer': mean_cancer,
                'mean_control': mean_control,
                'fold_change': fold_change,
                'log2_fold_change': np.log2(fold_change) if fold_change != np.inf else np.inf,
                'p_value': p_value,
                't_statistic': t_stat
            })
    
    # Convert to DataFrame
    de_df = pd.DataFrame(de_results)
    
    # Multiple testing correction
    if len(de_df) > 0:
        de_df['fdr'] = multipletests(de_df['p_value'], method='fdr')[1]
    else:
        de_df['fdr'] = []
    
    # Save results
    de_df.to_csv('differential_expression_results.tsv', sep='\t', index=False)
    
    # Generate summary
    with open('de_summary.txt', 'w') as f:
        f.write(f'Total genes tested: {len(de_df)}\\n')
        f.write(f'Significant genes (FDR < 0.05): {sum(de_df[\"fdr\"] < 0.05) if len(de_df) > 0 else 0}\\n')
        f.write(f'Upregulated genes: {sum(de_df[\"log2_fold_change\"] > 1) if len(de_df) > 0 else 0}\\n')
        f.write(f'Downregulated genes: {sum(de_df[\"log2_fold_change\"] < -1) if len(de_df) > 0 else 0}\\n')
    "
    """
}

process CALCULATE_IMMUNE_SIGNATURES {
    tag "immune_signatures"
    label 'process_medium'
    
    input:
    path normalized_counts
    path sample_metadata
    
    output:
    path("immune_signatures.tsv"), emit: immune_signatures
    path("immune_signature_summary.txt"), emit: immune_summary
    
    script:
    """
    # Calculate immune signatures using gene expression data
    python3 -c "
    import pandas as pd
    import numpy as np
    
    # Read normalized counts
    expr_data = pd.read_csv('${normalized_counts}', sep='\t', index_col=0)
    
    # Define immune signature gene sets
    immune_signatures = {
        'IFN_gamma_response': [
            'STAT1', 'IRF1', 'CXCL10', 'CXCL9', 'IDO1', 'GBP1', 'GBP2', 'GBP3',
            'GBP4', 'GBP5', 'GBP6', 'GBP7', 'IFI44', 'IFI44L', 'IFI6', 'IFI27',
            'IFI35', 'IFI16', 'IFIT1', 'IFIT2', 'IFIT3', 'IFIT5', 'ISG15', 'ISG20'
        ],
        'T_cell_exhaustion': [
            'PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT', 'CD244', 'CD160',
            'BTLA', 'CD200', 'CD200R1', 'CD276', 'CD274', 'PDCD1LG2'
        ],
        'Cytotoxic_T_cells': [
            'GZMA', 'GZMB', 'GZMH', 'GZMK', 'GZMM', 'PRF1', 'GNLY', 'NKG7',
            'KLRD1', 'KLRF1', 'KLRG1', 'KLRK1', 'KLRB1'
        ],
        'B_cells': [
            'CD19', 'CD20', 'CD22', 'CD79A', 'CD79B', 'MS4A1', 'PAX5', 'IGHM',
            'IGHA1', 'IGHA2', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4'
        ],
        'NK_cells': [
            'KIR2DL1', 'KIR2DL2', 'KIR2DL3', 'KIR2DL4', 'KIR2DL5A', 'KIR2DL5B',
            'KIR3DL1', 'KIR3DL2', 'KIR3DL3', 'KIR3DX1', 'KIR3DX1', 'KIR3DX1'
        ]
    }
    
    # Calculate signature scores for each sample
    signature_scores = {}
    
    for sample in expr_data.columns:
        sample_scores = {}
        
        for signature_name, genes in immune_signatures.items():
            # Find genes that are present in the expression data
            available_genes = [g for g in genes if g in expr_data.index]
            
            if available_genes:
                # Calculate mean expression of signature genes
                signature_expr = expr_data.loc[available_genes, sample]
                sample_scores[signature_name] = signature_expr.mean()
            else:
                sample_scores[signature_name] = 0
        
        signature_scores[sample] = sample_scores
    
    # Convert to DataFrame
    signature_df = pd.DataFrame(signature_scores).T
    signature_df.index.name = 'sample_id'
    
    # Save results
    signature_df.to_csv('immune_signatures.tsv', sep='\t')
    
    # Generate summary
    with open('immune_signature_summary.txt', 'w') as f:
        f.write(f'Total samples: {len(signature_df)}\\n')
        f.write(f'Total signatures: {len(signature_df.columns)}\\n')
        f.write('\\nSignature statistics:\\n')
        for sig in signature_df.columns:
            f.write(f'{sig}: mean={signature_df[sig].mean():.2f}, std={signature_df[sig].std():.2f}\\n')
    "
    """
}

// Workflow definition
workflow RNASEQ_PROCESSING {
    take:
    sample_channel    // channel: [sample_id, cancer_type, stage, sample_type, rnaseq_file]
    genome_fasta      // path: reference genome
    genome_index      // path: genome index
    annotation_gtf    // path: gene annotation
    target_reads      // val: target read count
    
    main:
    // Quality control and trimming
    QUALITY_CONTROL_RNA(sample_channel)
    
    // Alignment
    ALIGN_RNA(QUALITY_CONTROL_RNA.out.trimmed_reads, genome_index, annotation_gtf)
    
    // Gene quantification
    QUANTIFY_GENES(ALIGN_RNA.out.aligned_bam, annotation_gtf)
    
    // Normalize expression data
    NORMALIZE_EXPRESSION(QUANTIFY_GENES.out.gene_counts.collect(), sample_channel)
    
    // Differential expression analysis
    DIFFERENTIAL_EXPRESSION(NORMALIZE_EXPRESSION.out.normalized_counts, sample_channel)
    
    // Calculate immune signatures
    CALCULATE_IMMUNE_SIGNATURES(NORMALIZE_EXPRESSION.out.normalized_counts, sample_channel)
    
    emit:
    expression_features = QUANTIFY_GENES.out.gene_counts
    normalized_expression = NORMALIZE_EXPRESSION.out.normalized_counts
    differential_expression = DIFFERENTIAL_EXPRESSION.out.de_results
    immune_signatures = CALCULATE_IMMUNE_SIGNATURES.out.immune_signatures
    qc_reports = QUALITY_CONTROL_RNA.out.qc_report
}
