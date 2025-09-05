#!/usr/bin/env nextflow

/*
 * WES Processing Workflow
 * Processes whole exome sequencing data for variant calling, TMB calculation, and CNV analysis
 */

nextflow.enable.dsl = 2

// Process definitions
process PREPARE_WES_DATA {
    tag "${sample_id}"
    label 'process_low'
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(wes_tumor_file), path(wes_normal_file)
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_tumor.bam"), path("${sample_id}_normal.bam"), emit: wes_bams
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_wes_info.txt"), emit: wes_info
    
    script:
    """
    # Copy and rename files for consistency
    cp ${wes_tumor_file} ${sample_id}_tumor.bam
    cp ${wes_normal_file} ${sample_id}_normal.bam
    
    # Index BAM files
    samtools index ${sample_id}_tumor.bam
    samtools index ${sample_id}_normal.bam
    
    # Generate WES info file
    echo "Sample ID: ${sample_id}" > ${sample_id}_wes_info.txt
    echo "Cancer Type: ${cancer_type}" >> ${sample_id}_wes_info.txt
    echo "Stage: ${stage}" >> ${sample_id}_wes_info.txt
    echo "Sample Type: ${sample_type}" >> ${sample_id}_wes_info.txt
    echo "Tumor BAM: ${sample_id}_tumor.bam" >> ${sample_id}_wes_info.txt
    echo "Normal BAM: ${sample_id}_normal.bam" >> ${sample_id}_wes_info.txt
    """
}

process CALL_SOMATIC_VARIANTS {
    tag "${sample_id}"
    label 'process_high'
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(tumor_bam), path(normal_bam)
    path genome_fasta
    path target_intervals
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_somatic.vcf"), emit: somatic_vcf
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_variant_stats.txt"), emit: variant_stats
    
    script:
    """
    # Call somatic variants with Mutect2
    gatk Mutect2 \
        -R ${genome_fasta} \
        -I ${tumor_bam} \
        -I ${normal_bam} \
        -normal $(basename ${normal_bam} .bam) \
        -L ${target_intervals} \
        -O ${sample_id}_somatic.vcf \
        --germline-resource af-only-gnomad.hg38.vcf.gz \
        --af-of-alleles-not-in-resource 0.0000025 \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter
    
    # Filter variants
    gatk FilterMutectCalls \
        -V ${sample_id}_somatic.vcf \
        -O ${sample_id}_somatic_filtered.vcf \
        -R ${genome_fasta}
    
    # Generate variant statistics
    echo "Somatic variant calling completed for ${sample_id}" > ${sample_id}_variant_stats.txt
    echo "Total variants: $(grep -v '^#' ${sample_id}_somatic_filtered.vcf | wc -l)" >> ${sample_id}_variant_stats.txt
    echo "SNVs: $(grep -v '^#' ${sample_id}_somatic_filtered.vcf | grep -v 'SVTYPE=INS' | grep -v 'SVTYPE=DEL' | wc -l)" >> ${sample_id}_variant_stats.txt
    echo "Indels: $(grep -v '^#' ${sample_id}_somatic_filtered.vcf | grep -E 'SVTYPE=INS|SVTYPE=DEL' | wc -l)" >> ${sample_id}_variant_stats.txt
    """
}

process ANNOTATE_VARIANTS {
    tag "${sample_id}"
    label 'process_medium'
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(somatic_vcf)
    path annotation_db
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_annotated.vcf"), emit: annotated_vcf
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_annotation_summary.txt"), emit: annotation_summary
    
    script:
    """
    # Annotate variants with ANNOVAR
    annovar \
        ${somatic_vcf} \
        ${annotation_db} \
        -buildver hg38 \
        -out ${sample_id}_annotated \
        -protocol refGene,cosmic70,clinvar_20170905,avsnp147,dbnsfp30a \
        -operation g,f,f,f,f \
        -nastring . \
        -vcfinput
    
    # Generate annotation summary
    echo "Variant annotation completed for ${sample_id}" > ${sample_id}_annotation_summary.txt
    echo "Annotated variants: $(grep -v '^#' ${sample_id}_annotated.hg38_multianno.vcf | wc -l)" >> ${sample_id}_annotation_summary.txt
    """
}

process EXTRACT_MUTATIONAL_SIGNATURES {
    tag "${sample_id}"
    label 'process_medium'
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(annotated_vcf)
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_signatures.tsv"), emit: mutational_signatures
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_signature_plot.png"), emit: signature_plot
    
    script:
    """
    # Extract mutational signatures using SigProfilerExtractor
    python3 -c "
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # Read VCF file and extract mutations
    mutations = []
    
    with open('${annotated_vcf}', 'r') as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.strip().split('\\t')
                if len(fields) >= 5:
                    chrom = fields[0]
                    pos = int(fields[1])
                    ref = fields[3]
                    alt = fields[4]
                    
                    # Only process SNVs
                    if len(ref) == 1 and len(alt) == 1:
                        mutations.append({
                            'chromosome': chrom,
                            'position': pos,
                            'reference': ref,
                            'alternate': alt
                        })
    
    # Create mutation context (simplified)
    # In practice, you would use SigProfilerExtractor
    mutation_contexts = {
        'C>A': 0.15,
        'C>G': 0.05,
        'C>T': 0.35,
        'T>A': 0.10,
        'T>C': 0.20,
        'T>G': 0.15
    }
    
    # Create signature exposure (simplified)
    signatures = {
        'Signature_1': 0.3,
        'Signature_2': 0.2,
        'Signature_3': 0.15,
        'Signature_4': 0.1,
        'Signature_5': 0.1,
        'Signature_6': 0.15
    }
    
    # Save signature data
    sig_df = pd.DataFrame(list(signatures.items()), columns=['Signature', 'Exposure'])
    sig_df.to_csv('${sample_id}_signatures.tsv', sep='\t', index=False)
    
    # Create signature plot
    plt.figure(figsize=(12, 6))
    
    plt.subplot(1, 2, 1)
    plt.bar(mutation_contexts.keys(), mutation_contexts.values())
    plt.title('Mutation Context Distribution')
    plt.xlabel('Mutation Type')
    plt.ylabel('Frequency')
    plt.xticks(rotation=45)
    
    plt.subplot(1, 2, 2)
    plt.bar(sig_df['Signature'], sig_df['Exposure'])
    plt.title('Mutational Signature Exposure')
    plt.xlabel('Signature')
    plt.ylabel('Exposure')
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig('${sample_id}_signature_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f'Processed {len(mutations)} mutations for ${sample_id}')
    "
    """
}

process CALCULATE_TMB {
    tag "${sample_id}"
    label 'process_low'
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(annotated_vcf)
    path target_intervals
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_tmb.txt"), emit: tmb_result
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_tmb_summary.txt"), emit: tmb_summary
    
    script:
    """
    # Calculate Tumor Mutational Burden (TMB)
    python3 -c "
    import pandas as pd
    import numpy as np
    
    # Read VCF file
    mutations = []
    
    with open('${annotated_vcf}', 'r') as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.strip().split('\\t')
                if len(fields) >= 5:
                    chrom = fields[0]
                    pos = int(fields[1])
                    ref = fields[3]
                    alt = fields[4]
                    
                    # Only count SNVs and indels
                    if len(ref) <= 2 and len(alt) <= 2:
                        mutations.append({
                            'chromosome': chrom,
                            'position': pos,
                            'reference': ref,
                            'alternate': alt
                        })
    
    # Calculate TMB (mutations per megabase)
    # Assuming exome size of ~30 Mb
    exome_size_mb = 30
    tmb = len(mutations) / exome_size_mb
    
    # Save TMB result
    with open('${sample_id}_tmb.txt', 'w') as f:
        f.write(f'{tmb:.2f}')
    
    # Generate TMB summary
    with open('${sample_id}_tmb_summary.txt', 'w') as f:
        f.write(f'Sample: ${sample_id}\\n')
        f.write(f'Total mutations: {len(mutations)}\\n')
        f.write(f'Exome size (Mb): {exome_size_mb}\\n')
        f.write(f'TMB (mutations/Mb): {tmb:.2f}\\n')
        f.write(f'TMB category: {\"High\" if tmb > 10 else \"Low\"}\\n')
    
    print(f'TMB calculated for ${sample_id}: {tmb:.2f} mutations/Mb')
    "
    """
}

process ANALYZE_CNV {
    tag "${sample_id}"
    label 'process_medium'
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(tumor_bam), path(normal_bam)
    path genome_fasta
    path target_intervals
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_cnv.tsv"), emit: cnv_results
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_cnv_plot.png"), emit: cnv_plot
    
    script:
    """
    # Analyze copy number variations using cn.mops
    python3 -c "
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    
    # This is a simplified CNV analysis
    # In practice, you would use specialized tools like cn.mops, Sequenza, or FACETS
    
    # Create mock CNV data for demonstration
    cnv_data = {
        'chromosome': ['chr1', 'chr2', 'chr3', 'chr4', 'chr5'],
        'start': [1000000, 2000000, 3000000, 4000000, 5000000],
        'end': [2000000, 3000000, 4000000, 5000000, 6000000],
        'copy_number': [3, 1, 4, 2, 3],
        'log2_ratio': [0.58, -1.0, 1.0, 0.0, 0.58],
        'cnv_type': ['gain', 'loss', 'gain', 'normal', 'gain']
    }
    
    cnv_df = pd.DataFrame(cnv_data)
    cnv_df.to_csv('${sample_id}_cnv.tsv', sep='\t', index=False)
    
    # Create CNV plot
    plt.figure(figsize=(12, 6))
    
    plt.subplot(1, 2, 1)
    plt.bar(range(len(cnv_df)), cnv_df['copy_number'])
    plt.xlabel('CNV Region')
    plt.ylabel('Copy Number')
    plt.title('Copy Number Variations')
    plt.xticks(range(len(cnv_df)), cnv_df['chromosome'], rotation=45)
    
    plt.subplot(1, 2, 2)
    plt.bar(range(len(cnv_df)), cnv_df['log2_ratio'])
    plt.xlabel('CNV Region')
    plt.ylabel('Log2 Ratio')
    plt.title('Log2 Copy Number Ratio')
    plt.xticks(range(len(cnv_df)), cnv_df['chromosome'], rotation=45)
    plt.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('${sample_id}_cnv_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f'CNV analysis completed for ${sample_id}')
    "
    """
}

process GENERATE_VARIANT_FEATURES {
    tag "${sample_id}"
    label 'process_low'
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(annotated_vcf), path(tmb_result), path(cnv_results)
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_variant_features.tsv"), emit: variant_features
    
    script:
    """
    # Generate variant features for ML
    python3 -c "
    import pandas as pd
    import numpy as np
    
    # Read TMB
    with open('${tmb_result}', 'r') as f:
        tmb = float(f.read().strip())
    
    # Read CNV data
    cnv_df = pd.read_csv('${cnv_results}', sep='\t')
    
    # Read VCF and extract key mutations
    key_mutations = {
        'TP53': False,
        'KRAS': False,
        'PIK3CA': False,
        'APC': False,
        'EGFR': False,
        'BRAF': False
    }
    
    with open('${annotated_vcf}', 'r') as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.strip().split('\\t')
                if len(fields) >= 5:
                    # Extract gene name from annotation (simplified)
                    info = fields[7]
                    for gene in key_mutations.keys():
                        if gene in info:
                            key_mutations[gene] = True
    
    # Calculate CNV features
    cnv_gains = sum(cnv_df['cnv_type'] == 'gain')
    cnv_losses = sum(cnv_df['cnv_type'] == 'loss')
    total_cnv = len(cnv_df)
    
    # Create feature vector
    features = {
        'sample_id': '${sample_id}',
        'cancer_type': '${cancer_type}',
        'stage': '${stage}',
        'sample_type': '${sample_type}',
        'tmb': tmb,
        'tmb_category': 'High' if tmb > 10 else 'Low',
        'total_cnv': total_cnv,
        'cnv_gains': cnv_gains,
        'cnv_losses': cnv_losses,
        'cnv_gain_ratio': cnv_gains / total_cnv if total_cnv > 0 else 0,
        'cnv_loss_ratio': cnv_losses / total_cnv if total_cnv > 0 else 0
    }
    
    # Add key mutation features
    for gene, mutated in key_mutations.items():
        features[f'{gene}_mutated'] = mutated
    
    # Add mutation count features
    features['key_mutations_count'] = sum(key_mutations.values())
    features['key_mutations_ratio'] = sum(key_mutations.values()) / len(key_mutations)
    
    # Save features
    df = pd.DataFrame([features])
    df.to_csv('${sample_id}_variant_features.tsv', sep='\t', index=False)
    
    print(f'Variant features generated for ${sample_id}')
    "
    """
}

// Workflow definition
workflow WES_PROCESSING {
    take:
    sample_channel    // channel: [sample_id, cancer_type, stage, sample_type, wes_tumor_file, wes_normal_file]
    genome_fasta      // path: reference genome
    genome_index      // path: genome index
    target_coverage   // val: target coverage
    
    main:
    // Prepare WES data
    PREPARE_WES_DATA(sample_channel)
    
    // Call somatic variants
    CALL_SOMATIC_VARIANTS(PREPARE_WES_DATA.out.wes_bams, genome_fasta, target_intervals)
    
    // Annotate variants
    ANNOTATE_VARIANTS(CALL_SOMATIC_VARIANTS.out.somatic_vcf, annotation_db)
    
    // Extract mutational signatures
    EXTRACT_MUTATIONAL_SIGNATURES(ANNOTATE_VARIANTS.out.annotated_vcf)
    
    // Calculate TMB
    CALCULATE_TMB(ANNOTATE_VARIANTS.out.annotated_vcf, target_intervals)
    
    // Analyze CNV
    ANALYZE_CNV(PREPARE_WES_DATA.out.wes_bams, genome_fasta, target_intervals)
    
    // Generate variant features
    GENERATE_VARIANT_FEATURES(
        ANNOTATE_VARIANTS.out.annotated_vcf,
        CALCULATE_TMB.out.tmb_result,
        ANALYZE_CNV.out.cnv_results
    )
    
    emit:
    variant_features = GENERATE_VARIANT_FEATURES.out.variant_features
    tmb_features = CALCULATE_TMB.out.tmb_result
    cnv_features = ANALYZE_CNV.out.cnv_results
    mutational_signatures = EXTRACT_MUTATIONAL_SIGNATURES.out.mutational_signatures
    somatic_variants = CALL_SOMATIC_VARIANTS.out.somatic_vcf
}
