#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Simple WES test process
process WES_TEST {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path(wes_tumor_file), path(wes_normal_file)
    
    output:
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_wes_test.txt"), emit: test_output
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_somatic_variants.vcf"), emit: somatic_variants
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_tmb.txt"), emit: tmb_result
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_cnv.tsv"), emit: cnv_results
    tuple val(sample_id), val(cancer_type), val(stage), val(sample_type), path("${sample_id}_variant_features.tsv"), emit: variant_features
    
    script:
    """
    echo "Processing WES data for ${sample_id}" > ${sample_id}_wes_test.txt
    echo "Sample type: ${sample_type}" >> ${sample_id}_wes_test.txt
    echo "Cancer type: ${cancer_type}" >> ${sample_id}_wes_test.txt
    echo "Stage: ${stage}" >> ${sample_id}_wes_test.txt
    echo "Tumor file: ${wes_tumor_file}" >> ${sample_id}_wes_test.txt
    echo "Normal file: ${wes_normal_file}" >> ${sample_id}_wes_test.txt
    echo "Processing completed at: \$(date)" >> ${sample_id}_wes_test.txt
    
    # Create Python script file
    cat > wes_analysis.py << 'EOF'
import pandas as pd
import numpy as np
import random

# Define key cancer genes
cancer_genes = ['TP53', 'KRAS', 'PIK3CA', 'APC', 'EGFR', 'BRAF', 'MYC', 'TERT']

# Simulate somatic variants based on sample type
variants = []

if '${sample_type}' == 'Cancer':
    # Cancer samples: higher mutation burden
    num_variants = random.randint(15, 35)
    mutation_rate = 0.03
else:
    # Control samples: lower mutation burden
    num_variants = random.randint(2, 8)
    mutation_rate = 0.005

# Generate somatic variants
for i in range(num_variants):
    gene = random.choice(cancer_genes)
    mutation_type = random.choice(['SNV', 'INDEL'])
    
    if mutation_type == 'SNV':
        ref = random.choice(['A', 'T', 'G', 'C'])
        alt = random.choice([b for b in ['A', 'T', 'G', 'C'] if b != ref])
    else:
        ref = random.choice(['AT', 'GC', 'TA', 'CG'])
        alt = random.choice(['A', 'T', 'G', 'C'])
    
    variant = {
        'CHROM': f'chr{random.randint(1, 22)}',
        'POS': random.randint(1000000, 200000000),
        'REF': ref,
        'ALT': alt,
        'GENE': gene,
        'TYPE': mutation_type,
        'QUAL': random.randint(20, 60),
        'FILTER': 'PASS',
        'INFO': f'GENE={gene};TYPE={mutation_type}'
    }
    variants.append(variant)

# Create VCF file
vcf_content = '''##fileformat=VCFv4.2
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Mutation type">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
'''

for variant in variants:
    vcf_content += f"{variant['CHROM']}\\t{variant['POS']}\\t.\\t{variant['REF']}\\t{variant['ALT']}\\t{variant['QUAL']}\\t{variant['FILTER']}\\t{variant['INFO']}\\n"

with open('${sample_id}_somatic_variants.vcf', 'w') as f:
    f.write(vcf_content)

# Calculate TMB (mutations per megabase)
exome_size_mb = 30  # Approximate exome size
tmb = len(variants) / exome_size_mb

with open('${sample_id}_tmb.txt', 'w') as f:
    f.write(f'{tmb:.2f}')

# Simulate CNV analysis
cnv_data = {
    'chromosome': ['chr1', 'chr2', 'chr3', 'chr4', 'chr5'],
    'start': [1000000, 2000000, 3000000, 4000000, 5000000],
    'end': [2000000, 3000000, 4000000, 5000000, 6000000],
    'copy_number': [random.randint(1, 4) for _ in range(5)],
    'log2_ratio': [random.uniform(-1.0, 1.0) for _ in range(5)],
    'cnv_type': [random.choice(['gain', 'loss', 'normal']) for _ in range(5)]
}

cnv_df = pd.DataFrame(cnv_data)
cnv_df.to_csv('${sample_id}_cnv.tsv', sep='\\t', index=False)

# Generate variant features for ML
key_mutations = {gene: random.random() < 0.3 for gene in cancer_genes}

features = {
    'sample_id': '${sample_id}',
    'cancer_type': '${cancer_type}',
    'stage': '${stage}',
    'sample_type': '${sample_type}',
    'tmb': tmb,
    'tmb_category': 'High' if tmb > 10 else 'Low',
    'total_variants': len(variants),
    'snv_count': sum(1 for v in variants if v['TYPE'] == 'SNV'),
    'indel_count': sum(1 for v in variants if v['TYPE'] == 'INDEL'),
    'cnv_gains': sum(1 for cnv in cnv_data['cnv_type'] if cnv == 'gain'),
    'cnv_losses': sum(1 for cnv in cnv_data['cnv_type'] if cnv == 'loss'),
    'key_mutations_count': sum(key_mutations.values())
}

# Add key mutation features
for gene, mutated in key_mutations.items():
    features[f'{gene}_mutated'] = mutated

# Save features
features_df = pd.DataFrame([features])
features_df.to_csv('${sample_id}_variant_features.tsv', sep='\\t', index=False)

print(f'Generated WES analysis for ${sample_id}')
print(f'Total variants: {len(variants)}')
print(f'TMB: {tmb:.2f} mutations/Mb')
print(f'Key mutations: {sum(key_mutations.values())}')
EOF

    # Run Python script
    python3 wes_analysis.py
    """
}

// Test workflow
workflow {
    // Create test sample channel
    Channel.of(
        ["WES_TEST_001", "CRC", "Early", "Cancer", file("${baseDir}/data/test_data/wes/WES_TEST_001/WES_TEST_001_tumor_R1.fastq.gz"), file("${baseDir}/data/test_data/wes/WES_TEST_001/WES_TEST_001_normal_R1.fastq.gz")],
        ["WES_TEST_002", "CRC", "Early", "Control", file("${baseDir}/data/test_data/wes/WES_TEST_002/WES_TEST_002_tumor_R1.fastq.gz"), file("${baseDir}/data/test_data/wes/WES_TEST_002/WES_TEST_002_normal_R1.fastq.gz")],
        ["WES_TEST_003", "Lung", "Late", "Cancer", file("${baseDir}/data/test_data/wes/WES_TEST_003/WES_TEST_003_tumor_R1.fastq.gz"), file("${baseDir}/data/test_data/wes/WES_TEST_003/WES_TEST_003_normal_R1.fastq.gz")]
    )
    .set { sample_channel }
    
    // Run WES test
    WES_TEST(sample_channel)
    
    // Collect results
    WES_TEST.out.test_output
        .collect()
        .set { test_results }
    
    WES_TEST.out.somatic_variants
        .collect()
        .set { variants_results }
    
    WES_TEST.out.tmb_result
        .collect()
        .set { tmb_results }
    
    WES_TEST.out.cnv_results
        .collect()
        .set { cnv_results }
    
    WES_TEST.out.variant_features
        .collect()
        .set { features_results }
    
    // Print results
    test_results
        .map { file -> 
            """
            Test file: ${file}
            """
        }
        .view()
    
    variants_results
        .map { file -> 
            """
            Somatic variants file: ${file}
            """
        }
        .view()
    
    tmb_results
        .map { file -> 
            """
            TMB file: ${file}
            """
        }
        .view()
    
    features_results
        .map { file -> 
            """
            Variant features file: ${file}
            """
        }
        .view()
}

// Workflow completion
workflow.onComplete {
    log.info """
    ================================================
    CENTAUR WES Test completed successfully!
    ================================================
    """
}

// Workflow error handling
workflow.onError {
    log.error """
    ================================================
    CENTAUR WES Test failed!
    ================================================
    
    Error: ${workflow.errorMessage}
    """
    exit 1
}