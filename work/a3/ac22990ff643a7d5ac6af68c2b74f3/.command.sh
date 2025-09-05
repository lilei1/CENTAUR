#!/bin/bash -ue
echo "Processing RNA-seq data for RNA_TEST_001" > RNA_TEST_001_rnaseq_test.txt
echo "Sample type: Cancer" >> RNA_TEST_001_rnaseq_test.txt
echo "Cancer type: CRC" >> RNA_TEST_001_rnaseq_test.txt
echo "Stage: Early" >> RNA_TEST_001_rnaseq_test.txt
echo "Input file: RNA_TEST_001_R1.fastq.gz" >> RNA_TEST_001_rnaseq_test.txt
echo "Processing completed at: $(date)" >> RNA_TEST_001_rnaseq_test.txt

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
    if 'Cancer' == 'Cancer':
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
counts_df.to_csv('RNA_TEST_001_gene_counts.txt', sep='	', index=False)

# Calculate immune signatures
immune_signatures = {
    'sample_id': 'RNA_TEST_001',
    'IFN_gamma_response': np.mean([gene_counts[g] for g in ['STAT1', 'IRF1', 'CXCL10', 'CXCL9', 'IDO1']]),
    'T_cell_exhaustion': np.mean([gene_counts[g] for g in ['PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT']]),
    'Cytotoxic_T_cells': np.mean([gene_counts[g] for g in ['GZMA', 'GZMB', 'PRF1', 'GNLY', 'NKG7']]),
    'B_cells': np.mean([gene_counts[g] for g in ['CD19', 'CD20', 'CD22', 'CD79A', 'CD79B']]),
    'NK_cells': np.mean([gene_counts[g] for g in ['KIR2DL1', 'KIR2DL2', 'KIR2DL3', 'KIR2DL4']]),
    'Housekeeping': np.mean([gene_counts[g] for g in ['GAPDH', 'ACTB', 'TUBB', 'RPL13A', 'RPS18']])
}

# Create immune signatures DataFrame
immune_df = pd.DataFrame([immune_signatures])
immune_df.to_csv('RNA_TEST_001_immune_signatures.tsv', sep='	', index=False)

print(f'Generated gene expression data for RNA_TEST_001')
print(f'Total genes: {len(genes)}')
print(f'Immune signatures calculated')
"
