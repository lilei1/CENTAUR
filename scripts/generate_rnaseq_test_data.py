#!/usr/bin/env python3

"""
Generate simulated RNA-seq test data for CENTAUR pipeline testing
"""

import random
import gzip
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def generate_random_rna_sequence(length):
    """Generate a random RNA sequence"""
    bases = ['A', 'U', 'G', 'C']
    return ''.join(random.choices(bases, k=length))

def simulate_gene_expression(sequence, expression_level='medium'):
    """Simulate gene expression by adjusting sequence complexity"""
    if expression_level == 'high':
        # High expression: more complex sequence with higher GC content
        seq_list = list(sequence)
        for i in range(len(seq_list)):
            if random.random() < 0.3:  # Increase GC content
                if seq_list[i] == 'A':
                    seq_list[i] = random.choice(['G', 'C'])
                elif seq_list[i] == 'U':
                    seq_list[i] = random.choice(['G', 'C'])
        return ''.join(seq_list)
    elif expression_level == 'low':
        # Low expression: simpler sequence with lower GC content
        seq_list = list(sequence)
        for i in range(len(seq_list)):
            if random.random() < 0.2:  # Decrease GC content
                if seq_list[i] == 'G':
                    seq_list[i] = random.choice(['A', 'U'])
                elif seq_list[i] == 'C':
                    seq_list[i] = random.choice(['A', 'U'])
        return ''.join(seq_list)
    else:
        # Medium expression: keep original
        return sequence

def generate_fastq_record(sequence, quality_scores, read_id, read_number):
    """Generate a FASTQ record"""
    quality_string = ''.join([chr(q + 33) for q in quality_scores])
    return f"@{read_id}_{read_number}/1\n{sequence}\n+\n{quality_string}\n"

def create_rnaseq_samples():
    """Create test RNA-seq samples"""
    
    # Define gene expression patterns for different sample types
    gene_expression_patterns = {
        'Cancer': {
            'immune_genes': 'high',      # High immune gene expression
            'cell_cycle_genes': 'high',  # High cell cycle genes
            'apoptosis_genes': 'low',    # Low apoptosis genes
            'housekeeping_genes': 'medium' # Medium housekeeping genes
        },
        'Control': {
            'immune_genes': 'medium',    # Medium immune gene expression
            'cell_cycle_genes': 'medium', # Medium cell cycle genes
            'apoptosis_genes': 'medium', # Medium apoptosis genes
            'housekeeping_genes': 'medium' # Medium housekeeping genes
        }
    }
    
    # Test sample parameters
    samples = [
        {"id": "RNA_TEST_001", "type": "Cancer", "cancer_type": "CRC", "stage": "Early"},
        {"id": "RNA_TEST_002", "type": "Control", "cancer_type": "CRC", "stage": "Early"},
        {"id": "RNA_TEST_003", "type": "Cancer", "cancer_type": "Lung", "stage": "Late"}
    ]
    
    for sample in samples:
        print(f"Generating RNA-seq test data for {sample['id']} ({sample['type']})")
        
        # Create output directory
        os.makedirs(f"data/test_data/rnaseq/{sample['id']}", exist_ok=True)
        
        # Generate paired-end reads
        read_length = 150
        num_reads = 2000  # More reads for RNA-seq
        
        r1_file = f"data/test_data/rnaseq/{sample['id']}/{sample['id']}_R1.fastq.gz"
        r2_file = f"data/test_data/rnaseq/{sample['id']}/{sample['id']}_R2.fastq.gz"
        
        # Get expression pattern for this sample type
        expression_pattern = gene_expression_patterns[sample['type']]
        
        with gzip.open(r1_file, 'wt') as r1_out, gzip.open(r2_file, 'wt') as r2_out:
            for i in range(num_reads):
                # Generate random RNA sequence
                sequence = generate_random_rna_sequence(read_length)
                
                # Simulate gene expression based on sample type
                # Randomly assign gene type for this read
                gene_types = ['immune_genes', 'cell_cycle_genes', 'apoptosis_genes', 'housekeeping_genes']
                gene_type = random.choice(gene_types)
                expression_level = expression_pattern[gene_type]
                
                # Apply expression level to sequence
                expressed_seq = simulate_gene_expression(sequence, expression_level)
                
                # Generate quality scores (simulate good quality)
                quality_scores = [random.randint(30, 40) for _ in range(read_length)]
                
                # Generate reverse complement for R2
                r2_sequence = str(Seq(expressed_seq).reverse_complement())
                r2_quality = quality_scores[::-1]  # Reverse quality scores
                
                # Write R1
                r1_record = generate_fastq_record(expressed_seq, quality_scores, sample['id'], i+1)
                r1_out.write(r1_record)
                
                # Write R2
                r2_record = generate_fastq_record(r2_sequence, r2_quality, sample['id'], i+1)
                r2_out.write(r2_record)
        
        print(f"Generated {num_reads} paired-end reads for {sample['id']}")
        print(f"Files: {r1_file}, {r2_file}")

def create_rnaseq_sample_sheet():
    """Create a test sample sheet for RNA-seq"""
    sample_sheet_content = """sample_id,cancer_type,stage,sample_type,cfDNA_file,rnaseq_file,wes_tumor_file,wes_normal_file
RNA_TEST_001,CRC,Early,Cancer,,RNA_TEST_001/RNA_TEST_001_R1.fastq.gz,,
RNA_TEST_002,CRC,Early,Control,,RNA_TEST_002/RNA_TEST_002_R1.fastq.gz,,
RNA_TEST_003,Lung,Late,Cancer,,RNA_TEST_003/RNA_TEST_003_R1.fastq.gz,,"""
    
    with open("data/test_data/rnaseq_sample_sheet_test.csv", "w") as f:
        f.write(sample_sheet_content)
    
    print("Created RNA-seq test sample sheet: data/test_data/rnaseq_sample_sheet_test.csv")

def create_mock_annotation():
    """Create a mock gene annotation file"""
    annotation_content = """##gtf-version 2.2
chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript";
chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "1"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002234944"; exon_version "1";
chr1	HAVANA	exon	12613	12721	.	+	.	gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00003582793"; exon_version "1";
chr1	HAVANA	exon	13221	14409	.	+	.	gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "3"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002312635"; exon_version "1";
chr1	HAVANA	gene	14362	29370	.	-	.	gene_id "ENSG00000227232"; gene_version "5"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene";
chr1	HAVANA	transcript	14362	29370	.	-	.	gene_id "ENSG00000227232"; gene_version "5"; transcript_id "ENST00000438513"; transcript_version "1"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene"; transcript_name "WASH7P-201"; transcript_source "havana"; transcript_biotype "unprocessed_pseudogene";
chr1	HAVANA	exon	14362	14829	.	-	.	gene_id "ENSG00000227232"; gene_version "5"; transcript_id "ENST00000438513"; transcript_version "1"; exon_number "1"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene"; transcript_name "WASH7P-201"; transcript_source "havana"; transcript_biotype "unprocessed_pseudogene"; exon_id "ENSE00003623524"; exon_version "1";
chr1	HAVANA	exon	14970	15038	.	-	.	gene_id "ENSG00000227232"; gene_version "5"; transcript_id "ENST00000438513"; transcript_version "1"; exon_number "2"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene"; transcript_name "WASH7P-201"; transcript_source "havana"; transcript_biotype "unprocessed_pseudogene"; exon_id "ENSE00002154430"; exon_version "1";
chr1	HAVANA	exon	15796	15947	.	-	.	gene_id "ENSG00000227232"; gene_version "5"; transcript_id "ENST00000438513"; transcript_version "1"; exon_number "3"; gene_name "WASH7P"; gene_source "havana"; gene_biotype "unprocessed_pseudogene"; transcript_name "WASH7P-201"; transcript_source "havana"; transcript_biotype "unprocessed_pseudogene"; exon_id "ENSE00002001255"; exon_version "1";"""
    
    with open("data/test_data/mock_annotation.gtf", "w") as f:
        f.write(annotation_content)
    
    print("Created mock gene annotation: data/test_data/mock_annotation.gtf")

if __name__ == "__main__":
    print("Generating RNA-seq test data for CENTAUR pipeline...")
    create_rnaseq_samples()
    create_rnaseq_sample_sheet()
    create_mock_annotation()
    print("RNA-seq test data generation completed!")
    print("\nTest data structure:")
    print("data/test_data/")
    print("├── rnaseq/")
    print("│   ├── RNA_TEST_001/")
    print("│   │   ├── RNA_TEST_001_R1.fastq.gz")
    print("│   │   └── RNA_TEST_001_R2.fastq.gz")
    print("│   ├── RNA_TEST_002/")
    print("│   │   ├── RNA_TEST_002_R1.fastq.gz")
    print("│   │   └── RNA_TEST_002_R2.fastq.gz")
    print("│   └── RNA_TEST_003/")
    print("│       ├── RNA_TEST_003_R1.fastq.gz")
    print("│       └── RNA_TEST_003_R2.fastq.gz")
    print("├── rnaseq_sample_sheet_test.csv")
    print("└── mock_annotation.gtf")
