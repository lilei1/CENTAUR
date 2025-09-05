#!/usr/bin/env python3

"""
Generate simulated EM-seq test data for CENTAUR pipeline testing
"""

import random
import gzip
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def generate_random_sequence(length):
    """Generate a random DNA sequence"""
    bases = ['A', 'T', 'G', 'C']
    return ''.join(random.choices(bases, k=length))

def add_cpg_sites(sequence, cpg_density=0.1):
    """Add CpG sites to the sequence"""
    seq_list = list(sequence)
    for i in range(len(seq_list) - 1):
        if random.random() < cpg_density:
            if seq_list[i] == 'C' and seq_list[i+1] == 'G':
                continue  # Already a CpG
            elif seq_list[i] == 'G' and seq_list[i+1] == 'C':
                continue  # Already a CpG
            else:
                # Create a CpG site
                seq_list[i] = 'C'
                seq_list[i+1] = 'G'
    return ''.join(seq_list)

def simulate_methylation(sequence, methylation_rate=0.7):
    """Simulate methylation by converting some C to T (bisulfite conversion)"""
    seq_list = list(sequence)
    for i in range(len(seq_list)):
        if seq_list[i] == 'C':
            # Simulate methylation - if methylated, C stays C; if not methylated, C becomes T
            if random.random() > methylation_rate:
                seq_list[i] = 'T'
    return ''.join(seq_list)

def generate_fastq_record(sequence, quality_scores, read_id, read_number):
    """Generate a FASTQ record"""
    quality_string = ''.join([chr(q + 33) for q in quality_scores])
    return f"@{read_id}_{read_number}/1\n{sequence}\n+\n{quality_string}\n"

def create_test_samples():
    """Create test EM-seq samples"""
    
    # Test sample parameters
    samples = [
        {"id": "TEST_001", "type": "Cancer", "methylation_rate": 0.3},
        {"id": "TEST_002", "type": "Control", "methylation_rate": 0.7},
        {"id": "TEST_003", "type": "Cancer", "methylation_rate": 0.4}
    ]
    
    for sample in samples:
        print(f"Generating test data for {sample['id']} ({sample['type']})")
        
        # Create output directory
        os.makedirs(f"data/test_data/methylation/{sample['id']}", exist_ok=True)
        
        # Generate paired-end reads
        read_length = 150
        num_reads = 1000  # Small number for testing
        
        r1_file = f"data/test_data/methylation/{sample['id']}/{sample['id']}_R1.fastq.gz"
        r2_file = f"data/test_data/methylation/{sample['id']}/{sample['id']}_R2.fastq.gz"
        
        with gzip.open(r1_file, 'wt') as r1_out, gzip.open(r2_file, 'wt') as r2_out:
            for i in range(num_reads):
                # Generate random sequence
                sequence = generate_random_sequence(read_length)
                
                # Add CpG sites
                sequence = add_cpg_sites(sequence, cpg_density=0.15)
                
                # Simulate methylation
                methylated_seq = simulate_methylation(sequence, sample['methylation_rate'])
                
                # Generate quality scores (simulate good quality)
                quality_scores = [random.randint(30, 40) for _ in range(read_length)]
                
                # Generate reverse complement for R2
                r2_sequence = str(Seq(methylated_seq).reverse_complement())
                r2_quality = quality_scores[::-1]  # Reverse quality scores
                
                # Write R1
                r1_record = generate_fastq_record(methylated_seq, quality_scores, sample['id'], i+1)
                r1_out.write(r1_record)
                
                # Write R2
                r2_record = generate_fastq_record(r2_sequence, r2_quality, sample['id'], i+1)
                r2_out.write(r2_record)
        
        print(f"Generated {num_reads} paired-end reads for {sample['id']}")
        print(f"Files: {r1_file}, {r2_file}")

def create_sample_sheet():
    """Create a test sample sheet"""
    sample_sheet_content = """sample_id,cancer_type,stage,sample_type,cfDNA_file,rnaseq_file,wes_tumor_file,wes_normal_file
TEST_001,CRC,Early,Cancer,TEST_001/TEST_001_R1.fastq.gz,,,
TEST_002,CRC,Early,Control,TEST_002/TEST_002_R1.fastq.gz,,,
TEST_003,Lung,Late,Cancer,TEST_003/TEST_003_R1.fastq.gz,,,"""
    
    with open("data/test_data/sample_sheet_test.csv", "w") as f:
        f.write(sample_sheet_content)
    
    print("Created test sample sheet: data/test_data/sample_sheet_test.csv")

if __name__ == "__main__":
    print("Generating EM-seq test data for CENTAUR pipeline...")
    create_test_samples()
    create_sample_sheet()
    print("Test data generation completed!")
    print("\nTest data structure:")
    print("data/test_data/")
    print("├── methylation/")
    print("│   ├── TEST_001/")
    print("│   │   ├── TEST_001_R1.fastq.gz")
    print("│   │   └── TEST_001_R2.fastq.gz")
    print("│   ├── TEST_002/")
    print("│   │   ├── TEST_002_R1.fastq.gz")
    print("│   │   └── TEST_002_R2.fastq.gz")
    print("│   └── TEST_003/")
    print("│       ├── TEST_003_R1.fastq.gz")
    print("│       └── TEST_003_R2.fastq.gz")
    print("└── sample_sheet_test.csv")
