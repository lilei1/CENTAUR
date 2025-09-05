#!/usr/bin/env python3

"""
Generate simulated WES test data for CENTAUR pipeline testing
"""

import random
import gzip
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def generate_random_dna_sequence(length):
    """Generate a random DNA sequence"""
    bases = ['A', 'T', 'G', 'C']
    return ''.join(random.choices(bases, k=length))

def simulate_somatic_mutations(sequence, mutation_rate=0.01):
    """Simulate somatic mutations in tumor sequence"""
    seq_list = list(sequence)
    mutations = []
    
    # Process mutations in reverse order to avoid index issues
    positions_to_mutate = []
    for i in range(len(seq_list)):
        if random.random() < mutation_rate:
            positions_to_mutate.append(i)
    
    # Sort in reverse order to process from end to beginning
    positions_to_mutate.sort(reverse=True)
    
    for i in positions_to_mutate:
        if i < len(seq_list):
            original_base = seq_list[i]
            # Simulate different mutation types
            mutation_type = random.choice(['SNV', 'INDEL'])
            
            if mutation_type == 'SNV':
                # Single nucleotide variant
                new_base = random.choice([b for b in ['A', 'T', 'G', 'C'] if b != original_base])
                seq_list[i] = new_base
                mutations.append({
                    'pos': i,
                    'ref': original_base,
                    'alt': new_base,
                    'type': 'SNV'
                })
            else:
                # Small indel
                if random.random() < 0.5:
                    # Insertion
                    insertion = random.choice(['A', 'T', 'G', 'C'])
                    seq_list.insert(i, insertion)
                    mutations.append({
                        'pos': i,
                        'ref': original_base,
                        'alt': original_base + insertion,
                        'type': 'INS'
                    })
                else:
                    # Deletion
                    if i < len(seq_list) - 1:
                        deleted_base = seq_list[i+1]
                        seq_list.pop(i+1)
                        mutations.append({
                            'pos': i,
                            'ref': original_base + deleted_base,
                            'alt': original_base,
                            'type': 'DEL'
                        })
    
    return ''.join(seq_list), mutations

def generate_fastq_record(sequence, quality_scores, read_id, read_number):
    """Generate a FASTQ record"""
    quality_string = ''.join([chr(q + 33) for q in quality_scores])
    return f"@{read_id}_{read_number}/1\n{sequence}\n+\n{quality_string}\n"

def create_wes_samples():
    """Create test WES samples"""
    
    # Test sample parameters
    samples = [
        {"id": "WES_TEST_001", "type": "Cancer", "cancer_type": "CRC", "stage": "Early"},
        {"id": "WES_TEST_002", "type": "Control", "cancer_type": "CRC", "stage": "Early"},
        {"id": "WES_TEST_003", "type": "Cancer", "cancer_type": "Lung", "stage": "Late"}
    ]
    
    for sample in samples:
        print(f"Generating WES test data for {sample['id']} ({sample['type']})")
        
        # Create output directory
        os.makedirs(f"data/test_data/wes/{sample['id']}", exist_ok=True)
        
        # Generate paired-end reads
        read_length = 150
        num_reads = 1500  # WES typically has fewer reads than RNA-seq
        
        # Generate tumor and normal samples
        for sample_type in ['tumor', 'normal']:
            r1_file = f"data/test_data/wes/{sample['id']}/{sample['id']}_{sample_type}_R1.fastq.gz"
            r2_file = f"data/test_data/wes/{sample['id']}/{sample['id']}_{sample_type}_R2.fastq.gz"
            
            with gzip.open(r1_file, 'wt') as r1_out, gzip.open(r2_file, 'wt') as r2_out:
                for i in range(num_reads):
                    # Generate random DNA sequence
                    sequence = generate_random_dna_sequence(read_length)
                    
                    # Simulate mutations only in tumor samples
                    if sample_type == 'tumor' and sample['type'] == 'Cancer':
                        # Higher mutation rate for cancer samples
                        mutation_rate = random.uniform(0.02, 0.05)
                        mutated_seq, mutations = simulate_somatic_mutations(sequence, mutation_rate)
                        final_seq = mutated_seq
                    else:
                        # Normal samples or control samples have fewer/no mutations
                        mutation_rate = random.uniform(0.001, 0.005)
                        mutated_seq, mutations = simulate_somatic_mutations(sequence, mutation_rate)
                        final_seq = mutated_seq
                    
                    # Generate quality scores (simulate good quality)
                    quality_scores = [random.randint(30, 40) for _ in range(len(final_seq))]
                    
                    # Generate reverse complement for R2
                    r2_sequence = str(Seq(final_seq).reverse_complement())
                    r2_quality = quality_scores[::-1]  # Reverse quality scores
                    
                    # Write R1
                    r1_record = generate_fastq_record(final_seq, quality_scores, f"{sample['id']}_{sample_type}", i+1)
                    r1_out.write(r1_record)
                    
                    # Write R2
                    r2_record = generate_fastq_record(r2_sequence, r2_quality, f"{sample['id']}_{sample_type}", i+1)
                    r2_out.write(r2_record)
            
            print(f"Generated {num_reads} paired-end reads for {sample['id']} {sample_type}")
            print(f"Files: {r1_file}, {r2_file}")

def create_wes_sample_sheet():
    """Create a test sample sheet for WES"""
    sample_sheet_content = """sample_id,cancer_type,stage,sample_type,cfDNA_file,rnaseq_file,wes_tumor_file,wes_normal_file
WES_TEST_001,CRC,Early,Cancer,,,WES_TEST_001/WES_TEST_001_tumor_R1.fastq.gz,WES_TEST_001/WES_TEST_001_normal_R1.fastq.gz
WES_TEST_002,CRC,Early,Control,,,WES_TEST_002/WES_TEST_002_tumor_R1.fastq.gz,WES_TEST_002/WES_TEST_002_normal_R1.fastq.gz
WES_TEST_003,Lung,Late,Cancer,,,WES_TEST_003/WES_TEST_003_tumor_R1.fastq.gz,WES_TEST_003/WES_TEST_003_normal_R1.fastq.gz"""
    
    with open("data/test_data/wes_sample_sheet_test.csv", "w") as f:
        f.write(sample_sheet_content)
    
    print("Created WES test sample sheet: data/test_data/wes_sample_sheet_test.csv")

def create_mock_target_intervals():
    """Create mock target intervals for WES"""
    target_content = """chr1	10000	20000	exon1
chr1	30000	40000	exon2
chr1	50000	60000	exon3
chr2	10000	20000	exon4
chr2	30000	40000	exon5
chr3	10000	20000	exon6
chr3	30000	40000	exon7
chr4	10000	20000	exon8
chr4	30000	40000	exon9
chr5	10000	20000	exon10"""
    
    with open("data/test_data/mock_target_intervals.bed", "w") as f:
        f.write(target_content)
    
    print("Created mock target intervals: data/test_data/mock_target_intervals.bed")

if __name__ == "__main__":
    print("Generating WES test data for CENTAUR pipeline...")
    create_wes_samples()
    create_wes_sample_sheet()
    create_mock_target_intervals()
    print("WES test data generation completed!")
    print("\nTest data structure:")
    print("data/test_data/")
    print("├── wes/")
    print("│   ├── WES_TEST_001/")
    print("│   │   ├── WES_TEST_001_tumor_R1.fastq.gz")
    print("│   │   ├── WES_TEST_001_tumor_R2.fastq.gz")
    print("│   │   ├── WES_TEST_001_normal_R1.fastq.gz")
    print("│   │   └── WES_TEST_001_normal_R2.fastq.gz")
    print("│   ├── WES_TEST_002/")
    print("│   │   ├── WES_TEST_002_tumor_R1.fastq.gz")
    print("│   │   ├── WES_TEST_002_tumor_R2.fastq.gz")
    print("│   │   ├── WES_TEST_002_normal_R1.fastq.gz")
    print("│   │   └── WES_TEST_002_normal_R2.fastq.gz")
    print("│   └── WES_TEST_003/")
    print("│       ├── WES_TEST_003_tumor_R1.fastq.gz")
    print("│       ├── WES_TEST_003_tumor_R2.fastq.gz")
    print("│       ├── WES_TEST_003_normal_R1.fastq.gz")
    print("│       └── WES_TEST_003_normal_R2.fastq.gz")
    print("├── wes_sample_sheet_test.csv")
    print("└── mock_target_intervals.bed")
