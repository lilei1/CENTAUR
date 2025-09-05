# CENTAUR Multi-Omics Pipeline

A Scalable Platform for Early Cancer Biomarker Discovery

## Overview

CENTAUR is a comprehensive Nextflow-based pipeline for processing multi-omics data to identify early cancer biomarkers. The pipeline integrates cfDNA methylation, RNA-seq, and whole exome sequencing (WES) data to generate ML-ready features for biomarker discovery.

## Features

- **Multi-omics Integration**: Processes cfDNA methylation, RNA-seq, and WES data
- **Containerized Execution**: Docker-based execution for reproducibility
- **Scalable Architecture**: Supports local, HPC, and cloud execution
- **ML-ready Output**: Generates feature matrices for downstream machine learning
- **Comprehensive QC**: Quality control and reporting at each step

## Pipeline Components

### 1. Methylation Analysis (cfDNA EM-seq)
- Quality control and read trimming
- Alignment with Bismark
- Methylation calling with MethylDackel
- Fragmentomics feature extraction
- DMR (Differentially Methylated Regions) calling

### 2. RNA-seq Analysis (PBMC)
- Quality control and read trimming
- Alignment with STAR
- Gene quantification with featureCounts
- Expression normalization and differential analysis
- Immune signature calculation

### 3. WES Analysis (Tumor + Normal)
- Somatic variant calling with Mutect2
- Variant annotation with ANNOVAR
- TMB (Tumor Mutational Burden) calculation
- CNV (Copy Number Variation) analysis
- Mutational signature extraction

### 4. Multi-omics Integration
- Sample harmonization across modalities
- Biomarker feature extraction and selection
- Dimensionality reduction (PCA, t-SNE)
- ML-ready feature matrix generation

## Installation

### Prerequisites
- Docker
- Nextflow
- Python 3.8+
- R 4.0+

### Setup
```bash
# Clone the repository
git clone <repository-url>
cd CENTAUR

# Run setup script
./setup.sh
```

## Usage

### Basic Usage
```bash
nextflow run workflows/main.nf \
    --input_dir data/raw \
    --sample_sheet data/sample_sheet.csv \
    --genome_fasta references/hg38.fa \
    --genome_index references/hg38 \
    --annotation_gtf references/hg38.gtf \
    --output_dir data/results
```

### Configuration
The pipeline supports different execution environments:

- **Local execution**: `nextflow run workflows/main.nf -c configs/local.config`
- **HPC/SLURM**: `nextflow run workflows/main.nf -c configs/hpc.config`
- **AWS Batch**: `nextflow run workflows/main.nf -c configs/aws.config`

### Sample Sheet Format
The sample sheet should be a CSV file with the following columns:
- `sample_id`: Unique sample identifier
- `cancer_type`: Cancer type (CRC, Lung, Breast)
- `stage`: Cancer stage (Early, Late)
- `sample_type`: Sample type (Cancer, Control)
- `cfDNA_file`: Path to cfDNA sequencing file
- `rnaseq_file`: Path to RNA-seq file
- `wes_tumor_file`: Path to WES tumor BAM file
- `wes_normal_file`: Path to WES normal BAM file

## Output Structure

```
data/results/
├── methylation/
│   ├── methylation_profiles/
│   ├── fragmentomics_features/
│   └── dmr_results/
├── rnaseq/
│   ├── expression_counts/
│   ├── differential_expression/
│   └── immune_signatures/
├── wes/
│   ├── somatic_variants/
│   ├── tmb_results/
│   └── cnv_results/
├── integration/
│   ├── feature_matrix/
│   ├── biomarker_features/
│   └── final_features/
└── qc/
    ├── quality_reports/
    └── summary_statistics/
```

## Biomarker Features

The pipeline generates the following biomarker features:

### Methylation Biomarkers (DMRs)
- Differentially methylated regions
- Immune enhancer regions
- Cancer driver-adjacent enhancers
- Super-enhancers

### Fragmentomics Biomarkers
- Short fragment ratios
- Nucleosome footprinting (WPS)
- End motif frequencies
- Fragment entropy

### Expression Biomarkers
- Differentially expressed genes
- Immune signature scores
- Pathway enrichment scores

### Variant Biomarkers
- Somatic mutations (TP53, KRAS, PIK3CA, APC)
- Tumor mutational burden (TMB)
- Copy number variations
- Mutational signatures

## Configuration

### Resource Allocation
Modify the configuration files to adjust resource allocation:
- `configs/local.config`: Local execution
- `configs/hpc.config`: HPC/SLURM execution
- `configs/aws.config`: AWS Batch execution

### Parameters
Key parameters can be adjusted:
- `methylation_coverage`: Target coverage for methylation analysis
- `rnaseq_reads`: Target read count for RNA-seq
- `wes_coverage`: Target coverage for WES
- `dmr_fdr_threshold`: FDR threshold for DMR calling
- `methylation_diff_threshold`: Methylation difference threshold
- `tmb_threshold`: TMB threshold for filtering

## Docker Images

The pipeline uses specialized Docker images:
- `centaur-base`: Base image with common tools
- `centaur-methylation`: Methylation analysis tools
- `centaur-rnaseq`: RNA-seq analysis tools
- `centaur-wes`: WES analysis tools
- `centaur-integration`: Multi-omics integration tools

## Troubleshooting

### Common Issues
1. **Docker permission errors**: Ensure Docker is running and user has permissions
2. **Memory issues**: Increase memory allocation in configuration files
3. **Reference file errors**: Verify reference files are properly formatted
4. **Sample sheet errors**: Check CSV format and file paths

### Support
For issues and questions:
- Check the logs in `data/results/logs/`
- Review the QC reports in `data/results/qc/`
- Consult the Nextflow documentation

## Citation

If you use CENTAUR in your research, please cite:

```
CENTAUR: A Scalable Multi-Omics Platform for Early Cancer Biomarker Discovery
[Your citation information]
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## Acknowledgments

- Nextflow community
- Bioinformatics tool developers
- Cancer research community