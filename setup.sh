#!/bin/bash

# CENTAUR Multi-Omics Pipeline Setup Script
# This script sets up the pipeline environment and builds Docker images

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if Docker is installed
check_docker() {
    if ! command -v docker &> /dev/null; then
        print_error "Docker is not installed. Please install Docker first."
        exit 1
    fi
    print_status "Docker is installed"
}

# Check if Nextflow is installed
check_nextflow() {
    if ! command -v nextflow &> /dev/null; then
        print_warning "Nextflow is not installed. Installing Nextflow..."
        curl -s https://get.nextflow.io | bash
        sudo mv nextflow /usr/local/bin/
        print_status "Nextflow installed"
    else
        print_status "Nextflow is installed"
    fi
}

# Build Docker images
build_docker_images() {
    print_status "Building Docker images..."
    
    # Build base image
    print_status "Building base image..."
    docker build -t centaur-base:latest -f docker/Dockerfile.base .
    
    # Build methylation image
    print_status "Building methylation analysis image..."
    docker build -t centaur-methylation:latest -f docker/Dockerfile.methylation .
    
    # Build RNA-seq image
    print_status "Building RNA-seq analysis image..."
    docker build -t centaur-rnaseq:latest -f docker/Dockerfile.rnaseq .
    
    # Build WES image
    print_status "Building WES analysis image..."
    docker build -t centaur-wes:latest -f docker/Dockerfile.wes .
    
    # Build integration image
    print_status "Building integration analysis image..."
    docker build -t centaur-integration:latest -f docker/Dockerfile.integration .
    
    print_status "All Docker images built successfully"
}

# Create sample sheet template
create_sample_sheet() {
    print_status "Creating sample sheet template..."
    cat > data/sample_sheet_template.csv << EOF
sample_id,cancer_type,stage,sample_type,cfDNA_file,rnaseq_file,wes_tumor_file,wes_normal_file
SAMPLE_001,CRC,Early,Cancer,sample_001_cfdna_R1.fastq.gz,sample_001_rnaseq_R1.fastq.gz,sample_001_tumor.bam,sample_001_normal.bam
SAMPLE_002,CRC,Early,Control,sample_002_cfdna_R1.fastq.gz,sample_002_rnaseq_R1.fastq.gz,sample_002_tumor.bam,sample_002_normal.bam
SAMPLE_003,Lung,Late,Cancer,sample_003_cfdna_R1.fastq.gz,sample_003_rnaseq_R1.fastq.gz,sample_003_tumor.bam,sample_003_normal.bam
SAMPLE_004,Breast,Early,Cancer,sample_004_cfdna_R1.fastq.gz,sample_004_rnaseq_R1.fastq.gz,sample_004_tumor.bam,sample_004_normal.bam
EOF
    print_status "Sample sheet template created at data/sample_sheet_template.csv"
}

# Create reference directory structure
create_reference_structure() {
    print_status "Creating reference directory structure..."
    mkdir -p references/{genome,annotation,intervals}
    
    print_status "Reference directory structure created"
    print_warning "Please add your reference files to the references/ directory:"
    print_warning "  - references/genome/hg38.fa (reference genome)"
    print_warning "  - references/genome/hg38 (genome index)"
    print_warning "  - references/annotation/hg38.gtf (gene annotation)"
    print_warning "  - references/intervals/exome_targets.bed (target intervals)"
    print_warning "  - references/annotation/annovar_db (ANNOVAR database)"
}

# Create test data structure
create_test_structure() {
    print_status "Creating test data structure..."
    mkdir -p data/{raw,processed,results}/{methylation,rnaseq,wes,integration}
    mkdir -p data/{raw,processed,results}/qc
    mkdir -p data/{raw,processed,results}/reports
    
    print_status "Test data structure created"
}

# Main setup function
main() {
    print_status "Setting up CENTAUR Multi-Omics Pipeline..."
    
    # Check prerequisites
    check_docker
    check_nextflow
    
    # Create directory structure
    create_reference_structure
    create_test_structure
    
    # Create sample sheet template
    create_sample_sheet
    
    # Build Docker images
    build_docker_images
    
    print_status "CENTAUR Multi-Omics Pipeline setup completed!"
    print_status "Next steps:"
    print_status "1. Add your reference files to the references/ directory"
    print_status "2. Update the sample sheet with your data"
    print_status "3. Run the pipeline with: nextflow run workflows/main.nf --input_dir data/raw --sample_sheet data/sample_sheet.csv"
}

# Run main function
main "$@"
