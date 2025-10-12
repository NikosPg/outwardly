#!/bin/bash
# Quick setup script for RNA-seq analysis pipeline

set -e

echo "=========================================="
echo "RNA-seq Analysis Pipeline Setup"
echo "OutWardly Omics & Bioinformatics"
echo "=========================================="
echo ""

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "⚠️  Conda not found. Please install Miniconda or Anaconda first."
    echo "   Download from: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "✓ Conda found"

# Create conda environment
echo ""
echo "Creating conda environment from environment.yml..."
conda env create -f environment.yml

echo ""
echo "=========================================="
echo "Setup complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Activate the environment:"
echo "   conda activate rnaseq-analysis"
echo ""
echo "2. Edit config.yaml to specify:"
echo "   - Input data locations"
echo "   - Reference genome paths"
echo "   - Sample metadata"
echo ""
echo "3. Run the pipeline:"
echo "   snakemake --cores 8"
echo ""
echo "4. Or do a dry run first:"
echo "   snakemake -n"
echo ""
echo "5. Visualize the workflow:"
echo "   snakemake --dag | dot -Tpng > workflow.png"
echo ""
echo "For help: hello@outwardly.net"
echo "=========================================="
