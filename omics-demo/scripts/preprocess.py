#!/usr/bin/env python3
"""
Preprocessing script for RNA-seq count data
Converts featureCounts output to clean CSV format
"""

import pandas as pd
import sys

def preprocess_counts(input_file, output_file):
    """
    Read featureCounts output and convert to clean CSV

    Args:
        input_file: Path to featureCounts output
        output_file: Path to output CSV file
    """
    # Read featureCounts output (skip comment lines)
    df = pd.read_csv(input_file, sep='\t', comment='#')

    # Select gene ID and count columns (skip gene metadata columns)
    count_cols = [col for col in df.columns if col.endswith('.bam')]
    gene_cols = ['Geneid']

    counts = df[gene_cols + count_cols].copy()

    # Clean up column names (remove .bam extension)
    counts.columns = ['gene_id'] + [col.replace('.Aligned.sortedByCoord.out.bam', '')
                                     for col in count_cols]

    # Set gene_id as index
    counts.set_index('gene_id', inplace=True)

    # Filter low-count genes (at least 10 reads in at least 2 samples)
    min_reads = 10
    min_samples = 2
    counts_filtered = counts[(counts >= min_reads).sum(axis=1) >= min_samples]

    print(f"Input genes: {len(counts)}")
    print(f"Filtered genes: {len(counts_filtered)}")
    print(f"Removed: {len(counts) - len(counts_filtered)}")

    # Save to CSV
    counts_filtered.to_csv(output_file)
    print(f"Saved to: {output_file}")

if __name__ == "__main__":
    # Snakemake automatically provides these variables
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]

    preprocess_counts(input_file, output_file)
