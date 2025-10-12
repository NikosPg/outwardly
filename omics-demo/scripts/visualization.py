#!/usr/bin/env python3
"""
Visualization scripts for RNA-seq analysis results
Generates PCA plots, volcano plots, heatmaps, etc.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

def generate_pca_plot(counts_file, metadata_file, output_file):
    """
    Generate PCA plot from normalized counts

    Args:
        counts_file: Normalized count matrix CSV
        metadata_file: Sample metadata CSV
        output_file: Output PNG file
    """
    # Load data
    counts = pd.read_csv(counts_file, index_col=0)
    metadata = pd.read_csv(metadata_file, index_col=0)

    # Transpose (samples as rows)
    counts_t = counts.T

    # Select top 500 most variable genes
    gene_vars = counts.var(axis=1)
    top_genes = gene_vars.nlargest(500).index
    counts_subset = counts.loc[top_genes, :]

    # Standardize
    scaler = StandardScaler()
    counts_scaled = scaler.fit_transform(counts_subset.T)

    # PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(counts_scaled)

    # Create dataframe for plotting
    pca_df = pd.DataFrame(
        pca_result,
        columns=['PC1', 'PC2'],
        index=counts.columns
    )
    pca_df['condition'] = metadata.loc[pca_df.index, 'condition']

    # Plot
    fig, ax = plt.subplots(figsize=(10, 8))

    for condition in pca_df['condition'].unique():
        subset = pca_df[pca_df['condition'] == condition]
        ax.scatter(
            subset['PC1'],
            subset['PC2'],
            label=condition,
            s=100,
            alpha=0.7
        )

    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    ax.set_title('PCA of RNA-seq samples (top 500 variable genes)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"PCA plot saved to: {output_file}")

def generate_volcano_plot(deseq_results_file, output_file):
    """
    Generate volcano plot from DESeq2 results

    Args:
        deseq_results_file: DESeq2 results CSV
        output_file: Output PNG file
    """
    # Load results
    results = pd.read_csv(deseq_results_file)

    # Remove NA values
    results = results.dropna(subset=['padj', 'log2FoldChange'])

    # Calculate -log10(padj)
    results['neg_log10_padj'] = -np.log10(results['padj'])

    # Define significance thresholds
    pval_threshold = 0.05
    lfc_threshold = 1.0

    # Classify genes
    results['category'] = 'Not significant'
    results.loc[
        (results['padj'] < pval_threshold) & (results['log2FoldChange'] > lfc_threshold),
        'category'
    ] = 'Upregulated'
    results.loc[
        (results['padj'] < pval_threshold) & (results['log2FoldChange'] < -lfc_threshold),
        'category'
    ] = 'Downregulated'

    # Plot
    fig, ax = plt.subplots(figsize=(12, 10))

    colors = {
        'Not significant': 'lightgray',
        'Upregulated': 'red',
        'Downregulated': 'blue'
    }

    for category in ['Not significant', 'Upregulated', 'Downregulated']:
        subset = results[results['category'] == category]
        ax.scatter(
            subset['log2FoldChange'],
            subset['neg_log10_padj'],
            c=colors[category],
            label=f"{category} ({len(subset)})",
            alpha=0.6,
            s=20
        )

    # Add threshold lines
    ax.axhline(
        y=-np.log10(pval_threshold),
        color='gray',
        linestyle='--',
        linewidth=1,
        alpha=0.7,
        label=f'p-adj = {pval_threshold}'
    )
    ax.axvline(
        x=lfc_threshold,
        color='gray',
        linestyle='--',
        linewidth=1,
        alpha=0.7
    )
    ax.axvline(
        x=-lfc_threshold,
        color='gray',
        linestyle='--',
        linewidth=1,
        alpha=0.7,
        label=f'|LFC| = {lfc_threshold}'
    )

    ax.set_xlabel('log2 Fold Change', fontsize=12)
    ax.set_ylabel('-log10(adjusted p-value)', fontsize=12)
    ax.set_title('Volcano Plot: Treated vs Control', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', frameon=True, fancybox=True)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Volcano plot saved to: {output_file}")

    # Print summary
    print("\nVolcano plot summary:")
    print(f"Upregulated genes: {len(results[results['category'] == 'Upregulated'])}")
    print(f"Downregulated genes: {len(results[results['category'] == 'Downregulated'])}")
    print(f"Not significant: {len(results[results['category'] == 'Not significant'])}")

if __name__ == "__main__":
    # Determine which plot to generate based on output file
    output = snakemake.output[0]

    if 'pca' in output:
        generate_pca_plot(
            snakemake.input.counts,
            snakemake.input.metadata,
            output
        )
    elif 'volcano' in output:
        generate_volcano_plot(
            snakemake.input[0],
            output
        )
