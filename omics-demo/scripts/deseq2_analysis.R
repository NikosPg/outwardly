#!/usr/bin/env Rscript
# DESeq2 differential expression analysis

library(DESeq2)
library(tidyverse)

# Load count matrix
counts <- read.csv(snakemake@input$counts, row.names = 1)

# Load sample metadata
metadata <- read.csv(snakemake@input$metadata, row.names = 1)

# Ensure sample order matches between counts and metadata
metadata <- metadata[colnames(counts), , drop = FALSE]

cat("Samples:", nrow(metadata), "\n")
cat("Genes:", nrow(counts), "\n")

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)

# Run differential expression analysis
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("condition", "treated", "control"))

# Convert to dataframe
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj)

# Add significance flag
res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1.0

cat("Significant genes (padj < 0.05, |LFC| > 1):", sum(res_df$significant, na.rm = TRUE), "\n")

# Save results
write.csv(res_df, snakemake@output$results, row.names = FALSE)

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, snakemake@output$normalized)

# Print summary
cat("\nDESeq2 Analysis Summary:\n")
cat("========================\n")
summary(res)

# Top 10 upregulated genes
cat("\nTop 10 upregulated genes:\n")
print(head(res_df[res_df$log2FoldChange > 0, ], 10)[, c("gene_id", "log2FoldChange", "padj")])

# Top 10 downregulated genes
cat("\nTop 10 downregulated genes:\n")
print(head(res_df[res_df$log2FoldChange < 0, ], 10)[, c("gene_id", "log2FoldChange", "padj")])
