#!/usr/bin/env nextflow
/*
 * RNA-seq Differential Expression Analysis Pipeline
 * ΕΚΦΑΝΣΙΣ Bioinformatics
 *
 * This pipeline performs:
 * - Quality control (FastQC, MultiQC)
 * - Read alignment (STAR)
 * - Feature counting (featureCounts)
 * - Differential expression (DESeq2)
 * - Visualization and database storage
 */

nextflow.enable.dsl = 2

// Print pipeline info
log.info """
    ==========================================
    RNA-seq Analysis Pipeline
    ΕΚΦΑΝΣΙΣ Omics & Bioinformatics
    ==========================================
    reads      : ${params.reads}
    outdir     : ${params.outdir}
    genome     : ${params.genome}
    gtf        : ${params.gtf}
    ==========================================
    """.stripIndent()

/*
 * Process: FastQC - Quality control
 */
process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.{html,zip}"

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}

/*
 * Process: MultiQC - Aggregate QC reports
 */
process MULTIQC {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path('*')

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

/*
 * Process: STAR alignment
 */
process STAR_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path genome_index

    output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam")
    path "${sample_id}.Log.final.out"

    script:
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${genome_index} \\
        --readFilesIn ${reads} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample_id}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --quantMode GeneCounts
    """
}

/*
 * Process: Feature counting
 */
process FEATURECOUNTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path bam_files
    path gtf

    output:
    path "counts_matrix.txt"

    script:
    def bam_list = bam_files.collect().join(' ')
    """
    featureCounts \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -o counts_matrix.txt \\
        ${bam_list}
    """
}

/*
 * Process: Preprocess counts
 */
process PREPROCESS_COUNTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path counts

    output:
    path "counts_matrix.csv"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    df = pd.read_csv('${counts}', sep='\\t', comment='#')
    count_cols = [col for col in df.columns if col.endswith('.bam')]
    counts = df[['Geneid'] + count_cols].copy()
    counts.columns = ['gene_id'] + [col.replace('.Aligned.sortedByCoord.out.bam', '')
                                     for col in count_cols]
    counts.set_index('gene_id', inplace=True)

    # Filter low counts
    counts_filtered = counts[(counts >= ${params.min_reads}).sum(axis=1) >= ${params.min_samples}]
    counts_filtered.to_csv('counts_matrix.csv')
    """
}

/*
 * Process: DESeq2 analysis
 */
process DESEQ2 {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path counts
    path metadata

    output:
    path "deseq2_results.csv"
    path "normalized_counts.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(DESeq2)
    library(tidyverse)

    counts <- read.csv('${counts}', row.names = 1)
    metadata <- read.csv('${metadata}', row.names = 1)
    metadata <- metadata[colnames(counts), , drop = FALSE]

    dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = metadata,
        design = ~ condition
    )

    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", "treated", "control"))

    res_df <- as.data.frame(res) %>%
        rownames_to_column("gene_id") %>%
        arrange(padj)

    write.csv(res_df, 'deseq2_results.csv', row.names = FALSE)

    normalized_counts <- counts(dds, normalized = TRUE)
    write.csv(normalized_counts, 'normalized_counts.csv')
    """
}

/*
 * Process: Visualization
 */
process VISUALIZE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path deseq_results
    path normalized_counts
    path metadata

    output:
    path "pca_plot.png"
    path "volcano_plot.png"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    # Load data
    results = pd.read_csv('${deseq_results}')
    counts = pd.read_csv('${normalized_counts}', index_col=0)
    metadata = pd.read_csv('${metadata}', index_col=0)

    # PCA plot
    scaler = StandardScaler()
    counts_scaled = scaler.fit_transform(counts.T)
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(counts_scaled)

    fig, ax = plt.subplots(figsize=(10, 8))
    for condition in metadata['condition'].unique():
        mask = metadata['condition'] == condition
        ax.scatter(pca_result[mask, 0], pca_result[mask, 1],
                   label=condition, s=100, alpha=0.7)
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
    ax.set_title('PCA of RNA-seq samples')
    ax.legend()
    plt.savefig('pca_plot.png', dpi=300)

    # Volcano plot
    results = results.dropna(subset=['padj', 'log2FoldChange'])
    results['neg_log10_padj'] = -np.log10(results['padj'])

    fig, ax = plt.subplots(figsize=(12, 10))
    colors = ['lightgray', 'red', 'blue']
    results['color'] = 0
    results.loc[(results['padj'] < 0.05) & (results['log2FoldChange'] > 1), 'color'] = 1
    results.loc[(results['padj'] < 0.05) & (results['log2FoldChange'] < -1), 'color'] = 2

    for i, color in enumerate(colors):
        subset = results[results['color'] == i]
        ax.scatter(subset['log2FoldChange'], subset['neg_log10_padj'],
                   c=color, alpha=0.6, s=20)

    ax.axhline(y=-np.log10(0.05), color='gray', linestyle='--')
    ax.axvline(x=1, color='gray', linestyle='--')
    ax.axvline(x=-1, color='gray', linestyle='--')
    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('-log10(adjusted p-value)')
    ax.set_title('Volcano Plot: Treated vs Control')
    plt.savefig('volcano_plot.png', dpi=300)
    """
}

/*
 * Process: Store in database
 */
process STORE_DB {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path counts
    path deseq_results
    path metadata

    output:
    path "metadata.db"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import duckdb

    conn = duckdb.connect('metadata.db')

    counts = pd.read_csv('${counts}')
    results = pd.read_csv('${deseq_results}')
    metadata = pd.read_csv('${metadata}')

    conn.execute("CREATE TABLE expression_counts AS SELECT * FROM counts")
    conn.execute("CREATE TABLE deseq2_results AS SELECT * FROM results")
    conn.execute("CREATE TABLE sample_metadata AS SELECT * FROM metadata")

    conn.execute('''
        CREATE VIEW significant_genes AS
        SELECT * FROM deseq2_results
        WHERE padj < 0.05 AND ABS(log2FoldChange) > 1.0
    ''')

    conn.close()
    """
}

/*
 * Workflow
 */
workflow {
    // Input channel
    reads_ch = Channel
        .fromFilePairs(params.reads, checkIfExists: true)

    // Quality control
    fastqc_ch = FASTQC(reads_ch)
    MULTIQC(fastqc_ch.collect())

    // Alignment
    genome_index_ch = Channel.fromPath(params.star_index)
    aligned_ch = STAR_ALIGN(reads_ch, genome_index_ch)

    // Feature counting
    bam_ch = aligned_ch[0].map { sample_id, bam -> bam }.collect()
    gtf_ch = Channel.fromPath(params.gtf)
    counts_raw = FEATURECOUNTS(bam_ch, gtf_ch)

    // Preprocessing
    counts_clean = PREPROCESS_COUNTS(counts_raw)

    // Differential expression
    metadata_ch = Channel.fromPath(params.metadata)
    deseq_out = DESEQ2(counts_clean, metadata_ch)

    // Visualization
    VISUALIZE(deseq_out[0], deseq_out[1], metadata_ch)

    // Database
    STORE_DB(counts_clean, deseq_out[0], metadata_ch)
}

workflow.onComplete {
    log.info """
        ==========================================
        Pipeline completed!
        Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Duration  : ${workflow.duration}
        Output    : ${params.outdir}
        ==========================================
        """.stripIndent()
}
