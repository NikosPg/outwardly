#!/usr/bin/env python3
"""
Database utilities for storing RNA-seq metadata and results
Supports both DuckDB and PostgreSQL
"""

import pandas as pd
import duckdb
import yaml
from pathlib import Path

def create_duckdb_database(db_path, counts_file, deseq_file, metadata_file):
    """
    Create DuckDB database and store RNA-seq results

    Args:
        db_path: Path to DuckDB database file
        counts_file: Count matrix CSV
        deseq_file: DESeq2 results CSV
        metadata_file: Sample metadata CSV
    """
    # Connect to DuckDB
    conn = duckdb.connect(str(db_path))

    print(f"Creating DuckDB database at: {db_path}")

    # Load data
    counts = pd.read_csv(counts_file)
    deseq_results = pd.read_csv(deseq_file)
    metadata = pd.read_csv(metadata_file)

    # Create tables
    conn.execute("CREATE TABLE IF NOT EXISTS sample_metadata AS SELECT * FROM metadata")
    conn.execute("CREATE TABLE IF NOT EXISTS expression_counts AS SELECT * FROM counts")
    conn.execute("CREATE TABLE IF NOT EXISTS deseq2_results AS SELECT * FROM deseq_results")

    # Create indexed view for significant genes
    conn.execute("""
        CREATE OR REPLACE VIEW significant_genes AS
        SELECT *
        FROM deseq2_results
        WHERE padj < 0.05 AND ABS(log2FoldChange) > 1.0
        ORDER BY padj
    """)

    # Create summary statistics table
    conn.execute("""
        CREATE TABLE IF NOT EXISTS analysis_summary AS
        SELECT
            COUNT(*) as total_genes,
            SUM(CASE WHEN padj < 0.05 AND log2FoldChange > 1.0 THEN 1 ELSE 0 END) as upregulated,
            SUM(CASE WHEN padj < 0.05 AND log2FoldChange < -1.0 THEN 1 ELSE 0 END) as downregulated,
            SUM(CASE WHEN padj >= 0.05 OR ABS(log2FoldChange) <= 1.0 THEN 1 ELSE 0 END) as not_significant
        FROM deseq2_results
    """)

    # Print summary
    summary = conn.execute("SELECT * FROM analysis_summary").fetchdf()
    print("\nAnalysis Summary:")
    print("=" * 50)
    print(summary.to_string(index=False))

    print("\nTop 10 upregulated genes:")
    top_up = conn.execute("""
        SELECT gene_id, log2FoldChange, padj
        FROM significant_genes
        WHERE log2FoldChange > 0
        ORDER BY padj
        LIMIT 10
    """).fetchdf()
    print(top_up.to_string(index=False))

    print("\nTop 10 downregulated genes:")
    top_down = conn.execute("""
        SELECT gene_id, log2FoldChange, padj
        FROM significant_genes
        WHERE log2FoldChange < 0
        ORDER BY padj
        LIMIT 10
    """).fetchdf()
    print(top_down.to_string(index=False))

    # Create indexes for faster queries
    conn.execute("CREATE INDEX IF NOT EXISTS idx_gene_id ON deseq2_results(gene_id)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_padj ON deseq2_results(padj)")

    conn.close()
    print(f"\nDatabase created successfully: {db_path}")
    print("Tables: sample_metadata, expression_counts, deseq2_results, analysis_summary")
    print("Views: significant_genes")

def create_postgresql_database(config, counts_file, deseq_file, metadata_file):
    """
    Create PostgreSQL database and store RNA-seq results

    Args:
        config: Database configuration dict
        counts_file: Count matrix CSV
        deseq_file: DESeq2 results CSV
        metadata_file: Sample metadata CSV
    """
    try:
        import psycopg2
        from sqlalchemy import create_engine
    except ImportError:
        print("PostgreSQL support requires psycopg2 and sqlalchemy")
        print("Install with: pip install psycopg2-binary sqlalchemy")
        return

    # Create connection string
    conn_string = (
        f"postgresql://{config['user']}:{config['password']}"
        f"@{config['host']}:{config['port']}/{config['database']}"
    )

    # Create SQLAlchemy engine
    engine = create_engine(conn_string)

    print(f"Connecting to PostgreSQL: {config['host']}:{config['port']}/{config['database']}")

    # Load data
    counts = pd.read_csv(counts_file)
    deseq_results = pd.read_csv(deseq_file)
    metadata = pd.read_csv(metadata_file)

    # Store in PostgreSQL
    metadata.to_sql('sample_metadata', engine, if_exists='replace', index=False)
    counts.to_sql('expression_counts', engine, if_exists='replace', index=False)
    deseq_results.to_sql('deseq2_results', engine, if_exists='replace', index=False)

    print("Data stored in PostgreSQL successfully")

if __name__ == "__main__":
    # Get inputs from Snakemake
    counts_file = snakemake.input.counts
    deseq_file = snakemake.input.deseq
    metadata_file = snakemake.input.metadata
    output_db = snakemake.output[0]

    # Load config to determine database type
    with open("config.yaml", 'r') as f:
        config = yaml.safe_load(f)

    db_config = config.get('database', {})
    db_type = db_config.get('type', 'duckdb')

    if db_type == 'duckdb':
        create_duckdb_database(output_db, counts_file, deseq_file, metadata_file)
    elif db_type == 'postgresql':
        create_postgresql_database(db_config, counts_file, deseq_file, metadata_file)
        # Create empty file to satisfy Snakemake output
        Path(output_db).touch()
    else:
        print(f"Unknown database type: {db_type}")
