import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sqlite3
import os
import argparse

def read_database(db_path):
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query("SELECT * FROM variants", conn)
    conn.close()
    return df

def convert_chrom_to_int(chrom):
    try:
        return int(chrom)
    except ValueError:
        return None

def aggregate_genes_and_consequences(df):
    df['genes'] = df['gene_symbol'].apply(lambda x: ','.join(set(x.split(','))) if pd.notnull(x) else '')
    df['consequences'] = df['consequences'].apply(lambda x: ','.join(set(x.split(','))) if pd.notnull(x) else '')
    return df

def analyze_gene_frequency(df, output_dir, db_name):
    # Aggregate genes
    df = aggregate_genes_and_consequences(df)
    
    # Split genes and explode the DataFrame
    df_genes = df.assign(gene_symbol=df['gene_symbol'].str.split(',')).explode('gene_symbol')
    
    # Clean up whitespace
    df_genes['gene_symbol'] = df_genes['gene_symbol'].str.strip()
    
    # Count the frequency of each gene. Exclude the first line with an empty gene name
    gene_counts = df_genes['gene_symbol'].value_counts().head(20).drop('', errors='ignore')

    #print("gene count", gene_counts)
    
    # Plot the top 20 genes most frequently affected by indels
    plt.figure(figsize=(10, 6))
    gene_counts.plot(kind='bar', color='teal')
    plt.title('Top 20 Genes Most Frequently Affected by Indels')
    plt.xlabel('Gene')
    plt.ylabel('Count')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/top_genes_affected_by_indels_{db_name}.png")
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Analyze gene frequency from a SQLite database.")
    parser.add_argument('db_path', help="Path to the SQLite database")
    args = parser.parse_args()

    db_path = args.db_path
    output_dir = "analysis/figures"
    os.makedirs(output_dir, exist_ok=True)

    # Extract database name for labeling
    db_name = os.path.basename(db_path).replace('.db', '')

    df = read_database(db_path)
    
    analyze_gene_frequency(df, output_dir, db_name)

if __name__ == "__main__":
    main()
