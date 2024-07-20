import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sqlite3
import os

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
    df['genes'] = df['genes'].apply(lambda x: ','.join(set(x.split(','))))
    df['consequences'] = df['consequences'].apply(lambda x: ','.join(set(x.split(','))))
    return df

def analyze_gene_frequency(df, output_dir):
    # Aggregate genes
    df = aggregate_genes_and_consequences(df)
    
    # Split genes and explode the DataFrame
    df_genes = df.assign(genes=df['genes'].str.split(',')).explode('genes')
    
    # Count the frequency of each gene. Exclude the first line with an empty gene name
    gene_counts = df_genes['genes'].value_counts().head(20).drop('', errors='ignore')

    print("gene count", gene_counts)
    
    # Plot the top 20 genes most frequently affected by indels
    plt.figure(figsize=(10, 6))
    gene_counts.plot(kind='bar', color='teal')
    plt.title('Top 20 Genes Most Frequently Affected by Indels')
    plt.xlabel('Gene')
    plt.ylabel('Count')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/top_genes_affected_by_indels.png")
    plt.show()


def main():
    db_path = "output/database/GRCh38_indels_variant.db"
    output_dir = "analysis/figures"
    os.makedirs(output_dir, exist_ok=True)
    
    df = read_database(db_path)
    
    analyze_gene_frequency(df, output_dir)


if __name__ == "__main__":
    main()
