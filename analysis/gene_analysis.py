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

def analyze_coding_vs_noncoding(df, output_dir):
    # Classify regions as coding or non-coding
    coding_terms = ['coding_sequence_variant', 'missense_variant', 'synonymous_variant', 'frameshift_variant', 'stop_gained', 'stop_lost', 'start_lost']
    
    # Non-coding terms (added for reference)
    non_coding_terms = ['upstream_gene_variant', 'downstream_gene_variant', 'regulatory_region_variant',
    'TF_binding_site_variant', 'promoter_variant', 'non_coding_transcript_exon_variant',
    'non_coding_transcript_variant', 'splice_region_variant', 'splice_donor_variant',
    'splice_acceptor_variant', 'intron_variant', 'intergenic_variant']
    
    df['region'] = df['consequences'].apply(lambda x: 'coding' if any(term in x for term in coding_terms) else 'non-coding')
    
    region_counts = df['region'].value_counts()
    
    # Plot the distribution of indels in coding vs. non-coding regions
    plt.figure(figsize=(8, 5))
    region_counts.plot(kind='bar', color=['blue', 'orange'])
    plt.title('Distribution of Indels in Coding vs. Non-coding Regions')
    plt.xlabel('Region')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/indel_distribution_coding_vs_noncoding.png")
    plt.show()

def main():
    db_path = "output/database/GRCh38_indels_variant.db"
    output_dir = "analysis/figures"
    os.makedirs(output_dir, exist_ok=True)
    
    df = read_database(db_path)
    
    analyze_gene_frequency(df, output_dir)
    analyze_coding_vs_noncoding(df, output_dir)

if __name__ == "__main__":
    main()
