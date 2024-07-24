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



def analyze_coding_vs_noncoding(df, output_dir, db_name):
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
    plt.title(f'Distribution of Indels in Coding vs. Non-coding Regions')
    plt.xlabel('Region')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/indel_{db_name}distribution_coding_vs_noncoding.png")
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Analyze gene frequency from a SQLite database.")
    parser.add_argument('db_path', help="Path to the SQLite database")
    args = parser.parse_args()

    db_path = args.db_path
    
    output_dir = "analysis/figures"
    os.makedirs(output_dir, exist_ok=True)
    
    db_name = os.path.basename(db_path).replace('.db', '')
    
    df = read_database(db_path)
    

    analyze_coding_vs_noncoding(df, output_dir, db_name)

if __name__ == "__main__":
    main()
