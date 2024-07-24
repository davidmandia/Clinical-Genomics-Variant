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

def indel_distribution_analysis(df, output_dir, db_name):
    # Classify indels as Insertions or Deletions
    df['indel_type'] = np.where(df['ref'].str.len() > df['alt'].str.len(), 'Deletion', 'Insertion')
    indel_counts = df['indel_type'].value_counts()

    # Analyze the size distribution of indels
    indel_sizes = df.apply(lambda row: abs(len(row['ref']) - len(row['alt'])), axis=1)
    indel_size_distribution = indel_sizes.value_counts().sort_index()

    # Convert chromosome to integer where possible and discard others
    df['chrom_int'] = df['chrom'].apply(convert_chrom_to_int)
    df = df.dropna(subset=['chrom_int'])

    # Compare indel frequencies across different chromosomes
    chromosome_counts = df['chrom_int'].value_counts().sort_index()

    # Plotting
    plt.figure(figsize=(12, 8))

    # Plot frequency of Insertions vs. Deletions
    plt.subplot(3, 1, 1)
    indel_counts.plot(kind='bar', color=['blue', 'orange'])
    plt.title('Frequency of Insertions vs. Deletions')
    plt.xlabel('Indel Type')
    plt.ylabel('Count')
    plt.savefig(f"{output_dir}/indel_{db_name}frequency.png")

    # Plot size distribution of Indels
    plt.subplot(3, 1, 2)
    indel_size_distribution.plot(kind='bar', color='green')
    plt.title('Size Distribution of Indels')
    plt.xlabel('Indel Size (bp)')
    plt.ylabel('Count')
    plt.savefig(f"{output_dir}/indel_{db_name}size_distribution.png")

    # Plot indel frequencies across chromosomes
    plt.subplot(3, 1, 3)
    chromosome_counts.plot(kind='bar', color='purple')
    plt.title('Indel Frequency Across Chromosomes')
    plt.xlabel('Chromosome')
    plt.ylabel('Count')
    plt.savefig(f"{output_dir}/indel_{db_name}chromosome_frequency.png")

    plt.tight_layout(h_pad=3.0)
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
    indel_distribution_analysis(df, output_dir, db_name)

if __name__ == "__main__":
    main()
