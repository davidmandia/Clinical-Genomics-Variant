import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
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
    sns.set(style="whitegrid")

    fig, axs = plt.subplots(3, 1, figsize=(15, 18))

    # Plot frequency of Insertions vs. Deletions
    sns.barplot(x=indel_counts.index, y=indel_counts.values, ax=axs[0], palette='Blues', alpha=0.7)
    axs[0].set_title('Frequency of Insertions vs. Deletions', fontsize=32, pad=20)
    axs[0].set_xlabel('Indel Type', fontsize=28)
    axs[0].set_ylabel('Count', fontsize=28)
    axs[0].tick_params(axis='both', labelsize=24)
    for p in axs[0].patches:
        axs[0].annotate(format(p.get_height(), '.0f'), 
                        (p.get_x() + p.get_width() / 2., p.get_height()), 
                        ha = 'center', va = 'center', 
                        xytext = (0, 9), 
                        textcoords = 'offset points', fontsize=24)

    # Plot size distribution of Indels
    sns.barplot(x=indel_size_distribution.index, y=indel_size_distribution.values, ax=axs[1], palette='Greens', alpha=0.7)
    axs[1].set_title('Size Distribution of Indels', fontsize=32, pad=20)
    axs[1].set_xlabel('Indel Size (bp)', fontsize=28)
    axs[1].set_ylabel('Count', fontsize=28)
    axs[1].tick_params(axis='both', labelsize=24)
    for p in axs[1].patches:
        axs[1].annotate(format(p.get_height(), '.0f'), 
                        (p.get_x() + p.get_width() / 2., p.get_height()), 
                        ha = 'center', va = 'center', 
                        xytext = (0, 9), 
                        textcoords = 'offset points', fontsize=24)

    # Plot indel frequencies across chromosomes
    sns.barplot(x=chromosome_counts.index.astype(int).astype(str), y=chromosome_counts.values, ax=axs[2], palette='Purples', alpha=0.7)
    axs[2].set_title('Indel Frequency Across Chromosomes', fontsize=32, pad=20)
    axs[2].set_xlabel('Chromosome', fontsize=28)
    axs[2].set_ylabel('Count', fontsize=28)
    axs[2].tick_params(axis='both', labelsize=24)
    for p in axs[2].patches:
        axs[2].annotate(format(p.get_height(), '.0f'), 
                        (p.get_x() + p.get_width() / 2., p.get_height()), 
                        ha = 'center', va = 'center', 
                        xytext = (0, 9), 
                        textcoords = 'offset points', fontsize=24)

    plt.tight_layout(h_pad=3.0)
    plt.savefig(f"{output_dir}/indel_{db_name}_distribution_analysis.png")
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Analyze indel distribution from a SQLite database.")
    parser.add_argument('db_path', help="Path to the SQLite database")
    args = parser.parse_args()

    db_path = args.db_path
    output_dir = "analysis/figures"
    os.makedirs(output_dir, exist_ok=True)
    
    db_name = os.path.basename(db_path).replace('.db', '')

    df = read_database(db_path)
    indel_distribution_analysis(df, output_dir, db_name)

if __name__ == "__main__":
    main()
