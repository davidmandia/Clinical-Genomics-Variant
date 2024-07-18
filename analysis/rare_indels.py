import pandas as pd
import matplotlib.pyplot as plt
import sqlite3
import os

def read_database(db_path):
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query("SELECT * FROM variants", conn)
    conn.close()
    return df

def classify_indels(df, rare_threshold=0.01):
    populations = [
        'gnomadg', 'gnomadg_eas', 'gnomadg_nfe', 'gnomadg_fin', 'gnomadg_amr', 
        'gnomadg_afr', 'gnomadg_asj', 'gnomadg_oth', 'gnomadg_sas', 'gnomadg_mid', 'gnomadg_ami'
    ]
    
    # Convert population frequency columns to numeric, forcing errors to NaN
    for pop in populations:
        df[pop] = pd.to_numeric(df[pop], errors='coerce')
    
    # Filter rare indels based on the defined threshold
    rare_indels = df[df[populations].le(rare_threshold).any(axis=1)]
    
    return rare_indels

def analyze_rare_indels(rare_indels, output_dir):
    # Aggregate genes affected by rare indels
    rare_indels = rare_indels.assign(genes=rare_indels['genes'].str.split(',')).explode('genes')
    
    # Filter out empty gene names
    rare_indels = rare_indels[rare_indels['genes'] != '']
    
    # Count the frequency of each gene affected by rare indels
    gene_counts = rare_indels['genes'].value_counts()
    
    # Plotting the top 20 genes most frequently affected by rare indels
    top_genes = gene_counts.head(20)
    plt.figure(figsize=(12, 8))
    top_genes.plot(kind='bar', color='teal')
    plt.title('Top 20 Genes Most Frequently Affected by Rare Indels')
    plt.xlabel('Gene')
    plt.ylabel('Count')
    plt.xticks(rotation=45, ha='right')
    
    # Adding table to the plot
    plt.table(cellText=top_genes.reset_index().values,
              colLabels=['Gene', 'Count'],
              cellLoc = 'center', loc='bottom', bbox=[0, -0.3, 1, 0.25])
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/rare_indels_gene_counts.png")
    plt.show()

def main():
    db_path = "output/database/GRCh38_indels_variant.db"
    output_dir = "analysis/figures"
    os.makedirs(output_dir, exist_ok=True)
    
    df = read_database(db_path)
    
    rare_indels = classify_indels(df)
    analyze_rare_indels(rare_indels, output_dir)

if __name__ == "__main__":
    main()
