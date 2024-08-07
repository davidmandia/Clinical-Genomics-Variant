import pandas as pd
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

def classify_indels(df, rare_threshold=0.01):
    if "gnomadg" in df.columns:
        populations = [
            'gnomadg', 'gnomadg_eas', 'gnomadg_nfe', 'gnomadg_fin', 'gnomadg_amr', 
            'gnomadg_afr', 'gnomadg_asj', 'gnomadg_oth', 'gnomadg_sas', 'gnomadg_mid', 'gnomadg_ami'
        ]
    elif "af" in df.columns:
        populations = ["af", "af_eas", "af_nfe", "af_fin", "af_amr", "af_afr", "af_asj", "af_oth", "af_sas", "af_mid", "af_ami"]
    else:
        print("No known population frequency columns found in the database.")
        return pd.DataFrame()

    available_populations = [pop for pop in populations if pop in df.columns]

    for pop in available_populations:
        df[pop] = pd.to_numeric(df[pop], errors='coerce')

    if not available_populations:
        print("No population frequency columns available for analysis.")
        return pd.DataFrame()

    rare_indels = df[df[available_populations].le(rare_threshold).any(axis=1)]
    return rare_indels

def analyze_rare_indels(rare_indels, output_dir, db_name):
    if rare_indels.empty:
        print("No rare indels found for the analysis.")
        return

    rare_indels = rare_indels.assign(genes=rare_indels['genes'].str.split(',')).explode('genes')
    rare_indels = rare_indels[rare_indels['genes'] != '']

    gene_counts = rare_indels['genes'].value_counts()

    top_genes = gene_counts.head(20)

    plt.figure(figsize=(18, 10))
    sns.barplot(x=top_genes.index, y=top_genes.values, palette="viridis")
    plt.title(f'Top 20 Genes Most Frequently Affected by Rare Indels in {db_name}', fontsize=22)
    plt.xlabel('Gene', fontsize=20)
    plt.ylabel('Count', fontsize=20)
    plt.xticks(rotation=45, ha='right', fontsize=18)
    plt.yticks(fontsize=18)
    
    table_data = top_genes.reset_index().values
    cell_text = []
    for row in table_data:
        cell_text.append([row[0], int(row[1])])

    plt.table(cellText=cell_text, colLabels=['Gene', 'Count'], cellLoc='center', loc='right', bbox=[1.05, 0, 0.3, 1], edges='horizontal', fontsize=16)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/rare_indels_{db_name}_gene_counts.png", bbox_inches='tight')
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Analyze rare indels from a SQLite database.")
    parser.add_argument('db_path', help="Path to the SQLite database")
    args = parser.parse_args()

    db_path = args.db_path
    output_dir = "analysis/figures"
    os.makedirs(output_dir, exist_ok=True)
    
    db_name = os.path.basename(db_path).replace('.db', '')

    df = read_database(db_path)
    
    rare_indels = classify_indels(df)
    analyze_rare_indels(rare_indels, output_dir, db_name)

if __name__ == "__main__":
    main()
