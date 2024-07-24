import pandas as pd
import matplotlib.pyplot as plt
import sqlite3
import os
import argparse

def read_database(db_path):
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query("SELECT * FROM variants", conn)
    conn.close()
    return df

def print_database_columns(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("PRAGMA table_info(variants)")
    columns = cursor.fetchall()
    conn.close()
    print("Columns in the database:")
    for col in columns:
        print(col[1])  # col[1] contains the column name

def classify_indels(df, rare_threshold=0.01):
    # Define populations based on the columns present in the database
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
    #print("Available populations for analysis:", available_populations)

    # Convert population frequency columns to numeric, forcing errors to NaN
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
    plt.figure(figsize=(12, 8))
    top_genes.plot(kind='bar', color='teal')
    plt.title(f'Top 20 Genes Most Frequently Affected by Rare Indels in {db_name}')
    plt.xlabel('Gene')
    plt.ylabel('Count')
    plt.xticks(rotation=45, ha='right')

    plt.table(cellText=top_genes.reset_index().values,
              colLabels=['Gene', 'Count'],
              cellLoc = 'center', loc='bottom', bbox=[0, -0.3, 1, 0.25])

    plt.tight_layout()
    plt.savefig(f"{output_dir}/rare_indels_gene_counts.png")
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
