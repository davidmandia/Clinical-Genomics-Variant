import pandas as pd
import matplotlib.pyplot as plt
import sqlite3
import os

def read_database(db_path):
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query("SELECT * FROM variants", conn)
    conn.close()
    return df

def compare_population_frequencies(df, output_dir):
    populations = [
        'gnomadg', 'gnomadg_eas', 'gnomadg_nfe', 'gnomadg_fin', 'gnomadg_amr', 
        'gnomadg_afr', 'gnomadg_asj', 'gnomadg_oth', 'gnomadg_sas', 'gnomadg_mid', 'gnomadg_ami'
    ]
    
    # Convert population frequency columns to numeric, forcing errors to NaN
    for pop in populations:
        df[pop] = pd.to_numeric(df[pop], errors='coerce')
    
    # Filter out rows with no frequency data
    df_freq = df.dropna(subset=populations, how='all')
    
    # Plotting the distribution of allele frequencies across different populations
    plt.figure(figsize=(12, 8))
    df_freq[populations].plot(kind='box', vert=False)
    plt.title('Distribution of Indel Frequencies Across Populations')
    plt.xlabel('Allele Frequency')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/indel_frequencies_across_populations.png")
    plt.show()

def main():
    db_path = "output/database/GRCh38_indels_variant.db"
    output_dir = "analysis/figures"
    os.makedirs(output_dir, exist_ok=True)
    
    df = read_database(db_path)
    
    compare_population_frequencies(df, output_dir)

if __name__ == "__main__":
    main()
