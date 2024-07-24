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

def compare_population_frequencies(df, output_dir, db_name):
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
    
    
    # Filter out rows with no frequency data
    df_freq = df.dropna(subset=populations, how='all')
    
    # Plotting the distribution of allele frequencies across different populations
    plt.figure(figsize=(12, 8))
    df_freq[populations].plot(kind='box', vert=False)
    plt.title('Distribution of Indel Frequencies Across Populations in')
    plt.xlabel('Allele Frequency')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/indel_{db_name}frequencies_across_populations.png")
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
    
    compare_population_frequencies(df, output_dir, db_name)

if __name__ == "__main__":
    main()
