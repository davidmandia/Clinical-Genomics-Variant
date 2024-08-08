import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sqlite3
import os
import argparse

def read_database(db_path):
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query("SELECT * FROM variants", conn)
    conn.close()
    return df

def plot_allele_frequency_histogram(df, output_dir, db_name):
    # Define the general population allele frequency column
    population_col = 'af'
    
    if population_col not in df.columns:
        print(f"No {population_col} column found in the database.")
        return pd.DataFrame()

    # Convert allele frequency column to numeric, forcing errors to NaN
    df[population_col] = pd.to_numeric(df[population_col], errors='coerce')

    # Filter out rows with no frequency data
    df_freq = df.dropna(subset=[population_col], how='all')

    plt.figure(figsize=(18, 12))  # Increased figure size

    # Plot histogram
    sns.histplot(df_freq, x=population_col, kde=True, color='blue', bins=50)

    plt.title(f'Distribution of Allele Frequencies in the General Population in {db_name}', fontsize=24)
    plt.xlabel('Allele Frequency', fontsize=20)
    plt.ylabel('Count', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.tight_layout()

    # Save plot as a PNG file
    output_path = os.path.join(output_dir, f'allele_frequency_distribution_{db_name}.png')
    plt.savefig(output_path, dpi=300)
    print(f"Plot saved to {output_path}")

    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Analyze allele frequency distribution from a SQLite database.")
    parser.add_argument('db_path', help="Path to the SQLite database")
    args = parser.parse_args()

    db_path = args.db_path
    output_dir = "analysis/figures"
    os.makedirs(output_dir, exist_ok=True)
    
    db_name = os.path.basename(db_path).replace('.db', '')
    
    df = read_database(db_path)
    
    plot_allele_frequency_histogram(df, output_dir, db_name)

if __name__ == "__main__":
    main()
