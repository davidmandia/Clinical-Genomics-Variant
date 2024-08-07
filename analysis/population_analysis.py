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

def compare_population_frequencies(df, output_dir, db_name):
    # Define populations based on the columns present in the database
    if "gnomadg" in df.columns:
        populations = [
            'gnomadg', 'gnomadg_eas', 'gnomadg_nfe', 'gnomadg_fin', 'gnomadg_amr', 
            'gnomadg_afr', 'gnomadg_asj', 'gnomadg_oth', 'gnomadg_sas', 'gnomadg_mid', 'gnomadg_ami'
        ]
        population_labels = [
            'Global', 'East Asian', 'Non-Finnish European', 'Finnish', 'American', 
            'African', 'Ashkenazi Jewish', 'Other', 'South Asian', 'Middle Eastern', 'Amish'
        ]
    elif "af" in df.columns:
        populations = ["af", "af_eas", "af_nfe", "af_fin", "af_amr", "af_afr", "af_asj", "af_oth", "af_sas", "af_mid", "af_ami"]
        population_labels = [
            'Global', 'East Asian', 'Non-Finnish European', 'Finnish', 'American', 
            'African', 'Ashkenazi Jewish', 'Other', 'South Asian', 'Middle Eastern', 'Amish'
        ]
    else:
        print("No known population frequency columns found in the database.")
        return pd.DataFrame()

    available_populations = [pop for pop in populations if pop in df.columns]

    # Convert population frequency columns to numeric, forcing errors to NaN
    for pop in available_populations:
        df[pop] = pd.to_numeric(df[pop], errors='coerce')

    if not available_populations:
        print("No population frequency columns available for analysis.")
        return pd.DataFrame()
    
    # Filter out rows with no frequency data across all population columns
    df_freq = df.dropna(subset=available_populations, how='all')

    # Debugging: Check for NaN values
    print("NaN values check:")
    print(df_freq[available_populations].isna().sum())

    # Preparing data for Seaborn
    df_melted = df_freq.melt(value_vars=available_populations, var_name='Population', value_name='Allele Frequency')
    df_melted.dropna(subset=['Allele Frequency'], inplace=True)

    # Create a custom palette with transparent colors
    base_palette = sns.color_palette("pastel", len(available_populations))
    transparent_palette = [(r, g, b, 0.5) for r, g, b in base_palette]

    plt.figure(figsize=(18, 12))  # Increased figure size
    sns.boxplot(x='Allele Frequency', y='Population', data=df_melted, palette=transparent_palette)

    # Annotate with counts of non-null values in the center of the plot
    for i, pop in enumerate(available_populations):
        count = df_freq[pop].dropna().count()
        plt.text(1.05, i, f'Count: {count}', va='center', ha='left', color='blue', fontsize=14)

    # Add legend
    handles = [plt.Line2D([0], [0], color=base_palette[i], lw=4) for i in range(len(available_populations))]
    plt.legend(handles, population_labels, title="Population", bbox_to_anchor=(1.25, 1), loc='upper left')

    plt.title(f'Distribution of Indel Frequencies Across Populations in {db_name}', fontsize=24)
    plt.xlabel('Allele Frequency', fontsize=20)
    plt.ylabel('Population', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.tight_layout()

    # Save plot as a PNG file
    output_path = os.path.join(output_dir, f'indel_{db_name}_frequencies_across_populations.png')
    plt.savefig(output_path, dpi=300)
    print(f"Plot saved to {output_path}")

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
