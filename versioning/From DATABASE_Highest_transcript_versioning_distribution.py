import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sqlite3
import argparse
import os

def read_database(db_path):
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query("SELECT DISTINCT transcript_ref FROM variants", conn)
    conn.close()
    return df

def plot_reference_types(df, output_dir, db_name):
    # Extract reference types from transcript_ref column
    df['ref_type'] = df['transcript_ref'].str.extract(r'(\w+_)')
    
    ref_type_counts = df['ref_type'].value_counts()
    
    plt.figure(figsize=(12, 8))
    sns.barplot(x=ref_type_counts.values, y=ref_type_counts.index, palette="pastel")

    # Add count labels on each bar
    for i, v in enumerate(ref_type_counts.values):
        plt.text(v + 0.02, i, str(v), color='black', va='center', fontsize=16)

    plt.title(f'Distribution of Reference Types in {db_name}', fontsize=24)
    plt.xlabel('Count', fontsize=20)
    plt.ylabel('Reference Type', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.tight_layout()

    # Save plot as a PNG file
    output_path = os.path.join(output_dir, f'reference_types_distribution_{db_name}.png')
    plt.savefig(output_path, dpi=300)
    print(f"Plot saved to {output_path}")

    plt.show()


def get_highest_transcript_versions(df, output_dir, db_name):
    # Extract transcript base and version
    df['transcript_base'] = df['transcript_ref'].str.extract(r'(\w+_\d+)')
    df['transcript_version'] = df['transcript_ref'].str.extract(r'\.(\d+)$').astype(int)
    
    # Identify the highest version for each transcript base
    idx = df.groupby('transcript_base')['transcript_version'].idxmax()
    highest_versions_df = df.loc[idx]

    output_path = os.path.join(output_dir, f'highest_transcript_versions_{db_name}.csv')
    highest_versions_df.to_csv(output_path, index=False)
    print(f"Highest transcript versions saved to {output_path}")

    return highest_versions_df

def main():
    parser = argparse.ArgumentParser(description="Analyze reference types and transcript versions from a SQLite database.")
    parser.add_argument('db_path', help="Path to the SQLite database")
    args = parser.parse_args()

    db_path = args.db_path
    output_dir = "analysis/results"
    os.makedirs(output_dir, exist_ok=True)
    
    db_name = os.path.basename(db_path).replace('.db', '')
    
    df = read_database(db_path)
    
    plot_reference_types(df, output_dir, db_name)
    
    #highest_versions_df = get_highest_transcript_versions(df, output_dir, db_name)
    
    #print("Rows with the highest transcript versions:")
    #print(highest_versions_df)

if __name__ == "__main__":
    main()
