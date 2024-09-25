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

def extract_transcript_info(df):
    # Extract the transcript version
    df['transcript_version'] = df['transcript_ref'].str.extract(r'\.(\d+)$').astype(int)
    return df

def plot_transcript_versions(df, output_dir, db_name):
    # Count the occurrences of each transcript version
    version_counts = df['transcript_version'].value_counts().sort_index()
    print("version count", version_counts)

    # Plotting
    plt.figure(figsize=(12, 8))
    sns.barplot(x=version_counts.index, y=version_counts.values, palette="pastel", alpha=0.7)
    plt.title(f'Number of Transcript Versions in {db_name}', fontsize=20)
    plt.xlabel('Transcript Version', fontsize=16)
    plt.ylabel('Count of Transcripts', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()

    output_path = os.path.join(output_dir, f'transcript_versions_{db_name}.png')
    plt.savefig(output_path, dpi=300)
    print(f"Plot saved to {output_path}")

    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Analyze transcript versions from a SQLite database.")
    parser.add_argument('db_path', help="Path to the SQLite database")
    args = parser.parse_args()

    db_path = args.db_path
    output_dir = "analysis/figures"
    os.makedirs(output_dir, exist_ok=True)
    
    db_name = os.path.basename(db_path).replace('.db', '')
    
    df = read_database(db_path)
    
    df = extract_transcript_info(df)
    
    plot_transcript_versions(df, output_dir, db_name)

if __name__ == "__main__":
    main()
