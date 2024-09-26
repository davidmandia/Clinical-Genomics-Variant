import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse

def read_database(db_path):
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query("SELECT ref, alt FROM variants", conn)
    conn.close()
    return df

def label_indel(row):
    if len(row['alt']) > len(row['ref']):
        # Insertion
        return f"+{row['alt'][len(row['ref']):]}"  # Label with "+"
    elif len(row['alt']) < len(row['ref']):
        # Deletion
        return f"-{row['ref'][len(row['alt']):]}"  # Label with "-"
    else:
        return None

def analyze_indels(df):
    # Apply the labeling function to each row
    df['indel'] = df.apply(label_indel, axis=1)
    
    # Remove rows where ref == alt (no indel)
    df = df.dropna(subset=['indel'])
    
    # Group by the indel sequence and count occurrences
    indel_counts = df['indel'].value_counts().reset_index()
    indel_counts.columns = ['indel', 'count']

    return indel_counts

def plot_indels(indel_counts, output_dir, db_name):
    plt.figure(figsize=(12, 8))
    sns.barplot(x='count', y='indel', data=indel_counts.head(20), palette='Blues_d')
    plt.title(f'Most Common Indels in {db_name}', fontsize=32)
    plt.xlabel('Count', fontsize=28)
    plt.ylabel('Indel Sequence', fontsize=28)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.tight_layout()

    # Save plot as a PNG file
    output_path = f"{output_dir}/common_indels_{db_name}.png"
    plt.savefig(output_path, dpi=300)
    plt.show()
    print(f"Plot saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Analyze common indels from a SQLite database.")
    parser.add_argument('db_path', help="Path to the SQLite database")
    args = parser.parse_args()

    db_path = args.db_path
    output_dir = "analysis/figures"
    os.makedirs(output_dir, exist_ok=True)
    
    db_name = os.path.basename(db_path).replace('.db', '')
    
    df = read_database(db_path)
    indel_counts = analyze_indels(df)
    
    plot_indels(indel_counts, output_dir, db_name)

if __name__ == "__main__":
    main()
