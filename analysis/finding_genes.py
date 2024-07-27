import sqlite3
import pandas as pd
import argparse

def get_clinically_relevant_gene_counts(db_path):
    """
    Reads the database and returns a count of distinct clinically relevant genes labeled 'Green'.

    Args:
        db_path (str): Path to the SQLite database.

    Returns:
        pd.DataFrame: DataFrame with genes and their occurrence counts.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    query = """
    SELECT gene_symbol FROM variants
    WHERE clinical_label = 'Green' AND gene_symbol IS NOT NULL
    """

    cursor.execute(query)
    rows = cursor.fetchall()
    conn.close()

    # Extract genes, split them and count occurrences
    gene_counts = {}
    for row in rows:
        genes = row[0].split(',')
        for gene in genes:
            if gene:
                gene_counts[gene] = gene_counts.get(gene, 0) + 1

    # Convert to DataFrame for easier handling
    gene_counts_df = pd.DataFrame(list(gene_counts.items()), columns=['Gene', 'Count'])
    gene_counts_df = gene_counts_df.sort_values(by='Count', ascending=False)

    return gene_counts_df

def fetch_data(db_path):
    """
    Fetch data from the SQLite database and return as a DataFrame.
    """
    conn = sqlite3.connect(db_path)
    query = "SELECT * FROM variants"
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df

def calculate_metrics(df):
    """
    Calculate various metrics from the DataFrame.
    """
    total_genes = len(df['genes'].unique())
    genes_with_af = df[df['af'] != 'Na']['genes'].nunique()
    clinically_relevant_genes = df[~df['clinical_label'].isin(["Green", "Amber", "Red"])]['genes'].nunique()
    green_flag_genes = df[df['clinical_label'] == 'Green']['genes'].nunique()
    
    total_transcripts = len(df['transcript_ref'].unique())
    green_transcripts = df[df['clinical_label'] == 'Green']['transcript_ref'].nunique()

    return {
        'Total Genes': total_genes,
        'Genes with AF Value': genes_with_af,
        'Clinically Relevant Genes': clinically_relevant_genes,
        'Green Flag Genes': green_flag_genes,
        'Total Transcripts': total_transcripts,
        'Green Flag Transcripts': green_transcripts
    }
def main():
    parser = argparse.ArgumentParser(description="Count occurrences of clinically relevant genes labeled 'Green' in the database. Also, Analyze gene and transcript data from a SQLite database.")
    parser.add_argument('db_path', help="Path to the SQLite database")
    args = parser.parse_args()

    db_path = args.db_path
    gene_counts_df = get_clinically_relevant_gene_counts(db_path)

    print(gene_counts_df.head(10))  # Display top 10 genes for brevity
    gene_counts_df.to_csv("green_clinically_relevant_gene_counts.csv", index=False)  # Save the full list to a CSV file
    
        # Fetch data from the database
    df = fetch_data(db_path)

    # Calculate metrics
    metrics = calculate_metrics(df)

    # Print results
    print(f"Total number of rows: {len(df)}")
    print("Metrics Summary:")
    for key, value in metrics.items():
        print(f"{key}: {value}")
    
    print("\nTotal number of transcripts:", metrics['Total Transcripts'])
    print("Number of transcripts involved in clinically relevant genes (Green flag only):", metrics['Green Flag Transcripts'])
    # Find genes with "Green" clinical relevance and non-null AF values
    # Find genes with "Green" clinical relevance and non-null AF values (excluding "Na")
    green_clinically_relevant_genes_with_af = df[(df['clinical_label'] == "Green") & (df['af'] != "Na")]['genes'].unique()
    
    # Convert the 'af' column to numeric, forcing errors to NaN (for safe conversion)
    df['af_numeric'] = pd.to_numeric(df['af'], errors='coerce')
    
    # Filter the genes to include only those with AF values below 0.01
    green_clinically_relevant_genes_with_af_below_001 = df[(df['clinical_label'] == "Green") & (df['af_numeric'] < 0.01)]['genes'].unique()
    
    print("Clinical relevant genes with rare indels", len(green_clinically_relevant_genes_with_af_below_001))
    print("Genes with 'Green' clinical relevance and AF values below 0.01:", green_clinically_relevant_genes_with_af_below_001)


if __name__ == "__main__":
    main()

 

# Metrics Summary:
# Total Genes: 7395
# Genes with AF Value: 30483
# Clinically Relevant Genes: 7395
# Green Flag Genes: 1061
# Total Transcripts: 10423
# Green Flag Transcripts: 3134

# Total number of transcripts: 10423
# Number of transcripts involved in clinically relevant genes (Green flag only): 3134
