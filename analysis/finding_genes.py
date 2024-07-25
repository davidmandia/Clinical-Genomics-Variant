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

def main():
    parser = argparse.ArgumentParser(description="Count occurrences of clinically relevant genes labeled 'Green' in the database.")
    parser.add_argument('db_path', help="Path to the SQLite database")
    args = parser.parse_args()

    db_path = args.db_path
    gene_counts_df = get_clinically_relevant_gene_counts(db_path)

    print(gene_counts_df.head(10))  # Display top 10 genes for brevity
    gene_counts_df.to_csv("green_clinically_relevant_gene_counts.csv", index=False)  # Save the full list to a CSV file

if __name__ == "__main__":
    main()
