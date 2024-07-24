import sqlite3
import pandas as pd
import requests
import argparse
import time
from concurrent.futures import ThreadPoolExecutor

# Function to extract color label from evidence list
def get_color_from_evidence(evidence_list):
    """
    Determines the color label based on the evidence list.

    Args:
        evidence_list (list): List of evidence strings.

    Returns:
        str: The color label ('Green', 'Amber', 'Red') or 'No color label found'.
    """
    for evidence in evidence_list:
        if 'Green' in evidence:
            return 'Green'
        elif 'Amber' in evidence:
            return 'Amber'
        elif 'Red' in evidence:
            return 'Red'
    return 'No color label found'

# Function to get the clinically relevant gene color
def get_clinically_relevant_gene_color(gene, api_url_base, last_gene=None, last_color=None, max_retries=5, backoff_factor=1):
    """
    Fetches the color label for a given gene from an API, with retry logic.

    Args:
        gene (str): The gene symbol.
        api_url_base (str): Base URL for the API.
        last_gene (str): The last gene queried.
        last_color (str): The last color obtained from the API.
        max_retries (int): Maximum number of retries for the API call.
        backoff_factor (int): Factor for exponential backoff in case of API call failure.

    Returns:
        tuple: Gene and its associated color label.
    """
    if gene == last_gene:
        return last_gene, last_color

    api_url = f"{api_url_base}{gene}"
    print(f"Testing API for gene: {gene}")

    for attempt in range(max_retries):
        try:
            response = requests.get(api_url)
            if response.status_code == 200:
                data = response.json()
                if data["count"] > 0:
                    gene_data = data["results"][0]
                    evidence = gene_data.get("evidence", [])
                    color_label = get_color_from_evidence(evidence)
                    return gene, color_label
                else:
                    return gene, 'No Data'
            elif response.status_code == 404:
                return gene, 'No Data'
            else:
                print(f"API request failed for gene: {gene} with status code: {response.status_code}")
                print(response.text)
        except requests.exceptions.RequestException as e:
            if attempt < max_retries - 1:
                sleep_time = backoff_factor * (2 ** attempt)
                time.sleep(sleep_time)
            else:
                print(f"Failed to fetch data for gene: {gene} after {max_retries} attempts.")
                return gene, 'No Data'

# Function to check and update clinical relevance in the database
def check_clinical_relevance(db_path, api_url_base):
    """
    Updates a SQLite database with clinically relevant gene information from an API.

    Args:
        db_path (str): Path to the SQLite database.
        api_url_base (str): Base URL for the API.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Add columns for clinical relevance if they do not exist
    cursor.execute("PRAGMA table_info(variants)")
    columns = [column[1] for column in cursor.fetchall()]
    
    if 'clinically_relevant' not in columns:
        cursor.execute("ALTER TABLE variants ADD COLUMN clinically_relevant TEXT")
    if 'clinical_label' not in columns:
        cursor.execute("ALTER TABLE variants ADD COLUMN clinical_label TEXT")
    
    conn.commit()

    df = pd.read_sql_query("SELECT * FROM variants", conn)
    df['gene_symbol'] = df['gene_symbol'].fillna('')

    last_gene = None
    last_color = None
    
    # Function to get labels for each gene symbol
    def get_labels(gene_symbol):
        nonlocal last_gene, last_color
        if gene_symbol == '':
            return 'No Data', 'No'
        gene_list = gene_symbol.split(',')
        colors = []
        for gene in gene_list:
            if gene == last_gene:
                color = last_color
            else:
                last_gene, color = get_clinically_relevant_gene_color(gene, api_url_base, last_gene, last_color)
                last_color = color
            colors.append(color)
        return ','.join(set(colors)), 'Yes' if 'Green' in colors or 'Amber' in colors or 'Red' in colors else 'No'

    # Use ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor(max_workers=10) as executor:
        results = list(executor.map(get_labels, df['gene_symbol']))

    df['clinical_label'], df['clinically_relevant'] = zip(*results)

    # Update the database with the new information
    for index, row in df.iterrows():
        if row['gene_symbol'] != '':
            cursor.execute("""
                UPDATE variants 
                SET clinically_relevant = ?, clinical_label = ?
                WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ?
            """, (row['clinically_relevant'], row['clinical_label'], row['chrom'], row['pos'], row['ref'], row['alt']))
    
    conn.commit()
    conn.close()
    print("Database updated with clinically relevant genes and color labels.")

# Main function to parse arguments and execute the script
def main():
    parser = argparse.ArgumentParser(description="Update a SQLite database with clinically relevant genes from PanelApp.")
    parser.add_argument('db_path', help="Path to the SQLite database")
    args = parser.parse_args()
    
    db_path = args.db_path
    panel_app_api_url = 'https://panelapp.genomicsengland.co.uk/api/v1/genes/'
    
    check_clinical_relevance(db_path, panel_app_api_url)

if __name__ == "__main__":
    main()
