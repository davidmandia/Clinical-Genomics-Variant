import sqlite3
import pandas as pd
import requests
import argparse
import time
from concurrent.futures import ThreadPoolExecutor

def get_color_from_evidence(evidence_list):
    for evidence in evidence_list:
        if 'Green' in evidence:
            return 'Green'
        elif 'Amber' in evidence:
            return 'Amber'
        elif 'Red' in evidence:
            return 'Red'
    return 'No color label found'

def get_clinically_relevant_gene_color(gene, api_url_base, last_gene, last_color, max_retries=5, backoff_factor=1):
    if gene == last_gene:
        return last_color
    
    api_url = f"{api_url_base}{gene}"
    print(f"Testing API for gene: {gene}")
    print("api_url:", api_url)
    
    for attempt in range(max_retries):
        try:
            response = requests.get(api_url)
            if response.status_code == 200:
                data = response.json()
                if data["count"] > 0:
                    gene_data = data["results"][0]
                    evidence = gene_data.get("evidence", [])
                    color_label = get_color_from_evidence(evidence)
                    print(f"API is working correctly for gene: {gene}. Here's the color label found in the evidence:")
                    print("Color label:", color_label)
                    return gene, color_label
                else:
                    print(f"No results found for gene: {gene}")
                    return gene, 'No Data'
            elif response.status_code == 404:
                print(f"API request returned 404 for gene: {gene}")
                return gene, 'No Data'
            else:
                print(f"API request failed for gene: {gene} with status code: {response.status_code}")
                print(response.text)
        except requests.exceptions.RequestException as e:
            print(f"API request encountered an error: {e}")
            if attempt < max_retries - 1:
                sleep_time = backoff_factor * (2 ** attempt)
                print(f"Retrying in {sleep_time} seconds...")
                time.sleep(sleep_time)
            else:
                print(f"Failed to fetch data for gene: {gene} after {max_retries} attempts.")
                return gene, 'No Data'

def check_clinical_relevance(db_path, api_url_base):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    try:
        cursor.execute("ALTER TABLE variants ADD COLUMN clinically_relevant TEXT")
    except sqlite3.OperationalError:
        print("Column 'clinically_relevant' already exists.")
    
    try:
        cursor.execute("ALTER TABLE variants ADD COLUMN clinical_label TEXT")
    except sqlite3.OperationalError:
        print("Column 'clinical_label' already exists.")
    
    conn.commit()

    df = pd.read_sql_query("SELECT * FROM variants", conn)

    df['gene_symbol'] = df['gene_symbol'].fillna('')

    last_gene = None
    last_color = None
    
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

    with ThreadPoolExecutor(max_workers=10) as executor:
        results = list(executor.map(get_labels, df['gene_symbol']))

    df['clinical_label'], df['clinically_relevant'] = zip(*results)

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update a SQLite database with clinically relevant genes from PanelApp.")
    parser.add_argument('db_path', help="Path to the SQLite database")
    args = parser.parse_args()
    
    db_path = args.db_path
    panel_app_api_url = 'https://panelapp.genomicsengland.co.uk/api/v1/genes/'
    
    check_clinical_relevance(db_path, panel_app_api_url)
