import sqlite3
import pandas as pd
import requests

def get_clinically_relevant_genes(panel_app_api_url):
    response = requests.get(panel_app_api_url)
    if response.status_code == 200:
        data = response.json()
        relevant_genes = {}
        for entry in data['results']:
            gene_symbol = entry['gene_data'].get('gene_symbol')
            evidence = entry.get('evidence', [])
            if gene_symbol:
                # Determine the highest confidence level (color label)
                confidence_levels = [ev.split()[-1].lower() for ev in evidence if 'Expert Review' in ev]
                if 'green' in confidence_levels:
                    color = 'Green'
                elif 'amber' in confidence_levels:
                    color = 'Amber'
                elif 'red' in confidence_levels:
                    color = 'Red'
                else:
                    color = 'No Data'
                relevant_genes[gene_symbol] = color
        return relevant_genes
    else:
        print(f"Failed to fetch data from PanelApp. Status code: {response.status_code}")
        return {}

def check_clinical_relevance(db_path, panel_app_api_url):
    # Fetch clinically relevant genes from PanelApp
    relevant_genes = get_clinically_relevant_genes(panel_app_api_url)
    
    if not relevant_genes:
        print("No relevant genes found. Exiting.")
        return

    # Connect to the database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    try:
        # Add clinically_relevant and clinical_label columns if they don't exist
        cursor.execute("ALTER TABLE variants ADD COLUMN clinically_relevant TEXT")
    except sqlite3.OperationalError:
        print("Column 'clinically_relevant' already exists.")
    
    try:
        cursor.execute("ALTER TABLE variants ADD COLUMN clinical_label TEXT")
    except sqlite3.OperationalError:
        print("Column 'clinical_label' already exists.")
    
    conn.commit()

    # Read the database into a DataFrame
    df = pd.read_sql_query("SELECT * FROM variants", conn)

    # Handle None values in gene_symbol column
    df['gene_symbol'] = df['gene_symbol'].fillna('')

    # Check for clinical relevance and assign color labels
    df['clinically_relevant'] = df['gene_symbol'].apply(lambda x: 'Yes' if any(gene in relevant_genes for gene in x.split(',')) else 'No')
    df['clinical_label'] = df['gene_symbol'].apply(lambda x: ','.join(set(relevant_genes[gene] for gene in x.split(',') if gene in relevant_genes)) if any(gene in relevant_genes for gene in x.split(',')) else 'No Data')

    # Update the database
    for index, row in df.iterrows():
        cursor.execute("""
            UPDATE variants 
            SET clinically_relevant = ?, clinical_label = ?
            WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ?
        """, (row['clinically_relevant'], row['clinical_label'], row['chrom'], row['pos'], row['ref'], row['alt']))
    
    conn.commit()
    conn.close()

    print("Database updated with clinically relevant genes and color labels.")

if __name__ == "__main__":
    db_path = 'output/database/GRCh38_indels_variant.db'
    panel_app_api_url = 'https://panelapp.genomicsengland.co.uk/api/v1/genes/'
    
    check_clinical_relevance(db_path, panel_app_api_url)
