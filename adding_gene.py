import sqlite3
import pandas as pd
import argparse

# Function to parse the GFF file and extract gene information
def parse_gff(gff_file):
    genes = []
    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if parts[2] == 'gene':
                chrom = parts[0]
                if chrom.startswith('NC_'):
                    chrom = chrom.replace('NC_', '').split('.')[0]  # Simplified chromosome format
                start = int(parts[3])
                end = int(parts[4])
                info = parts[8]
                gene_id = None
                gene_name = None
                for entry in info.split(';'):
                    if entry.startswith('Dbxref=GeneID:'):
                        gene_id = entry.split(':')[1].split(',')[0]  # Correctly extract GeneID
                    elif entry.startswith('Name='):
                        gene_name = entry.split('=')[1]
                genes.append({'chrom': chrom, 'start': start, 'end': end, 'gene_id': gene_id, 'gene_name': gene_name})
    return pd.DataFrame(genes)

def convert_chrom(chrom):
    if chrom.isdigit():
        return int(chrom)
    elif chrom == 'X':
        return 23
    elif chrom == 'Y':
        return 24
    else:
        return None

# Function to update the database with gene information
def update_database(db_path, genes_df):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Check if the columns already exist
    cursor.execute("PRAGMA table_info(variants)")
    columns = [column[1] for column in cursor.fetchall()]
    
    if 'gene_symbol' not in columns:
        cursor.execute("ALTER TABLE variants ADD COLUMN gene_symbol TEXT")
    if 'gene_id' not in columns:
        cursor.execute("ALTER TABLE variants ADD COLUMN gene_id TEXT")
    
    conn.commit()

    # Update the database with gene information
    variants_df = pd.read_sql_query("SELECT * FROM variants", conn)
    variants_df['chrom_int'] = variants_df['chrom'].apply(convert_chrom)
    
    for index, variant in variants_df.iterrows():
        chrom = variant['chrom_int']
        pos = variant['pos']
        
        matching_genes = genes_df[(genes_df['chrom'] == chrom) & (genes_df['start'] <= pos) & (genes_df['end'] >= pos)]
        if not matching_genes.empty:
            gene_names = ','.join(matching_genes['gene_name'].unique())
            gene_ids = ','.join(filter(None, matching_genes['gene_id'].unique()))
            cursor.execute("""
                UPDATE variants 
                SET gene_symbol = ?, gene_id = ?
                WHERE chrom = ? AND pos = ?
            """, (gene_names, gene_ids, variant['chrom'], pos))
    
    conn.commit()
    conn.close()

def main():
    parser = argparse.ArgumentParser(description="Update a SQLite database with gene information from a GFF file.")
    parser.add_argument('db_path', help="Path to the SQLite database")
    parser.add_argument('gff', help="Path to the GFF File")
    args = parser.parse_args()

    db_path = args.db_path
    gff_file = args.gff

    genes_df = parse_gff(gff_file)
    genes_df['chrom'] = genes_df['chrom'].apply(convert_chrom)
    genes_df = genes_df.dropna(subset=['chrom'])  # Only keep rows with valid chromosome numbers
    print(genes_df.head())

    update_database(db_path, genes_df)
    print("Database updated with gene symbols and IDs.")

if __name__ == "__main__":
    main()
