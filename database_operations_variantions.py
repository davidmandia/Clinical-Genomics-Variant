import sqlite3
import vobject
import argparse

# Define the function to create the SQLite database and table
def create_database(db_name):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    # Creating a table for variants with specified columns
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS variants (
            chrom TEXT,
            pos INTEGER,
            id TEXT,
            ref TEXT,
            alt TEXT,
            qual REAL,
            filter TEXT,
            info TEXT
        )
    ''')
    conn.commit()
    conn.close()

# Define the function to parse VCF and insert data into the database
def insert_data_from_vcf(vcf_file, db_name):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    
    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  # Skip header lines
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue  # Skip malformed lines
            
            chrom, pos, var_id, ref, alt, qual, filter_, info = parts[:8]
            qual = float(qual) if qual != '.' else None
            
            cursor.execute('''
                INSERT INTO variants (chrom, pos, id, ref, alt, qual, filter, info)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ''', (chrom, int(pos), var_id, ref, alt, qual, filter_, info))
    
    conn.commit()
    conn.close()

# Define the main function to call the above functions
def main():
    parser = argparse.ArgumentParser(description="Process a VCF file to store variants in a SQLlite database.")
    parser.add_argument('vcf',  help="Provides the VCF file to add to the database")

    args = parser.parse_args()
    db_name = 'variants.db'
    vcf_file = args.vcf
    
    create_database(db_name)
    insert_data_from_vcf(vcf_file, db_name)
    print(f'Data from {vcf_file} has been inserted into {db_name}')

if __name__ == '__main__':
    main()
