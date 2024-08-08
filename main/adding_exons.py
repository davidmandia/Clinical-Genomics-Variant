import sqlite3
import pandas as pd
import argparse

def update_exons_in_db(csv_file, db_file):
    # Load the exons data from the CSV file
    exons_df = pd.read_csv(csv_file)

    # Connect to the SQLite database
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # Ensure the 'exon_id' column exists in the database
    cursor.execute("PRAGMA table_info(variants)")
    columns = [column[1] for column in cursor.fetchall()]
    if 'exon_id' not in columns:
        cursor.execute("ALTER TABLE variants ADD COLUMN exon_id TEXT")

    # Index the columns used for matching if not already indexed
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_variants_chrom_pos ON variants (chrom, pos)")

    # Create the exons table and load data if not exists
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS exons (
            Chromosome TEXT,
            Start INTEGER,
            End INTEGER,
            Gene_ID TEXT,
            Transcript_ID TEXT,
            Exon_ID TEXT
        )
    """)

    # Load exons data into the database for faster access if not already there
    exons_df.to_sql('exons', conn, if_exists='replace', index=False)

    # Create index on the exons table
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_exons_chrom_start_end ON exons (Chromosome, Start, End)")

    # Batch the updates
    batch_size = 1000
    updates = []

    # Fetch existing variants from the database
    variants_df = pd.read_sql_query("SELECT chrom, pos FROM variants", conn)

    # Update each variant with matching exon information
    for index, variant in variants_df.iterrows():
        chrom = variant['chrom']
        pos = variant['pos']
        
        # Match on chromosome and position range
        cursor.execute("""
            SELECT Exon_ID FROM exons 
            WHERE Chromosome = ? AND Start <= ? AND End >= ?
        """, (chrom, pos, pos))
        
        matching_exons = cursor.fetchall()
        if matching_exons:
            exon_ids = ','.join([exon[0] for exon in matching_exons])
            updates.append((exon_ids, chrom, pos))
            
            # Commit in batches
            if len(updates) >= batch_size:
                cursor.executemany("""
                    UPDATE variants 
                    SET exon_id = ?
                    WHERE chrom = ? AND pos = ?
                """, updates)
                conn.commit()
                updates = []

    # Commit remaining updates
    if updates:
        cursor.executemany("""
            UPDATE variants 
            SET exon_id = ?
            WHERE chrom = ? AND pos = ?
        """, updates)

    conn.commit()
    conn.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update SQLite database with exon information from a CSV file.")
    parser.add_argument('db_file', help="Path to the SQLite database file")

    parser.add_argument('csv_file', help="Path to the input CSV file containing exon data")
    args = parser.parse_args()

    update_exons_in_db(args.csv_file, args.db_file)
