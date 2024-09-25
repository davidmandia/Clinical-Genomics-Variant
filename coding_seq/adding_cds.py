import sqlite3
import pandas as pd
import argparse


def update_cds_in_db(csv_file, db_file):
    # Load the cds data from the CSV file
    cds_df = pd.read_csv(csv_file)

    # Connect to the SQLite database
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # Ensure the 'excds_info' column exists in the database
    cursor.execute("PRAGMA table_info(variants)")
    columns = [column[1] for column in cursor.fetchall()]
    if 'cds_info' not in columns:
        cursor.execute("ALTER TABLE variants ADD COLUMN cds_info TEXT")
    
    if 'cds_tag' not in columns:
        cursor.execute("ALTER TABLE variants ADD COLUMN cds_tag TEXT")

    # Index the columns used for matching if not already indexed
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_variants_chrom_pos ON variants (chrom, pos)")

    # Create the cds table and load data if not exists
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS cds (
            chrom TEXT,
            start INTEGER,
            end INTEGER,
            cds TEXT,
            mRNA_id TEXT,
            tag TEXT
        )
    """)

    # Load cds data into the database for faster access if not already there
    cds_df.to_sql('cds', conn, if_exists='replace', index=False)

    # Create index on the cds table
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_cds_chrom_start_end ON cds (chrom, start, end)")

    # Batch the updates
    batch_size = 1000
    updates = []

    # Fetch existing variants from the database
    variants_df = pd.read_sql_query("SELECT chrom, pos, transcript_ref FROM variants", conn)

    # Update each variant with matching cds information
    for index, variant in variants_df.iterrows():
        chrom = variant['chrom']
        pos = variant['pos']
        transcipt_id = variant['transcript_ref']
        
        # Match on chromosome and position range
        cursor.execute("""
            SELECT cds, tag FROM cds 
            WHERE chrom = ? AND start <= ? AND end >= ? AND mRNA_id = ?
        """, (chrom, pos, pos, transcipt_id))
        
        matching_cds = cursor.fetchall()
        #print("matchin_cds", matching_cds)
        if matching_cds:
            matching_cds_cds = matching_cds[0][0]
           # print("matching_cds_cds", matching_cds_cds)
            matching_cds_tag = matching_cds[0][1] if len(matching_cds[0]) > 1 else "Na"
            #print("matching_cds_tag", matching_cds_tag)
            
            updates.append((matching_cds_cds, matching_cds_tag, chrom, pos, transcipt_id))
            
            # Commit in batches
          #  print("len updates", len(updates))  
            if len(updates) >= batch_size:
               # print("updates", updates)
                cursor.executemany("""
                    UPDATE variants 
                    SET cds_info = ?, cds_tag = ?
                    WHERE chrom = ? AND pos = ? AND transcript_ref = ?
                """, updates)
                conn.commit()
                updates = []

            # Commit remaining updates
            if updates:
                #print("updates_1", updates)
                cursor.executemany("""
                    UPDATE variants 
                    SET cds_info = ?, cds_tag = ?
                    WHERE chrom = ? AND pos = ? AND transcript_ref = ?
                """, updates)

    conn.commit()
    conn.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update SQLite database with cds information from a CSV file.")
    parser.add_argument('db_file', help="Path to the SQLite database file")

    parser.add_argument('csv_file', help="Path to the input CSV file containing cds data")
    args = parser.parse_args()

    update_cds_in_db(args.csv_file, args.db_file)
    print("CDS information updated successfully.")

