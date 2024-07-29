import pandas as pd
import sqlite3

# Load the exons data from the CSV file
exons_df = pd.read_csv('main/exons.csv')
# Connect to the SQLite database
conn = sqlite3.connect('output/database/GRCh38_indels_variant.db')
cursor = conn.cursor()

# Ensure the 'exon_id' column exists in the database
cursor.execute("PRAGMA table_info(variants)")
columns = [column[1] for column in cursor.fetchall()]
if 'exon_id' not in columns:
    cursor.execute("ALTER TABLE variants ADD COLUMN exon_id TEXT")

# Index the columns used for matching if not already indexed
cursor.execute("CREATE INDEX IF NOT EXISTS idx_variants_chrom_pos ON variants (chrom, pos)")

# Batch the updates
batch_size = 1000
updates = []

# Fetch existing variants from the database
variants_df = pd.read_sql_query("SELECT * FROM variants", conn)

# Update each variant with matching exon information
for index, variant in variants_df.iterrows():
    chrom = variant['chrom']
    pos = variant['pos']
    
    # Match on chromosome and position range
    matching_exons = exons_df[(exons_df['Chromosome'] == chrom) & (exons_df['Start'] <= pos) & (exons_df['End'] >= pos)]
    if not matching_exons.empty:
        exon_ids = ','.join(matching_exons['Exon ID'].unique())
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