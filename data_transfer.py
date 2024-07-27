import sqlite3
import pandas as pd
import psycopg2

# Step 1: Extract Data from SQLite
def extract_data_from_sqlite(sqlite_db_path):
    # Connect to SQLite database
    sqlite_conn = sqlite3.connect(sqlite_db_path)
    # Extract data from SQLite database
    df = pd.read_sql_query("SELECT * FROM variants", sqlite_conn)
    # Close the SQLite connection
    sqlite_conn.close()
    return df

# Step 2: Prepare PostgreSQL Database
def prepare_postgresql_database(pg_conn):
    cursor = pg_conn.cursor()
    # Create table if it doesn't exist
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS variants (
            chrom TEXT,
            pos INTEGER,
            ref TEXT,
            alt TEXT,
            qual REAL,
            filter TEXT,
            genomic_ref TEXT,
            operation TEXT,
            transcript_ref TEXT,
            transcript_pos TEXT,
            af REAL,
            af_eas REAL,
            af_nfe REAL,
            af_fin REAL,
            af_amr REAL,
            af_afr REAL,
            af_asj REAL,
            af_oth REAL,
            af_sas REAL,
            af_mid REAL,
            af_ami REAL,
            genes TEXT,
            consequences TEXT,
            clinically_relevant TEXT,
            clinical_label TEXT
        );
    """)
    pg_conn.commit()
    cursor.close()

# Step 3: Transfer Data to PostgreSQL
def transfer_data_to_postgresql(pg_conn, df):
    cursor = pg_conn.cursor()
    # Insert data into PostgreSQL
    for _, row in df.iterrows():
        cursor.execute("""
            INSERT INTO variants (
                chrom, pos, ref, alt, qual, filter, genomic_ref,
                operation, transcript_ref, transcript_pos,
                af, af_eas, af_nfe, af_fin,
                af_amr, af_afr, af_asj, af_oth,
                af_sas, af_mid, af_ami, genes, consequences,
                clinically_relevant, clinical_label
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """, tuple(row))
    # Commit the transaction
    pg_conn.commit()
    cursor.close()

def main():
    # Database paths and credentials
    sqlite_db_path = "output/database/GRCh38_indels_variant.db"

    pg_conn = psycopg2.connect(
    host="grch38-indels.cjcmykm4g8ti.us-east-1.rds.amazonaws.com",
    database="grch38_indels",
    user="postgres",
    password="Cavia2014!",
    port=5432
)

    
    # Extract data from SQLite
    df = extract_data_from_sqlite(sqlite_db_path)
    
    # Prepare PostgreSQL database
    prepare_postgresql_database(pg_conn)
    
    # Transfer data to PostgreSQL
    transfer_data_to_postgresql(pg_conn, df)
    
    # Close the PostgreSQL connection
    pg_conn.close()
    print("Data transfer complete.")

if __name__ == "__main__":
    main()
