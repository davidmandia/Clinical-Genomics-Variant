import pandas as pd
import sqlite3
import argparse

def ensure_latest_column(conn, table_name, column_name):
    """ Ensure that the specified column exists in the table. """
    cursor = conn.cursor()
    
    # Check if the column exists
    cursor.execute(f"PRAGMA table_info({table_name})")
    columns = [col[1] for col in cursor.fetchall()]
    
    if column_name not in columns:
        # If the column does not exist, create it
        cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {column_name} TEXT")
        print(f"Column '{column_name}' added to table '{table_name}'.")

def mark_latest_transcripts(db_path, csv_file):
    """ 
    Marks transcripts as latest or not in a SQLite database based on a CSV file.
    
    Args:
        db_path (str): Path to the SQLite database.
        csv_file (str): Path to the CSV file with the latest version from SAM.
    """
    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)

    # Ensure the 'transcript_is_latest' column exists
    ensure_latest_column(conn, 'variants', 'transcript_is_latest')

    # Load the CSV file into a DataFrame
    df = pd.read_csv(csv_file)
    latest_versions = df["transcript_ref"].tolist()

    # Reset all entries in transcript_is_latest to 'No'
    with conn:
        conn.execute("UPDATE variants SET transcript_is_latest = 'No'")

    # Update the is_latest column for each transcript reference
    with conn:
        for tran in latest_versions:
            conn.execute("""
                UPDATE variants 
                SET transcript_is_latest = 'Yes' 
                WHERE transcript_ref = ?;
            """, (tran,))
    
    print("Transcript versions updated in the database.")
    
    # Close the connection
    conn.close()

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Update transcript versions in a SQLite database.")
    parser.add_argument('db_path', help="Path to the SQLite database.")
    parser.add_argument('csv_file', help="Path to the CSV file containing transcript references and versions.")

    args = parser.parse_args()

    # Run the marking function
    mark_latest_transcripts(args.db_path, args.csv_file)

if __name__ == "__main__":
    main()
    print("Script completed successfully.")
