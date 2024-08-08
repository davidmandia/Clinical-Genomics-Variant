import sqlite3
import os

def remove_duplicates(cursor):
    cursor.execute('''
        DELETE FROM variants
        WHERE rowid NOT IN (
            SELECT MIN(rowid)
            FROM variants
            GROUP BY chrom, pos, ref, alt
        )
    ''')

def query_database(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    remove_duplicates(cursor)

    # Count total rows
    cursor.execute("SELECT COUNT(*) FROM variants")
    total_rows = cursor.fetchone()[0]

    # Count rows where clinical label is different from "No Data"
    cursor.execute("SELECT COUNT(*) FROM variants WHERE clinical_label != 'No Data'")
    relevant_rows = cursor.fetchone()[0]

    # Calculate percentage
    relevant_percentage = (relevant_rows / total_rows) * 100 if total_rows > 0 else 0

    print(f"Total rows: {total_rows}")
    print(f"Rows with clinical relevance: {relevant_rows}")
    print(f"Percentage of rows with clinical relevance: {relevant_percentage:.2f}%")
    print("\n")

    #Execute a SELECT * query for rows where clinical label is different from "No Data"
    cursor.execute("SELECT * FROM variants  WHERE af < 0.1 LIMIT 10;")
    rows = cursor.fetchall()

    # Get column names
    column_names = [description[0] for description in cursor.description]

    # Print column names
    print("\t".join(column_names))
    print("-" * 100)  # Print a separator line

    # Print each row
    for row in rows:
        print("\t".join(str(value) if value is not None else "None" for value in row))

    conn.close()

# Specify the path to your database
db_name = "GRCh37_indels_variant.db"  # Replace with the database name you are querying
output_dir = "output"
dbs_dir = os.path.join(output_dir, "database")
db_path = os.path.join(dbs_dir, db_name)

query_database(db_path)
