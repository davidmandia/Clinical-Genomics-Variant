import re
import os
import argparse

def modify_sqlite_to_postgresql(sqlite_sql_file, postgresql_sql_file):
    with open(sqlite_sql_file, 'r') as file:
        sql_content = file.read()
    
    # Replace AUTOINCREMENT with SERIAL
    sql_content = re.sub(r'\bAUTOINCREMENT\b', '', sql_content, flags=re.IGNORECASE)
    
    # Replace TEXT with VARCHAR or keep as TEXT in PostgreSQL
    sql_content = re.sub(r'\bTEXT\b', 'VARCHAR', sql_content, flags=re.IGNORECASE)
    
    # Remove or comment out PRAGMA statements and other SQLite specific syntax
    sql_content = re.sub(r'PRAGMA.*?;', '', sql_content, flags=re.IGNORECASE)
    sql_content = re.sub(r'--.*?\n', '', sql_content, flags=re.IGNORECASE)  # Remove comments
    
    # Replace 'Na' with NULL (case-sensitive)
    sql_content = re.sub(r"'Na'", 'NULL', sql_content)
    
    # For GRCh37
    if 'GRCh37' in sqlite_sql_file:
        sql_content = re.sub(r'\bCREATE TABLE variants\b', 'CREATE TABLE variant_GRCH37', sql_content, flags=re.IGNORECASE)
    
    with open(postgresql_sql_file, 'w') as file:
        file.write(sql_content)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert SQLite SQL file to PostgreSQL compatible SQL file.")
    parser.add_argument('sqlite_sql_file', help="Path to the input SQLite SQL file")
    args = parser.parse_args()

    sqlite_sql_file = args.sqlite_sql_file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    input_filename = os.path.basename(sqlite_sql_file)
    postgresql_sql_file = os.path.join(script_dir, f"postgre_{input_filename}")

    modify_sqlite_to_postgresql(sqlite_sql_file, postgresql_sql_file)
    print(f"Converted SQL file saved as {postgresql_sql_file}")
