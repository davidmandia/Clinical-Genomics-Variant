import re

def modify_sqlite_to_postgresql(sqlite_sql_file, postgresql_sql_file):
    with open(sqlite_sql_file, 'r') as file:
        sql_content = file.read()
    
    # Replace AUTOINCREMENT with SERIAL
    sql_content = re.sub(r'\bAUTOINCREMENT\b', '', sql_content, flags=re.IGNORECASE)
    
    # Replace TEXT with VARCHAR or keep as TEXT in PostgreSQL
    # You can decide to convert TEXT to VARCHAR or keep as TEXT, below keeps as TEXT
    sql_content = re.sub(r'\bTEXT\b', 'VARCHAR', sql_content, flags=re.IGNORECASE)
    # or you can leave TEXT as it is in PostgreSQL since it's compatible
    # Remove or comment out PRAGMA statements and other SQLite specific syntax
    sql_content = re.sub(r'PRAGMA.*?;', '', sql_content, flags=re.IGNORECASE)
    sql_content = re.sub(r'--.*?\\n', '', sql_content, flags=re.IGNORECASE)  # Remove comments
    
    # Replace 'Na' with NULL (case-sensitive)
    sql_content = re.sub(r"'Na'", 'NULL', sql_content)
    sql_content = re.sub(r'\bCREATE TABLE variants\b', 'CREATE TABLE variant_GRCH37', sql_content, flags=re.IGNORECASE)
    
    # Remove or modify other SQLite specific syntax if necessary
    # This part can be expanded depending on specific requirements or known differences
    
    with open(postgresql_sql_file, 'w') as file:
        file.write(sql_content)
        


# Example usage
sqlite_sql_file = "data_transfer_AWS/GRCh37_indels_variant.sql" 
postgresql_sql_file = 'data_transfer_AWS/postgre_GRCh37_indels_variant.sql'

modify_sqlite_to_postgresql(sqlite_sql_file, postgresql_sql_file)
