import psycopg2
import argparse
import os
import re
import configparser

def execute_sql_file(cursor, sql_file_path):
    with open(sql_file_path, 'r', encoding='utf-8') as file:
        sql_content = file.read()

    # Split the SQL content by semicolon to handle multiple queries
    sql_commands = sql_content.split(';')

    for command in sql_commands:
        if command.strip():
            try:
                cursor.execute(command)
            except psycopg2.Error as e:
                if 'already exists' in e.pgerror:
                    table_name = extract_table_name(command)
                    if table_name:
                        print(f"Table {table_name} already exists. Dropping and recreating.")
                        cursor.execute(f"DROP TABLE IF EXISTS {table_name} CASCADE;")
                        cursor.execute(command)
                else:
                    print(f"Error executing command: {command}")
                    print(f"Error message: {e}")
                    cursor.execute("ROLLBACK;")
                    raise e

def extract_table_name(command):
    match = re.search(r'CREATE TABLE\s+(\w+)', command, re.IGNORECASE)
    if match:
        return match.group(1)
    return None

def main():
    parser = argparse.ArgumentParser(description="Execute a SQL file on a PostgreSQL database.")
    parser.add_argument('sql_file_path', help="Path to the SQL file")
    args = parser.parse_args()

    sql_file_path = args.sql_file_path

    # Read database configuration from the configuration file
    config = configparser.ConfigParser()
    config.read('data_transfer_AWS/db_config.ini')

    host = config['postgresql']['host']
    port = config['postgresql']['port']
    database = config['postgresql']['database']
    user = config['postgresql']['user']
    password = config['postgresql']['password']

    try:
        # Establishing the connection
        connection = psycopg2.connect(
            host=host,
            port=port,
            database=database,
            user=user,
            password=password
        )
        
        # Create a cursor object
        cursor = connection.cursor()

        # Print a message indicating that the connection was successful
        print("Connected to the database.")

        try:
            # Wrap the execution in a transaction
            cursor.execute("BEGIN;")
            
            # Execute the SQL file
            execute_sql_file(cursor, sql_file_path)

            # Commit the transaction
            cursor.execute("COMMIT;")
            
            print("SQL file executed successfully.")

        except Exception as e:
            # If any exception occurs, rollback the transaction
            connection.rollback()
            print(f"Error during transaction: {e}")

        # Execute a simple query to verify the table creation
        cursor.execute("SELECT version();")
        version_record = cursor.fetchone()
        print("PostgreSQL version - ", version_record, "\n")
        
        # Query to check the new table
        cursor.execute("SELECT * FROM variants LIMIT 10;")
        records = cursor.fetchall()
        print("First 10 records from variants table:")
        for record in records:
            print(record)
        
        # Close the cursor and connection
        cursor.close()
        connection.close()
        
    except Exception as error:
        print(f"Error while connecting to PostgreSQL: {error}")

if __name__ == "__main__":
    main()
