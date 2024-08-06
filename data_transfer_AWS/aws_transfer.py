import psycopg2

sql_file_path = "data_transfer_AWS/postgre_GRCh38_indels_variant.sql"

# Connection details
host = "clinical-variant.cjcmykm4g8ti.us-east-1.rds.amazonaws.com"
port = "5432"
database = "postgres"  # Replace with your database name
user = "postgres"  # Replace with your username
password = "admin123"  # Replace with your password

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
    
    # Drop the table if it exists
    cursor.execute("DROP TABLE IF EXISTS variants;")
    
    # Commit the changes
    connection.commit()
    
    print("Table variants deleted successfully.")
    
    # Open and read the SQL file
    with open(sql_file_path, 'r', encoding='utf-8') as file:
        sql_content = file.read()
    
    # Split the SQL content by semicolon to handle multiple queries
    sql_commands = sql_content.split(';')

    # Wrap the execution in a transaction
    cursor.execute("BEGIN;")
    
    for command in sql_commands:
        if command.strip():
            try:
                cursor.execute(command)
            except Exception as e:
                print(f"Error executing command: {command}")
                print(f"Error message: {e}")
                cursor.execute("ROLLBACK;")
                raise e

    # Commit the transaction
    cursor.execute("COMMIT;")
    
    print("SQL file executed successfully.")
    
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
