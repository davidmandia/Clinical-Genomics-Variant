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
    
    # Open and read the SQL file
    with open(sql_file_path, 'r', encoding='utf-8') as file:
        sql_content = file.read()
    
    # Split the SQL content by semicolon to handle multiple queries
    sql_commands = sql_content.split(';')
    
    # Execute each SQL command separately
    for command in sql_commands:
        if command.strip():
            cursor.execute(command)
    
    # Commit the changes
    connection.commit()
    
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
