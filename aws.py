import psycopg2

sql_file_path = "postgre_GRCh38_indels_variant.sql"

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
    
    # Execute the SQL commands
    cursor.execute(sql_content)
    
    # Commit the changes
    connection.commit()
    
    print("SQL file executed successfully.")
    

    # Execute a simple query
    cursor.execute("SELECT version();")
    
    cursor.execute("SELECT * FROM variants LIMIT 10;")

    
    
    
    # Fetch and print the result of the query
    record = cursor.fetchone()
    print("You are connected to - ", record, "\n")
    
    # Close the cursor and connection
    cursor.close()
    connection.close()
    
except Exception as error:
    print(f"Error while connecting to PostgreSQL: {error}")
