import psycopg2

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
    
    # Execute the query
    cursor.execute("SELECT * FROM variants LIMIT 10;")
    
    # Fetch the results
    rows = cursor.fetchall()
    
    # Print the results
    for row in rows:
        print(row)
    
    # Close the cursor and connection
    cursor.close()
    connection.close()
    
except psycopg2.DatabaseError as error:
    print(f"Error while connecting to PostgreSQL: {error}")
