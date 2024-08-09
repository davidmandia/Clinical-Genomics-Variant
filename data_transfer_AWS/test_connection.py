import psycopg2
import os

def fetch_and_print_table_data(connection, table_name, limit=5):
    try:
        # Create a cursor object
        cursor = connection.cursor()
        
        # Query to fetch the first 'limit' rows from the specified table
        query = f"SELECT * FROM {table_name} LIMIT {limit};"
        cursor.execute(query)
        
        # Fetch the results
        rows = cursor.fetchall()
        
        # Print the column names
        column_names = [desc[0] for desc in cursor.description]
        print(f"\nColumns in {table_name}: {', '.join(column_names)}")
        
        # Print the rows
        print(f"First {limit} rows from {table_name}:")
        for row in rows:
            print(row)
        
        # Close the cursor
        cursor.close()
        
    except Exception as error:
        print(f"Error fetching data from {table_name}: {error}")

try:
    # Connection details from environment variables
    host = os.getenv("DB_HOST", "clinical-variant.cjcmykm4g8ti.us-east-1.rds.amazonaws.com")
    port = os.getenv("DB_PORT", "5432")
    database = os.getenv("DB_NAME", "postgres")
    user = os.getenv("DB_USER", "postgres")
    password = os.getenv("DB_PASSWORD", "admin123")
    
    # Establishing the connection
    connection = psycopg2.connect(
        host=host,
        port=port,
        database=database,
        user=user,
        password=password
    )
    
    print("Connected to the PostgreSQL database.")

    # Fetch and print data from 'variants' table
    fetch_and_print_table_data(connection, "variants")
    
    # Fetch and print data from 'variant_GRCH37' table
    fetch_and_print_table_data(connection, "variant_GRCH37")
    
    # Close the connection
    connection.close()
    print("Database connection closed.")
    
except Exception as error:
    print(f"Error while connecting to PostgreSQL: {error}")
