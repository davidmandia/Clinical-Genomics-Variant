import psycopg2
import os

# Database connection details
host = "clinical-variant.cjcmykm4g8ti.us-east-1.rds.amazonaws.com"
port = "5432"
database = "postgres"  # Replace with your database name
user = "postgres"  # Replace with your username
password = "admin123"  # Replace with your password

def delete_table():
    try:
        # Connect to your postgres DB
        conn = psycopg2.connect(
            host=host,
            port=port,
            database=database,
            user=user,
            password=password
        )

        # Open a cursor to perform database operations
        cursor = conn.cursor()

        # Drop the table variant_grch37 if it exists
        drop_table_query = "DROP TABLE IF EXISTS variants"
        cursor.execute(drop_table_query)
        print("Table variants deleted successfully")

        # Commit the changes to the database
        conn.commit()

        # Close the cursor and connection
        cursor.close()
        conn.close()

    except Exception as error:
        print(f"Error deleting table: {error}")

# Run the delete_table function
delete_table()
