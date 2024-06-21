import psycopg2
from psycopg2 import Error

endpoint = "variants.cjcmykm4g8ti.us-east-1.rds.amazonaws.com"
database = "variants"
user="postgres"
password="postgres"

# Function to connect to PostgreSQL database
def connect_to_database():
    try:
        connection = psycopg2.connect(
            user="your_username",
            password="your_password",
            host="your_postgresql_endpoint",  # Replace with your PostgreSQL endpoint
            port="5432",                      # Replace with your PostgreSQL port
            database="your_database_name"      # Replace with your PostgreSQL database name
        )
        return connection
    except (Exception, Error) as error:
        print("Error while connecting to PostgreSQL:", error)
        return None

def main():
    # Connect to the database
    connection = connect_to_database()

    if connection:
        print("Successfully connected to the database!")

        # Further database operations can be performed here

        # Remember to close the connection when done
        connection.close()
        print("Connection closed.")

if __name__ == "__main__":
    main()
