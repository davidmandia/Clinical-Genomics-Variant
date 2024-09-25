import pandas as pd 

import sqlite3 
import argparse
 

#  Load the database of transcripts 

def load_transcripts_from_db(db_path): 

    conn = sqlite3.connect(db_path) 

    query = "SELECT DISTINCT transcript_ref FROM variants"  

    df_transcripts = pd.read_sql_query(query, conn) 
    
    transc_list = df_transcripts['transcript_ref'].to_list()

    conn.close() 

    return transc_list

 

# Load the CSV file 

def load_transcripts_from_csv(sam_path): 

    df_sam = pd.read_csv(sam_path) 

    return df_sam

 

# Identify transcript 

def get_latest_transcripts(df_sam): 
    #print("sam", df_sam["transcript_ref"])
    sam_transc_latest = df_sam['transcript_ref'].tolist()
    
    return sam_transc_latest

 

# Step 4: Count how many transcripts from the database are in the latest version of the CSV 

def count_matching_transcripts(db_transcripts, latest_transcripts): 

    matching_count = sum(1 for transcript in db_transcripts if transcript in latest_transcripts) 

    return matching_count 

 

# Main function 

def main(db_path, sam_path):
       
    db_transcripts = load_transcripts_from_db(db_path) 

    df_sam = load_transcripts_from_csv(sam_path) 

    latest_transcripts = get_latest_transcripts(df_sam) 

    matching_count = count_matching_transcripts(db_transcripts, latest_transcripts) 

     

    print(f"Number of transcripts in the database that are also in the latest version from the CSV file: {matching_count}") 


parser = argparse.ArgumentParser(description="count the number of latest version in the database")

parser.add_argument('db', type=str, help="Path to the BLAST database") 

parser.add_argument('csv_from_sam_transc', type=str, help="Path to the SAM file")

args = parser.parse_args()

db_path = args.db   # Path to  SQLite database 

sam_path = args.csv_from_sam_transc  # Path to  CSV file 
main(db_path, sam_path) 