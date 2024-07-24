import sqlite3
import pandas as pd
import argparse

def query_clinically_relevant_variants(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    query = """
    SELECT * FROM variants
    WHERE clinical_label IN ('Green', 'Amber', 'Red')
    LIMIT 10
    """

    cursor.execute(query)
    rows = cursor.fetchall()

    conn.close()
    return rows

def print_variants(variants):
    if variants:
        columns = ["chrom", "pos", "ref", "alt", "qual", "filter", "genomic_ref", "operation", 
                   "transcript_ref", "transcript_pos", "af", "af_eas", "af_nfe", "af_fin", 
                   "af_amr", "af_afr", "af_asj", "af_oth", "af_sas", "af_mid", "af_ami", 
                   "genes", "consequences", "clinically_relevant", "clinical_label"]
        
        df = pd.DataFrame(variants, columns=columns)
        print(df.head(10))
    else:
        print("No clinically relevant variants found.")

def main():
    parser = argparse.ArgumentParser(description="Print 10 lines from the database where the clinical relevance flag is green, amber, or red.")
    parser.add_argument('db_path', help="Path to the SQLite database")
    args = parser.parse_args()

    db_path = args.db_path

    variants = query_clinically_relevant_variants(db_path)
    print(variants)
    #print_variants(variants)

if __name__ == "__main__":
    main()
