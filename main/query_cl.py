import sqlite3
import pandas as pd
import argparse

def query_clinically_relevant_variants(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    query = """
    SELECT count(*) FROM variants
    WHERE clinical_label IN ('Green', 'Amber', 'Red');
    """

    cursor.execute(query)
    rows = cursor.fetchall()
    
    query = """
    SELECT count(*) FROM variants;
    """

    cursor.execute(query)
    rows2 = cursor.fetchall()
    
    query = """
    SELECT count(*) FROM variants 
    WHERE clinical_label IN ('Green', 'Amber', 'Red') and af != "Na" and gene = "MAPT";
    """

    cursor.execute(query)
    rows3 = cursor.fetchall()

    conn.close()
    return rows, rows2, rows3

def query_db(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    
    
    query = """
    SELECT count(*) FROM variants 
    WHERE clinical_label IN ('Green')  and gene_symbol = "MAPT";
    """

    cursor.execute(query)
    rows3 = cursor.fetchall()

    conn.close()
    return  rows3

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

    #variants = query_clinically_relevant_variants(db_path)
    variants = query_db(db_path)
    print(variants)
    #print_variants(variants)

if __name__ == "__main__":
    main()
