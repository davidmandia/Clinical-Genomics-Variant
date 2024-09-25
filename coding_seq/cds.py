import sqlite3
import pandas as pd
import argparse
import os
import re

def extract_chromosome_number(genome_ref):
    """
    Extracts the chromosome number from a reference genome identifier.

    Args:
        genome_ref (str): The reference genome identifier.

    Returns:
        str: The chromosome number as a string.
    """
    match = re.search(r'NC_(\d+)', genome_ref)
    if match:
        return str(int(match.group(1)))  # Remove leading zeros
    return genome_ref

# Function to parse the GFF file and extract CDS information
def parse_gff(gff_file):
    """
    Parses a GFF file and extracts CDS information into a DataFrame.
    
    Args:
        gff_file (str): Path to the GFF file.

    Returns:
        pd.DataFrame: DataFrame containing CDS information with columns for chromosome, start, end, mRNA_id, cds, and tag.
    """
    cds_entries = []
    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if parts[2] == 'CDS':
                chrom = extract_chromosome_number(parts[0])
                start = int(parts[3])
                end = int(parts[4])
                info = parts[8]
                mRNA_id = None
                cds = None
                tag = None
                for entry in info.split(';'):
                    if entry.startswith("ID=cds"):
                        cds = entry.split('=cds-')[1].strip()
                    if entry.startswith('Parent='):
                        mRNA_id = entry.split('-')[1].strip()
                    if entry.startswith('tag'):
                        tag = entry.split('=')[1].strip()
                    else:
                        tag = None
                
                cds_entries.append({'chrom': chrom, 'start': start, 'end': end, "cds": cds, 'mRNA_id': mRNA_id, "tag": tag})
    
    return pd.DataFrame(cds_entries)

# Function to save the extracted CDS information to a CSV file
def save_to_csv(cds_df, output_path):
    cds_df.to_csv(output_path, index=False)

def main():
    parser = argparse.ArgumentParser(description="Extract CDS information from a GFF file and save to CSV.")
    parser.add_argument('gff_file', help="Path to the GFF file")
    parser.add_argument('output_file', help="name output file")
    
    args = parser.parse_args()
    
    gff_file = args.gff_file
    output = args.output_file
    output_folder = "coding_seq/results"
    
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Parse the GFF file to extract CDS information
    cds_df = parse_gff(gff_file)
    
    # Save the CDS DataFrame to CSV
    output_path = os.path.join(output_folder, f'{output}_cds.csv')
    save_to_csv(cds_df, output_path)
    
    print(f"CDS information saved to {output_path}")

if __name__ == "__main__":
    main()
