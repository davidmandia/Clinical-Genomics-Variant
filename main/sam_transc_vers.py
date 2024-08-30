import os 

import pandas as pd 

import re 

import argparse 

 

def read_sam_file(sam_file): 

    """ 

    Reads a SAM file and extracts transcript information. 

 

    Args: 

        sam_file (str): Path to the SAM file. 

 

    Returns: 

        DataFrame: A DataFrame containing transcript references and related information. 

    """ 

    transcripts = [] 

     

    with open(sam_file, 'r') as file: 

        for line in file: 

            if line.startswith('@'): 

                continue  # Skip header lines 

            fields = line.split() 

            transcript_ref = fields[0] 

            transcripts.append(transcript_ref) 

     

    # Convert to DataFrame 

    df = pd.DataFrame(transcripts, columns=['transcript_ref']) 

     

    return df 

 

def get_highest_transcript_versions(df, output_dir, sam_name): 

    """ 

    Extracts transcript base and version, and identifies the highest version for each transcript base. 

 

    Args: 

        df (DataFrame): DataFrame containing transcript references. 

        output_dir (str): Directory to save the output file. 

        sam_name (str): Name of the SAM file (used for naming the output file). 

 

    Returns: 

        DataFrame: DataFrame containing the highest transcript versions. 

    """ 

    # Extract transcript base and version 

    df['transcript_base'] = df['transcript_ref'].str.extract(r'(\w+_\d+)') 

    df['transcript_version'] = df['transcript_ref'].str.extract(r'\.(\d+)$').astype(int) 

     

    # Identify the highest version for each transcript base 

    idx = df.groupby('transcript_base')['transcript_version'].idxmax() 

    highest_versions_df = df.loc[idx] 

 

    output_path = os.path.join(output_dir, f'highest_transcript_versions_{sam_name}.csv') 

    highest_versions_df.to_csv(output_path, index=False) 

    print(f"Highest transcript versions saved to {output_path}") 

     

    return highest_versions_df 

 

def main(): 

    parser = argparse.ArgumentParser(description="Extract and identify the highest transcript versions from a SAM file.") 

    parser.add_argument('sam_file', help="Path to the input SAM file") 

    args = parser.parse_args() 

 

    sam_file = args.sam_file 

    sam_name = os.path.basename(sam_file).replace('.sam', '') 

     

    output_dir = "analysis/results" 

    os.makedirs(output_dir, exist_ok=True) 

     

    # Read the SAM file 

    df = read_sam_file(sam_file) 

     

    # Get the highest transcript versions 

    highest_versions_df = get_highest_transcript_versions(df, output_dir, sam_name) 

     

    # Print the highest versions 

    print("Rows with the highest transcript versions:") 

    print(highest_versions_df) 

 

if __name__ == "__main__": 

    main() 