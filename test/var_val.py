import sys
import requests
import pandas as pd
from time import sleep

def get_vep_data_batch(hgvs_list, batch_size=200):
    """Fetch data from Ensembl VEP REST API in batches"""
    server = "https://rest.ensembl.org"
    ext = "/vep/human/hgvs"
    
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }
    
    all_results = []
    
    for i in range(0, len(hgvs_list), batch_size):
        batch = hgvs_list[i:i+batch_size]
        data = {"hgvs_notations": batch}
        
        response = requests.post(server + ext, headers=headers, json=data)
        
        if response.status_code != 200:
            print(f"Failed to retrieve data for batch {i//batch_size + 1}: {response.status_code}")
            continue
        
        all_results.extend(response.json())
        
        # Be nice to the API
        sleep(1)
    
    return all_results

def parse_vep_data(data):
    """Parse VEP API response for allele frequencies and SIFT score"""
    results = []
    for variant in data:
        hgvs = variant.get('input', '')
        frequencies = {}
        sift_score = None
        
        if 'colocated_variants' in variant:
            for cv in variant['colocated_variants']:
                if 'frequencies' in cv:
                    frequencies = cv['frequencies']
                    break
        
        if 'transcript_consequences' in variant:
            for tc in variant['transcript_consequences']:
                if 'sift_score' in tc:
                    sift_score = tc['sift_score']
                    break
        
        results.append({
            'HGVS': hgvs,
            'AF': ','.join([f"{pop}:{freq}" for pop, freq in frequencies.items()]) if frequencies else '',
            'SIFT': sift_score if sift_score is not None else ''
        })
    
    return results

def process_file(input_file):
    """Process the input file and return unique HGVS values"""
    unique_values = set()
    with open(input_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) > 11:
                    hgvs = fields[11]  # HGVS_Genomic_GRCh38 is at index 11
                    unique_values.add(hgvs)
    return list(unique_values)

def main(input_file, output_file):
    # Get unique HGVS values
    unique_hgvs = process_file(input_file)
    print(f"Found {len(unique_hgvs)} unique HGVS values")

    # Get VEP data
    vep_data = get_vep_data_batch(unique_hgvs)
    print("VEP data retrieved")

    # Parse VEP data
    parsed_data = parse_vep_data(vep_data)
    print("VEP data parsed")

    # Create a DataFrame with VEP data
    df_vep = pd.DataFrame(parsed_data)
    df_vep.to_csv(output_file, sep='\t', index=False)

    # # Read the original file into a DataFrame
    # df = pd.read_csv(input_file, sep='\t', comment='#')

    # # Merge the original DataFrame with VEP data
    # df_merged = df.merge(df_vep, left_on='HGVS_Genomic_GRCh38', right_on='HGVS', how='left')

    # # Drop the redundant 'HGVS' column
    # df_merged = df_merged.drop('HGVS', axis=1)

    # # Write the merged DataFrame to the output file
    # df_merged.to_csv(output_file, sep='\t', index=False)
    print(f"Results written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)