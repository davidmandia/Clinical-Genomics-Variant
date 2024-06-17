import re
import argparse
import os
import requests

def identify_short_indels(cigar_string):
    # Regular expression to extract length and operation from CIGAR string
    cigar_pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    
    # List to store short indels and mismatches
    short_indels = []
    position = 0
    
    # Iterate through each operation in the CIGAR string
    for length, op in re.findall(cigar_pattern, cigar_string):
        # Convert length to integer
        length = int(length)
        
        if op == 'I' or op == 'D':
            # Check if the length of the indel is less than or equal to 3
            if length < 3:
                # Append the operation, length, and position to the list
                short_indels.append((op, length, position))
        
        # Update the position for match/mismatch
        if op in 'M=X':
            position += length
        # Update the position to account for the indel or skipped region
        elif op in 'DN':
            position += length
    
    return short_indels



def get_ncbi_sequence(accession, start=None, end=None):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text"
    }
    
    # Add the sequence range if specified
    if start is not None and end is not None:
        params["seq_start"] = start
        params["seq_stop"] = end

    response = requests.get(base_url, params=params)

    if response.status_code == 200:
        return response.text
    else:
        response.raise_for_status()

def parse_sam(file_path):
    indels_obj = {}
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('@'):
                    # Skip header lines
                    continue
                fields = line.strip().split('\t')
                #print(fields)
                if len(fields) < 6:
                    raise ValueError("Invalid SAM format")
                read_id = fields[0]
                genome_ref = fields[2]
                pos = int(fields[3])
                cigar = fields[5]
                cigar_indels = identify_short_indels(cigar)
                
                # At least one short indel was found (op, length, position)
                if len(cigar_indels) > 0:
                    genomics_ref_sequence = get_ncbi_sequence(genome_ref, start=0, end=10)
                    #print(genomics_ref_sequence[cigar_indels[2] -3: cigar_indels[2] + 3])
                    print(genomics_ref_sequence)
                    indels_obj[read_id] = []
                    indels_obj[read_id].append(genome_ref)
                    indels_obj[read_id].append(cigar_indels)

                    
        return indels_obj
    except FileNotFoundError:
        print(f"Error: The file {file_path} does not exist.")
        return {}, {}
    except Exception as e:
        print(f"Error reading SAM file: {e}")
        return {}, {}   

def main():
    parser = argparse.ArgumentParser(description="Process a SAM file to identify short indels and mismatches.")
    parser.add_argument('sam_file', type=str, help="Path to the SAM file")
    
    args = parser.parse_args()
    
    # Check if the file has a .sam extension
    if not args.sam_file.lower().endswith('.sam'):
        print("Error: The provided file does not have a .sam extension.")
        return
    
    # Check if the file exists
    if not os.path.isfile(args.sam_file):
        print("Error: The provided file does not exist.")
        return
    
    indels_obj = parse_sam(args.sam_file)
    
    # Create the "outputs" directory if it doesn't exist
    output_dir = "outputs"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create output file name for indels
    output_file_name_indels = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(args.sam_file))[0]}_indels.txt")
    
    # Write the results to the output file for indels
    try:
        with open(output_file_name_indels, 'w') as output_file:
            for read_id, indels in indels_obj.items():
                output_file.write(f"{read_id}: {indels}\n")
        print(f"Results have been written to {output_file_name_indels}")
    except Exception as e:
        print(f"Error writing to output file: {e}")
    
    # Write the results to the output file for mismatches


if __name__ == "__main__":
    main()
