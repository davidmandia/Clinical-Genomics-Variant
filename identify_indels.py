import re
import argparse
import os

def identify_short_indels(cigar_string):
    # Regular expression to extract length and operation from CIGAR string
    cigar_pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    
    # List to store short indels
    short_indels = []
    
    # Iterate through each operation in the CIGAR string
    for length, op in re.findall(cigar_pattern, cigar_string):
        # Check if the operation is an insertion (I) or deletion (D)
        if op == 'I' or op == 'D':
            # Convert length to integer
            length = int(length)
            
            # Check if the length of the indel is less than or equal to 3
            if length < 3:
                # Append the operation and length to the short indels list
                short_indels.append((op, length))
    
    return short_indels

def parse_sam(file_path):
    alignments = {}
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('@'):
                    # Skip header lines
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 6:
                    raise ValueError("Invalid SAM format")
                read_id = fields[0]
                genome_ref = fields[2]
                cigar = fields[5]
                cigar_indel = identify_short_indels(cigar)
                if len(cigar_indel) > 0:
                    alignments[read_id] = []
                    alignments[read_id].append(cigar_indel)
                    
        return alignments
    except Exception as e:
        print(f"Error reading SAM file: {e}")
        return {}

def main():
    parser = argparse.ArgumentParser(description="Process a SAM file to identify short indels.")
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
    
    alignments = parse_sam(args.sam_file)
    
    # Create the "outputs" directory if it doesn't exist
    output_dir = "outputs"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create output file name
    output_file_name = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(args.sam_file))[0]}_indels.txt")
    
    # Write the results to the output file
    try:
        with open(output_file_name, 'w') as output_file:
            for read_id, indels in alignments.items():
                output_file.write(f"{read_id}: {indels}\n")
        print(f"Results have been written to {output_file_name}")
    except Exception as e:
        print(f"Error writing to output file: {e}")

if __name__ == "__main__":
    main()
