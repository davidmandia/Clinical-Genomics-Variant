import re
import argparse
import os

def identify_short_indels_and_mismatches(cigar_string):
    # Regular expression to extract length and operation from CIGAR string
    cigar_pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    
    # List to store short indels and mismatches
    short_indels = []
    mismatches = []
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
        
        if op == 'X':
            # Add mismatches to the list
            mismatches.append((op, length, position))
        
        # Update the position for match/mismatch
        if op in 'M=X':
            position += length
        # Update the position to account for the indel or skipped region
        elif op in 'DN':
            position += length
    
    return short_indels, mismatches

def parse_sam(file_path):
    indels_obj = {}
    mismatches_obj = {}
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
                pos = int(fields[3])
                cigar = fields[5]
                cigar_indels, cigar_mismatches = identify_short_indels_and_mismatches(cigar)

                if len(cigar_indels) > 0:
                    indels_obj[read_id] = []
                    indels_obj[read_id].append(genome_ref)
                    indels_obj[read_id].append(cigar_indels)

                if len(cigar_mismatches) > 0:
                    mismatches_obj[read_id] = []
                    mismatches_obj[read_id].append(genome_ref)
                    mismatches_obj[read_id].append(cigar_mismatches)
                    
        return indels_obj, mismatches_obj
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
    
    indels_obj, mismatches_obj = parse_sam(args.sam_file)
    
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
    if mismatches_obj:
        output_file_name_mismatches = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(args.sam_file))[0]}_mismatches.txt")
        try:
            with open(output_file_name_mismatches, 'w') as output_file:
                for read_id, mismatches in mismatches_obj.items():
                    output_file.write(f"{read_id}: {mismatches}\n")
            print(f"Results have been written to {output_file_name_mismatches}")
        except Exception as e:
            print(f"Error writing to output file: {e}")

if __name__ == "__main__":
    main()
