#!/usr/bin/env python
import re
import argparse
import os
import subprocess

# Function to fetch the reference genome sequence using BLAST database
def get_sequence_blast_db(db, accession, start=None, end=None):
    # Construct the blastdbcmd command
    cmd = ["blastdbcmd", "-db", db, "-entry", accession]
    
    if start and end:
        cmd.extend(["-range", f"{start}-{end}"])
    
    try:
        # Execute the command and fetch the sequence
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        sequence = "".join(result.stdout.split("\n")[1:])  # Skip the first line which is the header
        return sequence
    except subprocess.CalledProcessError as e:
        raise Exception(f"Error fetching reference sequence from BLAST DB: {e}")

# Function to identify short indels and provide data formatted for VCF
def identify_short_indels(cigar_string):
    cigar_pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    short_indels = []
    genomic_position = 0
    transcriptomic_position = 0
    
    for length, op in re.findall(cigar_pattern, cigar_string):
        length = int(length)
        
        if op in ['I', 'D', 'N']:
            if length < 3:
                short_indels.append((op, length, genomic_position, transcriptomic_position, cigar_string))
        
        # Update the transcriptomic position
        if op in 'MIS=X':
            transcriptomic_position += length
        # Update the genomic position for match/mismatch
        if op in 'M=XDN':
            genomic_position += length
    
    return short_indels

# Function to extract chromosome number from the reference genome identifier
def extract_chromosome_number(genome_ref):
    match = re.search(r'NC_(\d+)', genome_ref)
    if match:
        return str(int(match.group(1)))  # Convert to integer and back to string to remove leading zeros
    return genome_ref

# Function to parse SAM file and format indels for VCF
def parse_sam(file_path, db, pseudo=False): # If pseudo == False we want the actual reference sequence 
    indels = []
    missing_sequences = []
    
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('@'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 6:
                    raise ValueError("Invalid SAM format")
                read_id = fields[0]
                genome_ref = fields[2]
                chrom = extract_chromosome_number(genome_ref)
                pos = int(fields[3])
                cigar = fields[5]
                cigar_indels = identify_short_indels(cigar)
                
                if len(cigar_indels) > 0:
                    for cigar_indel in cigar_indels:
                        op, length, genomic_pos_offset, transcriptomic_pos_offset, cigar_string = cigar_indel
                        genomic_pos = pos + genomic_pos_offset
                        transcript_sequence = fields[9]

                        try: # Will not run API call if pseudo is set to True - generate Pseudo VCF
                            if pseudo == True:
                                if op == 'D':
                                    ref_seq = "N"+"N"*length
                                    alt_seq = transcript_sequence[transcriptomic_pos_offset -1]
                                elif op == 'I':
                                    ref_seq = "N"
                                    alt_seq = ref_seq + transcript_sequence[transcriptomic_pos_offset-1:transcriptomic_pos_offset + length -1]
                            
                            elif pseudo == False: # Will use BLAST DB to fetch the reference sequence
                                if op == 'D':
                                    ref_seq = get_sequence_blast_db(db, genome_ref, genomic_pos, genomic_pos + length)
                                    alt_seq = ref_seq[0]  # Deletion
                                elif op == 'I':
                                    ref_seq = get_sequence_blast_db(db, genome_ref, genomic_pos, genomic_pos)
                                    alt_seq = ref_seq + transcript_sequence[transcriptomic_pos_offset-1:transcriptomic_pos_offset + length -1]

                                # Verify the sequences
                                if op == 'D':
                                    genomic_ref_check = get_sequence_blast_db(db, genome_ref, genomic_pos, genomic_pos + length)
                                    if ref_seq != genomic_ref_check:
                                        print(f"Error: Reference sequence does not match for {genome_ref} at position {genomic_pos}")
                                elif op == 'I':
                                    genomic_ref_check = get_sequence_blast_db(db, genome_ref, genomic_pos, genomic_pos)
                                    transcript_alt_check = transcript_sequence[transcriptomic_pos_offset-1:transcriptomic_pos_offset + length -1]
                                    if ref_seq != genomic_ref_check:
                                        print(f"Error: Reference sequence does not match for {genome_ref} at position {genomic_pos}")
                                    if alt_seq[1:] != transcript_alt_check:
                                        print(f"Error: ALT sequence does not match transcript sequence for {read_id} at transcript position {transcriptomic_pos_offset}")
                                
                            indels.append({
                                'CHROM': chrom,  # Use extracted chromosome number
                                'POS': genomic_pos,
                                'ID': '.',
                                'REF': ref_seq,
                                'ALT': alt_seq,
                                'QUAL': 99,
                                'FILTER': 'PASS',
                                'INFO': f'DP=100;LEN= {length} ;TYPE={"DEL" if op == "D" else "INS"};TRANSCRIPT= {read_id} ;TRANSCRIPT_POS= {transcriptomic_pos_offset} ;CIGAR: {cigar_string} ;GENOME_REF= {genome_ref}'
                            })
                        except Exception as e:
                            missing_sequences.append({
                                'TRANSCRIPT': read_id,
                                'GENOME_REF': genome_ref,
                                'TRANSCRIPT_POS': transcriptomic_pos_offset,
                                'GENOME_POS': genomic_pos,
                                'OP': op,
                                'LENGTH': length,
                                'CIGAR': cigar_string
                            })
        return indels, missing_sequences
    except FileNotFoundError:
        print(f"Error: The file {file_path} does not exist.")
        return [], []
    except Exception as e:
        print(f"Error reading SAM file: {e}")
        return [], []

# Function to write indels to a VCF file
def write_to_vcf(indels, output_file):
    with open(output_file, 'w') as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=myVariantCaller\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        for indel in indels:
            vcf.write(f"{indel['CHROM']}\t{indel['POS']}\t{indel['ID']}\t{indel['REF']}\t{indel['ALT']}\t{indel['QUAL']}\t{indel['FILTER']}\t{indel['INFO']}\n")

# Function to write missing sequences to a text file
def write_missing_sequences(missing_sequences, output_file):
    with open(output_file, 'w') as txt:
        txt.write("TRANSCRIPT\tGENOME_REF\tTRANSCRIPT_POS\tGENOME_POS\tOP\tLENGTH\tCIGAR\n")
        
        for seq in missing_sequences:
            txt.write(f"{seq['TRANSCRIPT']}\t{seq['GENOME_REF']}\t{seq['TRANSCRIPT_POS']}\t{seq['GENOME_POS']}\t{seq['OP']}\t{seq['LENGTH']}\t{seq['CIGAR']}\n")

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def main():
    parser = argparse.ArgumentParser(description="Process a SAM file to identify short indels and mismatches.")
    parser.add_argument('sam_file', type=str, help="Path to the SAM file")
    parser.add_argument('--db', type=str, help="Path to the BLAST database (not required if --pseudo is True)")
    parser.add_argument('--pseudo', type=str2bool, nargs='?', const=True, default=False, help="Provides or not the reference in the PseudoVCF")

    args = parser.parse_args()
    pseudo = args.pseudo

    # Check if the file has a .sam extension
    if not args.sam_file.lower().endswith('.sam'):
        print("Error: The provided file does not have a .sam extension.")
        return
    
    # Check if the file exists
    if not os.path.isfile(args.sam_file):
        print("Error: The provided file does not exist.")
        return
    
    # Ensure db argument is provided if pseudo is False
    if not pseudo and not args.db:
        print("Error: The --db argument is required unless --pseudo is set to True.")
        return
    
    output_dir = "output"
    vcfs_dir = os.path.join(output_dir, "VCFs")
    os.makedirs(vcfs_dir, exist_ok=True)

    # Parse the SAM file
    indels_obj, missing_sequences = parse_sam(args.sam_file, args.db, pseudo)
    
    # Write indels to a VCF file
    output_file_name_indels_vcf = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(args.sam_file))[0]}_{pseudo}_indels.vcf")
    write_to_vcf(indels_obj, output_file_name_indels_vcf)
    print(f"VCF file created: {output_file_name_indels_vcf}")
    
    # If there are missing sequences, write them to a text file
    if missing_sequences:
        output_file_name_missing_txt = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(args.sam_file))[0]}_missing_sequences.txt")
        write_missing_sequences(missing_sequences, output_file_name_missing_txt)
        print(f"Some sequences could not be fetched. Details saved to: {output_file_name_missing_txt}")

if __name__ == "__main__":
    main()
