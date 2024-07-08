#!/usr/bin/env python
import re
import argparse
import os
import subprocess
import time

# Global variable to store cached sequences
accession_cache = {}

# Function to fetch the reference genome sequence using BLAST database
def get_sequence_blast_db(db, accession, start=None, end=None, sleep_time=2):
    global accession_cache

    # Check if the accession sequence is already cached
    if accession in accession_cache:
        full_sequence = accession_cache[accession]
    else:
        # Construct the blastdbcmd command to fetch the full sequence for the accession
        cmd = ["blastdbcmd", "-db", db, "-entry", accession]

        try:
            # Execute the command and fetch the sequence
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            full_sequence = "".join(result.stdout.split("\n")[1:])  # Skip the first line which is the header
            accession_cache[accession] = full_sequence  # Cache the full sequence
            
            # Sleep for the specified amount of time
            time.sleep(sleep_time)
        except subprocess.CalledProcessError as e:
            raise Exception(f"Error fetching reference sequence from BLAST DB: {e}")

    # If start and end are specified, extract the precise region
    if start is not None and end is not None:
        try:
            # Full Sequence is a String so python 0_based index
            sequence = full_sequence[start:end+1]
            return sequence
        except IndexError:
            raise Exception(f"Invalid range: {start}-{end} is out of bounds for accession {accession}")
    else:
        return full_sequence

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
               # print("indels", (op, length, genomic_position, transcriptomic_position))
                short_indels.append((op, length, genomic_position, transcriptomic_position))

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
def parse_sam(file_path, db, pseudo=False, sleep_time=1):  # If pseudo == False we want the actual reference sequence
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
                cigar_string = fields[5]
                cigar_indels = identify_short_indels(cigar_string)

                if len(cigar_indels) > 0:
                    for cigar_indel in cigar_indels:
                        #print(cigar_indel)
                        op, length, genomic_pos_offset, transcriptomic_pos_offset = cigar_indel
                        # -2 because 1 for the staring position in the sam file + genomi position and the other 0_indexing
                        genomic_pos_0_based = pos + genomic_pos_offset - 2  # 0-based position
                        transcript_sequence = fields[9]
                        transcriptomic_pos_offset_0_based = transcriptomic_pos_offset - 1

                        try:
                            if pseudo:
                                if op == 'D':
                                    ref_seq = transcript_sequence[transcriptomic_pos_offset_0_based] + "N" * length
                                    alt_seq = transcript_sequence[transcriptomic_pos_offset_0_based]

                                elif op == 'I':
                                    ref_seq = "N"
                                    alt_seq = ref_seq + transcript_sequence[transcriptomic_pos_offset_0_based:transcriptomic_pos_offset_0_based + length]
                            elif pseudo == False:
                                # Fetch reference sequence using cached or direct BLAST DB query
                                if op == 'D':
                                   # print("deletion")
                                    ref_seq = get_sequence_blast_db(db, genome_ref, genomic_pos_0_based, genomic_pos_0_based + length, sleep_time)
                                    alt_seq = transcript_sequence[transcriptomic_pos_offset_0_based]  # Deletion, alt_seq should be empty or a placeholder

                                elif op == 'I':
                                    #print("insertion")
                                    ref_seq = get_sequence_blast_db(db, genome_ref, genomic_pos_0_based, genomic_pos_0_based, sleep_time)  # Fetch reference sequence
                                    alt_seq = transcript_sequence[transcriptomic_pos_offset_0_based:transcriptomic_pos_offset_0_based + length + 1]  # (python slicing)

                            # Append indel record to the list
                            indels.append({
                                'CHROM': chrom,
                                'POS': genomic_pos_0_based + 1,  # Convert to 1-based position for VCF
                                'ID': '.',
                                'REF': ref_seq,
                                'ALT': alt_seq,
                                'QUAL': 99,
                                'FILTER': 'PASS',
                                'INFO': f'DP=100;LEN={length};TYPE={"DEL" if op == "D" else "INS"};TRANSCRIPT= {read_id} ;TRANSCRIPT_POS= {transcriptomic_pos_offset} ;GENOME_REF= {genome_ref} ' # for testing i have added the CIGAR on the output INFO "CIGAR: {cigar_string}" 
                            })

                        except Exception as e:
                            missing_sequences.append({
                                'TRANSCRIPT': read_id,
                                'GENOME_REF': genome_ref,
                                'TRANSCRIPT_POS': transcriptomic_pos_offset,
                                'GENOME_POS': genomic_pos_0_based,
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
        
        #write VCF header 
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=myVariantCaller\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        vcf.write('##contig=<ID=1>\n')
        vcf.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
        vcf.write('##INFO=<ID=LEN,Number=1,Type=Integer,Description="Variant Length">\n')
        vcf.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant Type">\n')
        vcf.write('##INFO=<ID=TRANSCRIPT,Number=1,Type=String,Description="Transcript ID">\n')
        vcf.write('##INFO=<ID=TRANSCRIPT_POS,Number=1,Type=Integer,Description="Position in Transcript">\n')
        vcf.write('##INFO=<ID=CIGAR,Number=1,Type=String,Description="CIGAR string">\n')
        vcf.write('##INFO=<ID=GENOME_REF,Number=1,Type=String,Description="Genome Reference">\n')
        

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
    parser.add_argument('--db', type=str, help="Path to the BLAST database")
    parser.add_argument('--pseudo', type=str2bool, nargs='?', const=True, default=False, help="Generate pseudo VCF")
    parser.add_argument('--sleep', type=float, default=1, help="Time to sleep between requests (seconds)")

    args = parser.parse_args()
    pseudo = args.pseudo
    sleep_time = args.sleep

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
    indels_obj, missing_sequences = parse_sam(args.sam_file, args.db, pseudo, sleep_time)

    # Write indels to a VCF file
    output_file_name_indels_vcf = os.path.join(vcfs_dir, f"{os.path.splitext(os.path.basename(args.sam_file))[0]}_{pseudo}_indels.vcf")
    write_to_vcf(indels_obj, output_file_name_indels_vcf)
    print(f"VCF file created: {output_file_name_indels_vcf}")

    # If there are missing sequences, write them to a text file
    if missing_sequences:
        output_file_name_missing_txt = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(args.sam_file))[0]}_missing_sequences.txt")
        write_missing_sequences(missing_sequences, output_file_name_missing_txt)
        print(f"Some sequences could not be fetched. Details saved to: {output_file_name_missing_txt}")

if __name__ == "__main__":
    main()
