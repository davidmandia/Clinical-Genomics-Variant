import re
import argparse
import os
import requests
import time

# Function to fetch and cache the reference genome sequence
def get_ncbi_sequence(accession, start=None, end=None, retries=2):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text"
    }
    
    for attempt in range(retries):
        try:
            response = requests.get(base_url, params=params)
            response.raise_for_status()
            sequence = response.text.split('\n', 1)[1].replace('\n', '')
            if start and end:
                return sequence[start-1:end]
            return sequence
        except requests.exceptions.RequestException as e:
            if attempt < retries - 1:
                print(f"Error fetching reference sequence, retrying... ({attempt + 1}/{retries})")
                time.sleep(1)
            else:
                raise Exception(f"Failed to fetch sequence after {retries} retries: {e}")

# Cache the sequence
reference_sequences = {}

def get_sequence_cached(accession, start, end):
    if accession not in reference_sequences:
        reference_sequences[accession] = get_ncbi_sequence(accession)
    return reference_sequences[accession][start-1:end]

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
                short_indels.append((op, length, genomic_position, transcriptomic_position))
        
        if op in 'M=X':
            genomic_position += length
            transcriptomic_position += length
        elif op == 'D':
            genomic_position += length
        elif op == 'I':
            transcriptomic_position += length
        elif op == 'N':
            genomic_position += length
        elif op == 'S':
            transcriptomic_position += length
    
    return short_indels

# Function to parse SAM file and format indels for VCF
def parse_sam(file_path):
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
                pos = int(fields[3])
                cigar = fields[5]
                cigar_indels = identify_short_indels(cigar)
                
                if len(cigar_indels) > 0:
                    for cigar_indel in cigar_indels:
                        op, length, genomic_pos_offset, transcriptomic_pos_offset = cigar_indel
                        genomic_pos = pos + genomic_pos_offset
                        transcript_sequence = fields[9]
                        
                        try:
                            if op == 'D':
                                ref_seq = get_sequence_cached(genome_ref, genomic_pos, genomic_pos + length)
                                alt_seq = ref_seq[0]  # Deletion
                            elif op == 'I':
                                ref_seq = get_sequence_cached(genome_ref, genomic_pos, genomic_pos)
                                alt_seq = ref_seq + transcript_sequence[transcriptomic_pos_offset:transcriptomic_pos_offset + length]
                            
                            indels.append({
                                'CHROM': fields[2],
                                'POS': genomic_pos,
                                'ID': '.',
                                'REF': ref_seq,
                                'ALT': alt_seq,
                                'QUAL': 99,
                                'FILTER': 'PASS',
                                'INFO': f'DP=100;LEN={length};TYPE={"DEL" if op == "D" else "INS"};TRANSCRIPT={read_id};TRANSCRIPT_POS={transcriptomic_pos_offset}'
                            })
                        except Exception as e:
                            missing_sequences.append({
                                'TRANSCRIPT': read_id,
                                'GENOME_REF': genome_ref,
                                'TRANSCRIPT_POS': transcriptomic_pos_offset,
                                'GENOME_POS': genomic_pos,
                                'OP': op,
                                'LENGTH': length
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
        vcf.write("##reference=NC_000001.10\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        for indel in indels:
            vcf.write(f"{indel['CHROM']}\t{indel['POS']}\t{indel['ID']}\t{indel['REF']}\t{indel['ALT']}\t{indel['QUAL']}\t{indel['FILTER']}\t{indel['INFO']}\n")

# Function to write missing sequences to a text file
def write_missing_sequences(missing_sequences, output_file):
    with open(output_file, 'w') as txt:
        txt.write("TRANSCRIPT\tGENOME_REF\tTRANSCRIPT_POS\tGENOME_POS\tOP\tLENGTH\n")
        
        for seq in missing_sequences:
            txt.write(f"{seq['TRANSCRIPT']}\t{seq['GENOME_REF']}\t{seq['TRANSCRIPT_POS']}\t{seq['GENOME_POS']}\t{seq['OP']}\t{seq['LENGTH']}\n")

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
     
    # Create the "outputs" directory if it doesn't exist
    output_dir = "outputs"
    os.makedirs(output_dir, exist_ok=True)

    indels_obj, missing_sequences = parse_sam(args.sam_file)
    output_file_name_indels_VCF = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(args.sam_file))[0]}_indels.vcf")
    write_to_vcf(indels_obj, output_file_name_indels_VCF)
    print(f"VCF file created: {output_file_name_indels_VCF}")
    
    if missing_sequences:
        output_file_name_missing_txt = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(args.sam_file))[0]}_missing_sequences.txt")
        write_missing_sequences(missing_sequences, output_file_name_missing_txt)
        print(f"Some sequences could not be fetched. Details saved to: {output_file_name_missing_txt}")

if __name__ == "__main__":
    main()
