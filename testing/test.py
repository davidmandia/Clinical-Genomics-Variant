import re
import os

def parse_vcf(file_path):
    indels = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            indel = {
                'CHROM': fields[0],
                'POS': int(fields[1]),
                'ID': fields[2],
                'REF': fields[3],
                'ALT': fields[4],
                'QUAL': fields[5],
                'FILTER': fields[6],
                'INFO': fields[7]
            }
            indels.append(indel)
    return indels

def validate_vcf(vcf_file, reference_file, transcript_file, expected_indels):
    with open(reference_file, 'r') as file:
        reference_sequence = ''.join(file.read().split('\n')[1:])
    
    with open(transcript_file, 'r') as file:
        transcript_sequence = ''.join(file.read().split('\n')[1:])
    
    vcf_indels = parse_vcf(vcf_file)
    
    for indel in vcf_indels:
        chrom = indel['CHROM']
        pos = indel['POS']
        ref = indel['REF']
        alt = indel['ALT']
        
        if chrom != 'NC_000001.11':
            print(f"Error: Unexpected chromosome {chrom}")
            continue
        
        info_fields = indel['INFO'].split(';')
        transcript_pos = int([x.split('=')[1] for x in info_fields if 'TRANSCRIPT_POS' in x][0])
        transcript_id = [x.split('=')[1] for x in info_fields if 'TRANSCRIPT' in x][0]
        cigar_string = [x.split('=')[1] for x in info_fields if 'CIGAR' in x][0]
        
        for expected in expected_indels:
            if expected['POS'] == pos and expected['REF'] == ref and expected['ALT'] == alt:
                break
        else:
            print(f"Error: Indel at {pos} not found in expected indels")
            continue
        
        if ref != reference_sequence[pos-1:pos-1+len(ref)]:
            print(f"Error: Reference sequence mismatch at {pos}")
        
        if alt != reference_sequence[pos-1:pos-1+len(ref)] and not alt.startswith(reference_sequence[pos-1:pos-1+len(ref)]):
            print(f"Error: Alternate sequence mismatch at {pos}")
        
        expected_transcript_seq = transcript_sequence[transcript_pos-1:transcript_pos-1+len(ref)]
        if transcript_sequence[transcript_pos-1:transcript_pos-1+len(ref)] != ref:
            print(f"Error: Transcript sequence mismatch at transcript position {transcript_pos}")

def main():
    vcf_file = "outputs/test_indels.vcf"
    reference_file = "reference.fasta"
    transcript_file = "transcript.fasta"
    
    # Expected indels based on the test.sam and reference.fasta
    expected_indels = [
        {'POS': 8, 'REF': 'T', 'ALT': 'TA'},
        {'POS': 10, 'REF': 'A', 'ALT': 'A'},
        {'POS': 12, 'REF': 'A', 'ALT': 'AA'}
    ]
    
    if not os.path.isfile(vcf_file):
        print(f"Error: The VCF file {vcf_file} does not exist.")
        return
    
    if not os.path.isfile(reference_file):
        print(f"Error: The reference file {reference_file} does not exist.")
        return
    
    if not os.path.isfile(transcript_file):
        print(f"Error: The transcript file {transcript_file} does not exist.")
        return
    
    validate_vcf(vcf_file, reference_file, transcript_file, expected_indels)

if __name__ == "__main__":
    main()
