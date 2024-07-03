import subprocess
import os
import pytest

# Define pytest fixtures for your files and db
@pytest.fixture
def vcf_file():
    return "output/VCFs/sample_False_indels.vcf"

@pytest.fixture
def sam_file():
    return "input/sample.sam"

@pytest.fixture
def db():
    return "input/blastdb/GRCh38"

# Function to fetch sequence from BLAST DB
def get_sequence_blast_db(db, accession, start, end):
    cmd = ["blastdbcmd", "-db", db, "-entry", accession, "-range", f'{start}-{end}']
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    sequence = "".join(result.stdout.split("\n")[1:])
    return sequence

# Function to parse SAM file and retrieve transcriptomic sequence
def get_transcriptomic_sequence_from_sam(sam_file, read_id):
    with open(sam_file, 'r') as file:
        for line in file:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 11:
                raise ValueError("Invalid SAM format")
            if fields[0] == read_id:
                return fields[9]
    raise ValueError(f"Read ID '{read_id}' not found in SAM file")

# Test function using pytest fixtures
def test_vcf_records(vcf_file, sam_file, db):
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'):  # Skip header lines
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            ref_seq = fields[3]
            alt_seq = fields[4]

            # Extract information from INFO field
            info = fields[7]
            transcript_id = info.split()[1]
            transcript_pos = int(info.split()[3])
            genome_ref_accession = info.split()[-1]

            # Get transcriptomic sequence from SAM file
            transcriptomic_sequence = get_transcriptomic_sequence_from_sam(sam_file, transcript_id)[transcript_pos-1: transcript_pos-1 + len(alt_seq)]

            # Get genomic reference sequence from BLAST DB
            end_Nc = pos + len(ref_seq) - 1
            genomic_sequence = get_sequence_blast_db(db, genome_ref_accession, pos , end_Nc)

            # Validate sequences
            assert genomic_sequence == ref_seq, f"Genomic sequence mismatch for {chrom}:{pos}"
            assert transcriptomic_sequence == alt_seq, f"Transcriptomic sequence mismatch for transcript {transcript_id}"

# Entry point for running the test with pytest
if __name__ == "__main__":
    pytest.main()
