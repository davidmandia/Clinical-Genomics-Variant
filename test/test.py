import sys
import os
import pytest
import subprocess
from unittest.mock import patch, mock_open

# Ensure the main script can be imported
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import functions from the main script
from main import (
    get_sequence_blast_db,
    identify_short_indels,
    extract_chromosome_number,
    parse_sam,
    write_to_vcf,
    write_missing_sequences
)

# Test the get_sequence_blast_db function
@patch('subprocess.run')
def test_get_sequence_blast_db(mock_subproc_run):
    mock_subproc_run.return_value = subprocess.CompletedProcess(args=['blastdbcmd'], returncode=0, stdout=">header\nATGCATGCATGC")
    
    db = "test_db"
    accession = "test_accession"
    sequence = get_sequence_blast_db(db, accession, 0, 11, sleep_time=0)
    
    assert sequence == "ATGCATGCATGC"

# Test the identify_short_indels function


def test_identify_short_indels():
    cigar_string = "5M2I3M2D4M"
    short_indels = identify_short_indels(cigar_string)
    expected_indels = [('I', 2, 5, 5), ('D', 2, 8, 10)]
    print(short_indels, expected_indels)
    assert short_indels == expected_indels


# Test the extract_chromosome_number function
def test_extract_chromosome_number():
    genome_ref = "NC_000001.11"
    chrom_number = extract_chromosome_number(genome_ref)
    
    assert chrom_number == "1"

# Mock the open function to test parse_sam without reading an actual file
@patch('builtins.open', new_callable=mock_open, read_data="@header\nread1\t0\tNC_000001.11\t1\t255\t5M2I3M2D4M\t*\t0\t0\tACGTACGTACGT\t*")
@patch('main.get_sequence_blast_db')
def test_parse_sam(mock_get_sequence_blast_db, mock_open):
    mock_get_sequence_blast_db.side_effect = [
        "A",  # Return value for insertion
        "AT"  # Return value for deletion
    ]
    
    db = "test_db"
    pseudo = False
    indels, missing_sequences = parse_sam("test.sam", db, pseudo, sleep_time=0)
    
    assert len(indels) == 2
    assert indels[0]['REF'] == "A"
    assert indels[1]['REF'] == "AT"
    assert len(missing_sequences) == 0

# Test the write_to_vcf function
def test_write_to_vcf(tmp_path):
    indels = [
        {'CHROM': '1', 'POS': 1, 'ID': '.', 'REF': 'A', 'ALT': 'AC', 'QUAL': 99, 'FILTER': 'PASS', 'INFO': 'INFO1'},
        {'CHROM': '1', 'POS': 2, 'ID': '.', 'REF': 'T', 'ALT': 'TA', 'QUAL': 99, 'FILTER': 'PASS', 'INFO': 'INFO2'}
    ]
    output_file = tmp_path / "output.vcf"
    write_to_vcf(indels, output_file)
    
    with open(output_file, 'r') as f:
        lines = f.readlines()
    
    assert lines[3].strip() == "1\t1\t.\tA\tAC\t99\tPASS\tINFO1"
    assert lines[4].strip() == "1\t2\t.\tT\tTA\t99\tPASS\tINFO2"

# Test the write_missing_sequences function
def test_write_missing_sequences(tmp_path):
    missing_sequences = [
        {'TRANSCRIPT': 'read1', 'GENOME_REF': 'NC_000001.11', 'TRANSCRIPT_POS': 1, 'GENOME_POS': 1, 'OP': 'D', 'LENGTH': 2, 'CIGAR': '5M2D5M'},
        {'TRANSCRIPT': 'read2', 'GENOME_REF': 'NC_000002.11', 'TRANSCRIPT_POS': 2, 'GENOME_POS': 2, 'OP': 'I', 'LENGTH': 3, 'CIGAR': '5M3I5M'}
    ]
    output_file = tmp_path / "missing_sequences.txt"
    write_missing_sequences(missing_sequences, output_file)
    
    with open(output_file, 'r') as f:
        lines = f.readlines()
    
    assert lines[1].strip() == "read1\tNC_000001.11\t1\t1\tD\t2\t5M2D5M"
    assert lines[2].strip() == "read2\tNC_000002.11\t2\t2\tI\t3\t5M3I5M"
