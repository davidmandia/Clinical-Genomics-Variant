import csv
import re

def parse_sam(file_path):
    alignments = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('@'):
                # Skip header lines
                continue
            fields = line.strip().split('\t')
            read_id = fields[0]
            cigar = fields[5]
            alignments.append((read_id, cigar))
    return alignments

def analyze_cigar(cigar):
    operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    analysis = {}
    for length, op in operations:
        if op not in analysis:
            analysis[op] = 0
        analysis[op] += int(length)
    return analysis

# Path to your SAM file
sam_file_path = 'output.sam'

# Parse the SAM file
alignments = parse_sam(sam_file_path)

# Analyze the CIGAR strings
output_rows = []
for read_id, cigar in alignments:
    analysis = analyze_cigar(cigar)
    output_rows.append([read_id, cigar] + [analysis.get(op, 0) for op in '=ND'])

# Write the output to a CSV file
output_file_path = 'output.csv'
with open(output_file_path, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['Read ID', 'CIGAR', '=', 'N', 'D'])  # Write header
    csvwriter.writerows(output_rows)

print(f"Output written to {output_file_path}")
