import sqlite3
import pandas as pd
import argparse

def read_vcf_file(vcf_file):
    variants = []
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos, _, ref, alt, qual, filter_, info = fields[:8]
            variants.append({
                'chrom': chrom,
                'pos': int(pos),
                'ref': ref,
                'alt': alt,
                'qual': qual,
                'filter': filter_,
                'info': info
            })
    return variants

def check_expected_variations(db_path, vcf_variants):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    expected_variations = []
    unexpected_variations = []

    for variant in vcf_variants:
        cursor.execute("""
            SELECT * FROM variants
            WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ?
        """, (variant['chrom'], variant['pos'], variant['ref'], variant['alt']))

        result = cursor.fetchone()
        if result:
            expected_variations.append({
                'variant': variant,
                'clinical_label': result[23],  # Adjust the index based on your schema
                'clinically_relevant': result[22],  # Adjust the index based on your schema
                'extra_data': result  # Add any extra relevant data here
            })
        else:
            unexpected_variations.append(variant)

    conn.close()
    return expected_variations, unexpected_variations

def output_results(expected_variations, unexpected_variations):
    print("Expected Variations:")
    for var in expected_variations:
        print(f"Variant: {var['variant']}")
        print(f"Clinical Label: {var['clinical_label']}")
        print(f"Clinically Relevant: {var['clinically_relevant']}")
        print(f"Extra Data: {var['extra_data']}")
        print("")

    print("Unexpected Variations:")
    for var in unexpected_variations:
        print(f"Variant: {var}")
        print("")

def main():
    parser = argparse.ArgumentParser(description="Check for expected variations in a given VCF input and interface with the database.")
    parser.add_argument('db_path', help="Path to the SQLite database")
    parser.add_argument('vcf_file', help="Path to the input VCF file")
    args = parser.parse_args()

    db_path = args.db_path
    vcf_file = args.vcf_file

    vcf_variants = read_vcf_file(vcf_file)
    expected_variations, unexpected_variations = check_expected_variations(db_path, vcf_variants)

    output_results(expected_variations, unexpected_variations)

if __name__ == "__main__":
    main()
