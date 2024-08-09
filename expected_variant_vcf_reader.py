import os
import sqlite3
import argparse

def extract_chromosome_digit(chrom):
    """
    Extracts the numeric part of the chromosome identifier.

    Args:
        chrom (str): Chromosome identifier (e.g., 'chr1', 'chrX').

    Returns:
        str: The numeric part of the chromosome as a string, or 'X', 'Y', 'MT' for non-numeric chromosomes.
    """
    if chrom.startswith('chr'):
        chrom = chrom[3:]  # Remove 'chr' prefix if present
    
    return chrom  # Return the remaining part (e.g., '1', 'X', 'Y', 'MT')

def read_vcf_file(vcf_file):
    """
    Reads a VCF file and extracts variant information.

    Args:
        vcf_file (str): Path to the VCF file.

    Returns:
        list: A list of dictionaries, each containing variant information.
    """
    variants = []
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            fields = line.strip().split()
            if len(fields) < 8:
                print(f"Skipping malformed line: {line.strip()}")
                continue
            chrom, pos, _, ref, alt, qual, filter_, info = fields[:8]
            chrom_digit = extract_chromosome_digit(chrom)
            variants.append({
                'chrom': chrom_digit,
                'pos': int(pos),
                'ref': ref,
                'alt': alt,
                'qual': qual,
                'filter': filter_,
                'info': info
            })
    print("Variants extracted from VCF")
    return variants

def check_expected_variations(db_path, vcf_variants):
    """
    Checks for expected variations in the database against a list of VCF variants.

    Args:
        db_path (str): Path to the SQLite database.
        vcf_variants (list): List of variants extracted from a VCF file.

    Returns:
        tuple: Two lists, one for expected and one for unexpected variations.
    """
    print("Connecting to the database...")
    if not os.path.exists(db_path):
        print(f"Database file {db_path} does not exist.")
        return [], []

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
            gnomad_data = {
                'af': result[10]  # Adjust this index as needed based on your table structure
            }
            gene_affected = result[24]
            clinical_label = result[27]
            clinically_relevant = result[26]

            expected_variations.append({
                'variant': variant,
                'gnomad_data': gnomad_data,
                'affected_gene': gene_affected,
                'clinical_label': clinical_label,
                'clinically_relevant': clinically_relevant
            })
        else:
            unexpected_variations.append(variant)

    conn.close()
    print("Connection to database closed.")
    return expected_variations, unexpected_variations

def output_results(expected_variations, unexpected_variations):
    """
    Outputs the results of the expected and unexpected variations, and prints the count of expected variations.

    Args:
        expected_variations (list): List of expected variant information.
        unexpected_variations (list): List of unexpected variant information.
    """
    expected_count = len(expected_variations)
    #print(f"Number of Expected Variations: {expected_count}")

    print("\nExpected Variations:")
    for var in expected_variations:
        print(f"Variant: {var['variant']}")
        print("Affected Gene:", var['affected_gene'])
        print(f"Clinical Label: {var['clinical_label']}")
        print(f"Clinically Relevant: {var['clinically_relevant']}")
        print(f"gnomAD Data: {var['gnomad_data']}")
        print("")

    print(f"Number of Expected Variations: {expected_count}")
    print("Unexpected Variations:")
    #### I have commented the unexpected variants, the ones that are not linked to a gap in the genome)
    # for var in unexpected_variations:
    #     print(f"Variant: {var}")
    #     print("")

def main():
    """
    Main function to handle argument parsing and function execution.
    """
    parser = argparse.ArgumentParser(description="Check for expected variations in a given VCF input and interface with the database.")
    parser.add_argument('db_path', help="Path to the SQLite database")
    parser.add_argument('vcf_file', help="Path to the input VCF file")
    args = parser.parse_args()

    db_path = os.path.abspath(args.db_path)
    vcf_file = os.path.abspath(args.vcf_file)

    vcf_variants = read_vcf_file(vcf_file)
    expected_variations, unexpected_variations = check_expected_variations(db_path, vcf_variants)

    output_results(expected_variations, unexpected_variations)

if __name__ == "__main__":
    main()
