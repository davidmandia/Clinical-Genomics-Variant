import requests
import argparse
import json

#https://github.com/aws-samples/aws-genomics-datalake/blob/main/1000Genomes.ipynb
#https://gnomad.broadinstitute.org/variant/1-55051215-G-GA?dataset=gnomad_r4
#https://gnomad.broadinstitute.org/downloads#v4-resources
#https://registry.opendata.aws/broad-gnomad/ 
#https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/
#https://broadinstitute.github.io/gnomad_methods/api_reference/utils/vep.html#gnomad.utils.vep.CURRENT_VEP_VERSION
#https://rest.ensembl.org/#VEP
#https://rest.ensembl.org/documentation/info/vep_hgvs_get

# Function to fetch gnomAD data for a variant
def fetch_gnomad_data(chrom, pos, ref, alt):
    url = f"https://gnomad.broadinstitute.org/api/v3/variant/{chrom}-{pos}-{ref}-{alt}"
    headers = {
        "Content-Type": "application/json"
    }

    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching gnomAD data for variant {chrom}-{pos}-{ref}-{alt}: {e}")
        return None

# Function to parse the VCF file (example function, replace with your actual parsing logic)
def parse_vcf(file_path):
    # Example placeholder function, replace with actual parsing logic
    variants = [
        {"CHROM": "1", "POS": 136889, "REF": "GG", "ALT": "G"},
        {"CHROM": "1", "POS": 1052530, "REF": "C", "ALT": "CTG"}
        # Add more variants as needed
    ]
    return variants

# Main function to process VCF file and fetch gnomAD data
def main():
    parser = argparse.ArgumentParser(description="Process VCF file and fetch gnomAD data.")
    parser.add_argument('vcf_file', type=str, help="Path to the VCF file")

    args = parser.parse_args()

    # Parse the VCF file
    variants = parse_vcf(args.vcf_file)

    # Fetch gnomAD data for each variant and print results
    for variant in variants:
        chrom = variant['CHROM']
        pos = variant['POS']
        ref = variant['REF']
        alt = variant['ALT']

        print(f"Fetching gnomAD data for variant {chrom}-{pos}-{ref}-{alt}...")
        gnomad_data = fetch_gnomad_data(chrom, pos, ref, alt)

        if gnomad_data:
            print(json.dumps(gnomad_data, indent=2))  # Print gnomAD data

        print()  # Add a newline for clarity between variants

if __name__ == "__main__":
    main()
