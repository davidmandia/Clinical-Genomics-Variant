import requests
import argparse
import json
import time

def read_vcf(file_path):
    variants = []
    with open(file_path, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos, _, ref, alt = fields[:5]
            variants.append(f"{chrom} {pos} {ref} {alt}")
    return variants

def process_variants(variants, server, assembly, sleep_time=1):
    data = {
        "variants": variants,
        "variant_class": True,
        "allele_number": True,
        "regulatory": True,
        "canonical": True,
        "check_existing": True,
        "assembly": assembly
    }
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    response = requests.post(f"{server}/vep/human/region", headers=headers, json=data)
    
    time.sleep(sleep_time)  # Sleep to avoid overwhelming the API
    
    if response.status_code != 200:
        print(f"Error: API request failed with status code {response.status_code}")
        print(f"Response content: {response.text}")
        return None
    
    return response.json()

def parse_frequencies(variant_data):
    frequencies = {}
    if 'colocated_variants' in variant_data:
        for cv in variant_data['colocated_variants']:
            if 'frequencies' in cv:
                for allele, freq_data in cv['frequencies'].items():
                    frequencies[allele] = freq_data
    return frequencies

def main():
    parser = argparse.ArgumentParser(description="Test VEP API for variants in a VCF file")
    parser.add_argument('vcf', help="Path to the input VCF file")
    parser.add_argument('--assembly', choices=['GRCh37', 'GRCh38'], required=True, help="Genome assembly")
    parser.add_argument('--limit', type=int, default=10, help="Limit the number of variants to process")
    args = parser.parse_args()

    vcf_file = args.vcf
    assembly = args.assembly
    limit = args.limit

    if assembly == 'GRCh38':
        server = "https://rest.ensembl.org"
    elif assembly == 'GRCh37':
        server = "https://grch37.rest.ensembl.org"

    variants = read_vcf(vcf_file)[:limit]  # Limit the number of variants
    
    batch_size = 200  # VEP API limit
    for i in range(0, len(variants), batch_size):
        batch = variants[i:i+batch_size]
        print(f"Processing batch {i//batch_size + 1} ({len(batch)} variants)")
        
        result = process_variants(batch, server, assembly)
        
        if result:
            for variant_data in result:
                print("\nVariant:", variant_data.get('input'))
                print("Most severe consequence:", variant_data.get('most_severe_consequence'))
                
                frequencies = parse_frequencies(variant_data)
                if frequencies:
                    print("Frequencies:")
                    print(json.dumps(frequencies, indent=2))
                else:
                    print("No frequency data available")
                
                print("-" * 50)
        else:
            print("No result returned from API for this batch")
        
        print("\n" + "=" * 50 + "\n")

if __name__ == "__main__":
    main()