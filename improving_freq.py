import requests
import json
import argparse
import time

def annotate_with_vep(variant_lines, server="https://rest.ensembl.org", assembly="GRCh38"):
    annotated_variants = []
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    batch_size = 200

    for i in range(0, len(variant_lines), batch_size):
        batch = variant_lines[i:i+batch_size]
        data = {
            "variants": batch,
            "vcf_string": 1,
            "canonical": 1,
            "picks": 1,
            "variant_class": 1,
            "frequency": 1,
            "population_freqs": 1,
            "existing_variation": 1,
            "dbNSFP": 1,
            "all_consequences": 1,
            "assembly": assembly
        }
        
        response = requests.post(f"{server}/vep/human/region", headers=headers, json=data)
        
        if response.status_code != 200:
            print(f"Error: API request failed with status code {response.status_code}")
            print(response.text)
            return None
        
        annotated_variants.extend(response.json())
        print("annotated_variants", annotated_variants)
        time.sleep(3)  # To respect the rate limit

    return annotated_variants

def read_vcf_file(vcf_file):
    with open(vcf_file, 'r') as file:
        lines = file.readlines()
    return lines

def write_annotated_vcf(original_vcf_lines, annotated_variants, output_vcf):
    header = [line for line in original_vcf_lines if line.startswith('#')]
    variant_lines = [line for line in original_vcf_lines if not line.startswith('#')]

    with open(output_vcf, 'w') as out_file:
        out_file.writelines(header)
        for original_line, variant_data in zip(variant_lines, annotated_variants):
            fields = original_line.strip().split('\t')
            rs_id = next((v['id'] for v in variant_data.get('colocated_variants', []) if v.get('id')), '.')
            # If VEP returns a different coordinate, update the VCF line accordingly
            if 'start' in variant_data and int(fields[1]) != variant_data['start']:
                fields[1] = str(variant_data['start'])
            if 'allele_string' in variant_data:
                ref, alt = variant_data['allele_string'].split('/')
                fields[3] = ref
                fields[4] = alt
            fields[2] = rs_id  # Update the ID column with rsID
            out_file.write('\t'.join(fields) + '\n')

def main():
    parser = argparse.ArgumentParser(description="Annotate and validate variants using Ensembl VEP and create an annotated VCF file.")
    parser.add_argument('input_vcf', help="Path to the input VCF file")
    parser.add_argument('output_vcf', help="Path to the output annotated VCF file")
    args = parser.parse_args()

    original_vcf_lines = read_vcf_file(args.input_vcf)
    variant_lines = [line.strip().split('\t', 5)[0:5] for line in original_vcf_lines if not line.startswith('#')]
    variant_lines = [f"{fields[0]} {fields[1]} . {fields[3]} {fields[4]}" for fields in variant_lines]

    annotated_variants = annotate_with_vep(variant_lines)
    
    if annotated_variants:
        write_annotated_vcf(original_vcf_lines, annotated_variants, args.output_vcf)
        print(f"Annotated VCF has been written to {args.output_vcf}")

if __name__ == "__main__":
    main()
