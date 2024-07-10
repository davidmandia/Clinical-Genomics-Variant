import sys
import requests
import csv
from pysam import VariantFile

## Test with the batch results. See if you obtain something different from the previous results.


def get_vep_data(chrom, pos, ref, alt, strand=1):
    """Fetch data from Ensembl VEP REST API using VCF format"""
    server = "https://rest.ensembl.org"
    ext = "/vep/human/region"

    # Remove 'chr' prefix if present
    chrom = chrom.replace('chr', '')
    
    # Construct the variant in VCF-like format
    variant = f"{chrom} {pos} . {ref} {alt} . . ."
    
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }

    data = {
        "variants": [variant],
        "variant_class": True,
        "allele_number": True,
        "regulatory": True,
        "canonical": True,
        "hgvs": True
    }

    response = requests.post(server + ext, headers=headers, json=data)
   # print(response.json())
    
    if response.status_code != 200:
        print(f"Failed to retrieve data for {variant}: {response.status_code}")
        return None
    
    return response.json()

def parse_vep_data(data):
    """Parse VEP API response for allele frequencies and SIFT score"""
    if not data or len(data) == 0:
        return None, None, None
    print("data",data[0])
    variant_data = data[0]
    #print("keys",variant_data.keys())
    #print("variants id",variant_data["colocated_variants"][0]["id"])
    #print("variants keys",variant_data["colocated_variants"][0].keys())
    #next(iter(outer_object.values()))
    if "colocated_variants" in variant_data.keys():
       # print("colocated_variants",variant_data["colocated_variants"][0].keys())
        
        if "frequencies" in variant_data["colocated_variants"][0].keys():
            gnomad_fre = next(iter(variant_data["colocated_variants"][0]["frequencies"].values()))
          #  print("gnomad", type(gnomad_fre), gnomad_fre)
            for pop, freq in gnomad_fre.items():
                continue
             #   print(f"{pop}: {freq}")
    
    
    
    # Extract SIFT score
    sift_score = None
    if 'transcript_consequences' in variant_data:
        for tc in variant_data['transcript_consequences']:
            if 'sift_score' in tc:
                sift_score = tc['sift_score']
                break
    
    # Extract HGVS notation
    hgvs = variant_data.get('id', '')
    frequencies = {}
    return frequencies, sift_score, hgvs

def process_vcf(input_file, output_file):
    """Process VCF file, add frequencies and SIFT scores, and write to CSV file"""
    vcf_in = VariantFile(input_file)
    
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'HGVS', 'AF', 'SIFT']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        
        for record in vcf_in:
            chrom = record.chrom
            pos = record.pos
            ref = record.ref
            alt = record.alts[0]  # Assuming single alternative allele
            
            print(f"Processing variant: {chrom}:{pos}:{ref}>{alt}")
            
            vep_data = get_vep_data(chrom, pos, ref, alt)
            frequencies, sift_score, hgvs = parse_vep_data(vep_data)
            
            # Prepare row for CSV
            row = {
                'CHROM': chrom,
                'POS': pos,
                'ID': record.id,
                'REF': ref,
                'ALT': alt,
                'QUAL': record.qual,
                'FILTER': ','.join(record.filter.keys()),
                'INFO': ';'.join([f"{k}={v}" for k, v in record.info.items()]),
                'HGVS': hgvs,
                'AF': ','.join([f"{pop}:{freq}" for pop, freq in frequencies.items()]) if frequencies else '',
                'SIFT': sift_score if sift_score is not None else ''
            }
            
            writer.writerow(row)
    
    vcf_in.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_vcf> <output_csv>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_vcf(input_file, output_file)