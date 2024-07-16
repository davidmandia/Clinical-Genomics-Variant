from Bio import Entrez
import argparse
import json

# Set your email here
Entrez.email = "david.mandia@yahoo.com"

def fetch_gene_name(rsid):
    try:
        handle = Entrez.elink(dbfrom="snp", id=rsid.replace('rs', ''), db="gene")
        record = Entrez.read(handle)
        handle.close()
        gene_ids = record[0]['LinkSetDb'][0]['Link']
        gene_names = []
        for gene_id in gene_ids:
            handle = Entrez.esummary(db="gene", id=gene_id['Id'])
            uid_record = Entrez.read(handle)
            handle.close()
            gene_names.append(uid_record['DocumentSummarySet']['DocumentSummary'][0]['Name'])
        return gene_names
    except Exception as e:
        print(f"Error fetching gene name for rsID {rsid}: {e}")
        return None

def fetch_frequencies(rsid):
    try:
        handle = Entrez.efetch(db="snp", id=rsid.replace('rs', ''), rettype="docsum", retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        return record
    except Exception as e:
        print(f"Error fetching frequencies for rsID {rsid}: {e}")
        return None

def parse_frequencies(record):
    frequencies = {}
    if record and 'DocumentSummarySet' in record:
        for docsum in record['DocumentSummarySet']['DocumentSummary']:
            if 'GLOBAL_MAFS' in docsum:
                for maf in docsum['GLOBAL_MAFS']:
                    if 'STUDY' in maf and 'FREQ' in maf:
                        study = maf['STUDY']
                        freq = maf['FREQ']
                        allele_freq = freq.split('=')[1].split('/')[0]
                        frequencies[study] = float(allele_freq)
    return frequencies

def read_annotated_vcf(annotated_vcf_file):
    with open(annotated_vcf_file, 'r') as file:
        lines = file.readlines()
    return [line for line in lines if not line.startswith('#')]

def main():
    parser = argparse.ArgumentParser(description="Retrieve gene names and frequencies using rsIDs from an annotated VCF file.")
    parser.add_argument('annotated_vcf', help="Path to the annotated VCF file")
    args = parser.parse_args()

    vcf_lines = read_annotated_vcf(args.annotated_vcf)
    rsids = [line.split('\t')[2].strip() for line in vcf_lines if line.split('\t')[2].strip() != '.']
    
    for rsid in rsids:
        print(f"Fetching gene names for rsID: {rsid}")
        gene_names = fetch_gene_name(rsid)
        if gene_names:
            print(f"Gene names: {gene_names}")
        else:
            print(f"Failed to retrieve gene names for rsID {rsid}")
        
        print(f"Fetching frequencies for rsID: {rsid}")
        ncbi_response = fetch_frequencies(rsid)
        if ncbi_response:
            print(json.dumps(ncbi_response, indent=2))  # Debug: print the structure of the response
            frequencies = parse_frequencies(ncbi_response)
            print(f"Frequencies: {frequencies}")
        else:
            print(f"Failed to retrieve data for rsID {rsid}")
        print("\n")

if __name__ == "__main__":
    main()
