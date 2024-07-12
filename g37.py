import argparse
import requests
import sqlite3
import time
import os

def read_vcf(file_path):
    variants = []
    with open(file_path, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            variants.append((chrom, pos, ref, alt))
    return variants

def process_variant(variant, server, headers, sleep_time=1):
    chrom, pos, ref, alt = variant
    data = {
        "variants": [f"{chrom} {pos} . {ref} {alt} . . ."],
        "species": "human",
        "vcf_string": 1,
        "canonical": 1,
        "picks": 1,
        "variant_class": 1
    }
    
    response = requests.post(f"{server}/vep/human/region", headers=headers, json=data)
    time.sleep(sleep_time)
    
    if response.status_code != 200:
        print(f"Error: API request failed for {chrom}:{pos} {ref}>{alt} with status code {response.status_code}")
        return None
    
    return response.json()

def parse_vep_data(variant_data):
    frequencies = {}
    genes = set()
    
    if variant_data:
        for data in variant_data:
            if 'colocated_variants' in data:
                for cv in data['colocated_variants']:
                    if 'frequencies' in cv:
                        for allele, freqs in cv['frequencies'].items():
                            for pop, freq in freqs.items():
                                if pop in ['af', 'eas', 'amr', 'sas', 'afr', 'eur']:
                                    frequencies[pop] = freq
            if 'transcript_consequences' in data:
                for consequence in data['transcript_consequences']:
                    genes.add(consequence.get('gene_symbol', 'N/A'))
    
    return {
        'frequencies': frequencies,
        'genes': list(genes)
    }

def create_database(db_name):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS variants (
            chrom TEXT,
            pos INTEGER,
            ref TEXT,
            alt TEXT,
            genes TEXT,
            af REAL,
            eas REAL,
            amr REAL,
            sas REAL,
            afr REAL,
            eur REAL
        )
    ''')
    conn.commit()
    conn.close()

def insert_variant(db_name, variant_data, parsed_data):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    
    chrom, pos, ref, alt = variant_data
    genes = ','.join(parsed_data['genes'])
    frequencies = parsed_data['frequencies']
    
    cursor.execute('''
        INSERT INTO variants (
            chrom, pos, ref, alt, genes, af, eas, amr, sas, afr, eur
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', (
        chrom, pos, ref, alt, genes, 
        frequencies.get('af'), frequencies.get('eas'), frequencies.get('amr'), 
        frequencies.get('sas'), frequencies.get('afr'), frequencies.get('eur')
    ))
    
    conn.commit()
    conn.close()

def main():
    parser = argparse.ArgumentParser(description="Process GRCh37 variants from a VCF file, fetch VEP data, and store in a SQLite database.")
    parser.add_argument('vcf', help="Path to the input VCF file")
    parser.add_argument('--db', help="Path to the SQLite database", default="variants.db")
    parser.add_argument('--sleep', type=float, default=1, help="Sleep time between API calls in seconds")
    args = parser.parse_args()

    vcf_file = args.vcf
    db_name = args.db
    sleep_time = args.sleep
    
    if os.path.exists(db_name):
        os.remove(db_name)
    
    create_database(db_name)
    
    variants = read_vcf(vcf_file)
    server = "https://grch37.rest.ensembl.org"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    for variant in variants:
        vep_data = process_variant(variant, server, headers, sleep_time)
        if vep_data:
            parsed_data = parse_vep_data(vep_data)
            insert_variant(db_name, variant, parsed_data)
            chrom, pos, ref, alt = variant
            print(f"Variant: {chrom} {pos} {ref}>{alt}")
            print(f"Genes: {', '.join(parsed_data['genes'])}")
            print(f"Frequencies: {parsed_data['frequencies']}")
            print("")

if __name__ == "__main__":
    main()
