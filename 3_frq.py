import sys
import requests
import sqlite3
import argparse
import os
import time
from Bio import Entrez

# Set your email here for NCBI Entrez
Entrez.email = "david.mandia@yahoo.com"

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
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'qual': qual,
                'filter': filter_,
                'info': info
            })
    return variants

def process_variant(variant, server, headers, assembly, sleep_time=1):
    variant_str = f"{variant['chrom']} {variant['pos']} . {variant['ref']} {variant['alt']} . . ."
    data = {
        "variants": [variant_str],
        "variant_class": True,
        "allele_number": True,
        "regulatory": True,
        "canonical": True,
        "check_existing": True,
        "assembly": assembly
    } 
    response = requests.post(f"{server}/vep/human/variant", headers=headers, json=data)
    time.sleep(sleep_time)  # Sleep after each API call
    return response.json()

def parse_vep_data(variant_data):
    frequencies = {}
    genes = set()
    consequences = set()
    existing_variant = "Na"
    
    if 'colocated_variants' in variant_data:
        for cv in variant_data['colocated_variants']:
            if 'id' in cv:
                existing_variant = cv['id']
            if 'frequencies' in cv:
                gnomad_fre = next(iter(cv["frequencies"].values()), {})
                for pop, freq in gnomad_fre.items():
                    frequencies[pop] = freq

    if 'transcript_consequences' in variant_data:
        for tc in variant_data['transcript_consequences']:
            if 'gene_symbol' in tc:
                genes.add(tc['gene_symbol'])
            elif 'gene_id' in tc:
                genes.add(tc['gene_id'])
            if 'consequence_terms' in tc:
                consequences.update(tc['consequence_terms'])
 
    return frequencies, list(genes), list(consequences), existing_variant

def parse_info(info_str):
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
    return info_dict

def create_database(db_name):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS variants (
            chrom TEXT,
            pos INTEGER,
            ref TEXT,
            alt TEXT,
            qual REAL,
            filter TEXT,
            genomic_ref TEXT,
            operation TEXT,
            transcript_ref TEXT,
            transcript_pos TEXT,
            gnomadg REAL,
            gnomadg_eas REAL,
            gnomadg_nfe REAL,
            gnomadg_fin REAL,
            gnomadg_amr REAL,
            gnomadg_afr REAL,
            gnomadg_asj REAL,
            gnomadg_oth REAL,
            gnomadg_sas REAL,
            gnomadg_mid REAL,
            gnomadg_ami REAL,
            genes TEXT,
            consequences TEXT,
            existing_variant TEXT,
            rsid TEXT
        )
    ''')
    conn.commit()
    conn.close()

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

def process_results_and_insert(variant_data, original_variant, db_name):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    
    frequencies, genes, consequences, existing_variant = parse_vep_data(variant_data[0])
    
    info_dict = parse_info(original_variant['info'])
    
    variant_type = info_dict.get('TYPE', '')
    variant_len = info_dict.get('LEN', '')
    operation = f"{variant_type}:{variant_len}" if variant_type and variant_len else None
    
    transcript_ref = info_dict.get('TRANSCRIPT', '')
    transcript_pos = info_dict.get('TRANSCRIPT_POS', '')
    
    genomic_ref = info_dict.get('GENOME_REF', '')  
    
    # Check if frequencies are available, if not fetch using rsID
    if not frequencies and existing_variant != "Na":
        ncbi_response = fetch_frequencies(existing_variant)
        if ncbi_response:
            frequencies = parse_frequencies(ncbi_response)
    
    # Handle both GRCh37 and GRCh38 frequency keys
    if "af" in frequencies.keys():
        gnomad_fields = ['af', 'eas', 'eur', 'gnomadg_fin', 'amr', 'afr', 'gnomadg_asj', 'gnomadg_oth', 'sas', 'gnomadg_mid', 'gnomadg_ami']
    elif "gnomadg" in frequencies.keys():
        gnomad_fields = ['gnomadg', 'gnomadg_eas', 'gnomadg_nfe', 'gnomadg_fin', 'gnomadg_amr', 'gnomadg_afr', 'gnomadg_asj', 'gnomadg_oth', 'gnomadg_sas', 'gnomadg_mid', 'gnomadg_ami']
    else:
        gnomad_fields = ['gnomadg', 'gnomadg_eas', 'gnomadg_nfe', 'gnomadg_fin', 'gnomadg_amr', 'gnomadg_afr', 'gnomadg_asj', 'gnomadg_oth', 'gnomadg_sas', 'gnomadg_mid', 'gnomadg_ami']
    
    gnomad_values = [frequencies.get(field, "Na") for field in gnomad_fields]
    
    cursor.execute('''
        INSERT INTO variants (
            chrom, pos, ref, alt, qual, filter, genomic_ref,
            operation, transcript_ref, transcript_pos,
            gnomadg, gnomadg_eas, gnomadg_nfe, gnomadg_fin,
            gnomadg_amr, gnomadg_afr, gnomadg_asj, gnomadg_oth,
            gnomadg_sas, gnomadg_mid, gnomadg_ami, genes, consequences,
            existing_variant, rsid
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', (
        original_variant['chrom'],
        int(original_variant['pos']),
        original_variant['ref'],
        original_variant['alt'],
        float(original_variant['qual']) if original_variant['qual'] != '.' else None,
        original_variant['filter'],
        genomic_ref,
        operation,
        transcript_ref,
        transcript_pos,
        *gnomad_values,
        ','.join(genes),
        ','.join(consequences),
        existing_variant,
        existing_variant  # rsID column
    ))
    
    conn.commit()
    conn.close()

def main():
    parser = argparse.ArgumentParser(description="Process a VCF file, fetch VEP data, and store in a SQLite database.")
    parser.add_argument('vcf', help="Path to the input VCF file")
    parser.add_argument('--assembly', choices=['GRCh37', 'GRCh38'], required=True, help="Genome assembly (GRCh37 or GRCh38)")
    parser.add_argument('--sleep', type=float, default=1, help="Sleep time between API calls in seconds")
    args = parser.parse_args()

    vcf_file = args.vcf
    assembly = args.assembly
    sleep_time = args.sleep
    vcf_name = os.path.basename(vcf_file).split('.')[0]
    
    output_dir = "output"
    dbs_dir = os.path.join(output_dir, "database")
    os.makedirs(dbs_dir, exist_ok=True)
    
    db_name = os.path.join(dbs_dir, f'{vcf_name}_variant.db')
    
    # Remove existing database if it exists
    if os.path.exists(db_name):
        os.remove(db_name)
    
    create_database(db_name)
    
    vcf_data = read_vcf_file(vcf_file)
    
    ## The API is different based on the genome assembly used 
    if assembly == 'GRCh38':
        server = "https://rest.ensembl.org"
    elif assembly == 'GRCh37':
        server = "https://grch37.rest.ensembl.org"
        
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    # Process each variant individually
    for variant in vcf_data:
        result = process_variant(variant, server, headers, assembly, sleep_time)
        process_results_and_insert(result, variant, db_name)
    
    print(f'Data from {vcf_file} has been processed and inserted into {db_name}')

if __name__ == "__main__":
    main()
