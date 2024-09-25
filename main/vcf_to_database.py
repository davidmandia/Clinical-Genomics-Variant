#!/usr/bin/env python
import sys
import requests
import sqlite3
import argparse
import os
import multiprocessing
from functools import partial
import time


def choose_best_rsid(colocated_variants):
    """
    Chooses the best rsID based on the available frequency data and annotations.

    Args:
        colocated_variants (list): List of colocated variants from the VEP API response.

    Returns:
        str: The chosen rsID.
    """
    best_rsid = None
    best_frequency_count = 0

    for variant in colocated_variants:
        if 'frequencies' in variant:
            frequency_count = len(variant['frequencies'])
            if frequency_count > best_frequency_count:
                best_frequency_count = frequency_count
                best_rsid = variant['id']

    return best_rsid

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
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue  # Skip malformed lines
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
    print(f"Read {len(variants)} variants from the VCF file.")
    return variants

def process_batch(batch, server, headers, assembly, sleep_time=1):
    """
    Processes a batch of variants by querying the Ensembl VEP API.

    Args:
        batch (list): A list of variant dictionaries.
        server (str): The Ensembl VEP server URL.
        headers (dict): HTTP headers for the API request.
        assembly (str): The genome assembly (e.g., 'GRCh37', 'GRCh38').
        sleep_time (int): Time to sleep between API requests.

    Returns:
        list: The JSON response from the VEP API.
    """
    variants = [f"{v['chrom']} {v['pos']} . {v['ref']} {v['alt']} . . ." for v in batch]
    data = {
        "variants": variants,
        "variant_class": True,
        "allele_number": True,
        "regulatory": True,
        "canonical": True,
        "check_existing": True,
        "assembly": assembly
    }
    response = requests.post(f"{server}/vep/human/region", headers=headers, json=data)
    time.sleep(sleep_time)
    return response.json()

def parse_vep_data(variant_data):
    """
    Parses the response from the VEP API to extract relevant data.

    Args:
        variant_data (dict): The VEP API response for a variant.

    Returns:
        tuple: Frequencies, gene list, consequences list, and existing variant ID.
    """
    frequencies = {}
    genes = set()
    consequences = set()
    existing_variant = "Na"
    
    if 'colocated_variants' in variant_data:
        for cv in variant_data['colocated_variants']:             
                
            if 'frequencies' in cv:
                #rsid is assigned to the best rsid based on the available frequency data and annotations
                gnomad_fre = next(iter(cv["frequencies"].values()), {})
                for pop, freq in gnomad_fre.items():
                    frequencies[pop] = freq
                existing_variant = cv['id']
            elif "frequencies" not in cv: # if best (frequency) not available, the rsid is assigned to the first rsid in the list
                existing_variant = cv['id']

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
    """
    Parses the INFO field from the VCF file.

    Args:
        info_str (str): The INFO field as a semicolon-separated string.

    Returns:
        dict: A dictionary with key-value pairs extracted from the INFO field.
    """
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
    return info_dict

def create_database(db_name):
    """
    Creates a SQLite database with a specific schema for storing variant data.

    Args:
        db_name (str): The name of the SQLite database file.
    """

    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    ## remove qual and filter from the table, add latest version
    cursor.execute(f'''
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
            af REAL,
            af_eas REAL,
            af_nfe REAL,
            af_fin REAL,
            af_amr REAL,
            af_afr REAL,
            af_asj REAL,
            af_oth REAL,
            af_sas REAL,
            af_mid REAL,
            af_ami REAL,
            genes TEXT,
            consequences TEXT,
            existing_variant TEXT
        )
    ''')
    conn.commit()
    conn.close()

def process_results_and_insert(results, original_variants, db_name):
    """
    Processes VEP API results and inserts them into the SQLite database.

    Args:
        results (list): The VEP API results.
        original_variants (list): The original variants from the VCF file.
        db_name (str): The SQLite database file name.
    """
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    
    insert_data = []
    for variant_data, original_variant in zip(results, original_variants):
        frequencies, genes, consequences, existing_variant = parse_vep_data(variant_data)
        
        info_dict = parse_info(original_variant['info'])
        
        variant_type = info_dict.get('TYPE', '')
        variant_len = info_dict.get('LEN', '')
        operation = f"{variant_type}:{variant_len}" if variant_type and variant_len else None
        
        transcript_ref = info_dict.get('TRANSCRIPT', '')
        transcript_pos = info_dict.get('TRANSCRIPT_POS', '')
        
        genomic_ref = info_dict.get('GENOME_REF', '')  
        
        # Determine which gnomAD fields to use based on the presence of data
        #for 1000 genomes project results in GRCh37
        if "af" in frequencies.keys():
            gnomad_fields = ['af', 'eas', 'eur', 'af_fin', 'amr', 'afr', 'af_asj', 'af_oth', 'sas', 'af_mid', 'af_ami']
        
        elif "gnomadg" in frequencies.keys():
            gnomad_fields = ['gnomadg', 'gnomadg_eas', 'gnomadg_nfe', 'gnomadg_fin', 'gnomadg_amr', 'gnomadg_afr', 'gnomadg_asj', 'gnomadg_oth', 'gnomadg_sas', 'gnomadg_mid', 'gnomadg_ami']
        else:
            gnomad_fields = ['gnomadg', 'gnomadg_eas', 'gnomadg_nfe', 'gnomadg_fin', 'gnomadg_amr', 'gnomadg_afr', 'gnomadg_asj', 'gnomadg_oth', 'gnomadg_sas', 'gnomadg_mid', 'gnomadg_ami']   

        gnomad_values = [frequencies.get(field, None ) for field in gnomad_fields]
        
        insert_data.append((
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
            existing_variant
        ))
    
    cursor.executemany('''
        INSERT INTO variants (
            chrom, pos, ref, alt, qual, filter, genomic_ref,
            operation, transcript_ref, transcript_pos,
            af, af_eas, af_nfe, af_fin,
            af_amr, af_afr, af_asj, af_oth,
            af_sas, af_mid, af_ami, genes, consequences,
            existing_variant
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', insert_data)
    
    conn.commit()
    conn.close()
    print(f"Inserted {len(insert_data)} variants into the database.")

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
    
    if os.path.exists(db_name):
        os.remove(db_name)
    
    create_database(db_name)
    
    vcf_data = read_vcf_file(vcf_file)
    
    if assembly == 'GRCh38':
        server = "https://rest.ensembl.org"
    elif assembly == 'GRCh37':
        server = "https://grch37.rest.ensembl.org"
        
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    # Create batches of 200 variants
    batches = [vcf_data[i:i+200] for i in range(0, len(vcf_data), 200)]
    print(f"Created {len(batches)} batches for processing.")
    
    # Set up multiprocessing
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    process_func = partial(process_batch, server=server, headers=headers, assembly=assembly, sleep_time=sleep_time)
    
    # Process batches in parallel
    results = pool.map(process_func, batches)
    
    # Flatten results
    all_results = [item for sublist in results for item in sublist]
    print(f"Processed {len(all_results)} variants from VEP API.")
    
    # Process results and insert into database
    process_results_and_insert(all_results, vcf_data, db_name)
    
    print(f'Data from {vcf_file} has been processed and inserted into {db_name}')

if __name__ == "__main__":
    main()
