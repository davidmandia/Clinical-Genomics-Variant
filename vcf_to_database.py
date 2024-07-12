import sys
import requests
import sqlite3
import argparse
import os
import multiprocessing
from functools import partial
import time

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

def process_batch(batch, server, headers, assembly, sleep_time=1):
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
    print("info_dict", info_dict)
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
            existing_variant TEXT
        )
    ''')
    conn.commit()
    conn.close()

def process_results_and_insert(results, original_variants, db_name):
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
        
        
        # in CHG38 the eur gnomadg field is called gnomadg_nfe vs 37 is eur
        # in CHG38 the average gnomadg vs 37 is af
        # the rest in Ch38 they have a prefix gnomadg_ vs 37 they don't
        # frequencies = {'af': "Na", 'eas': "Na", 'amr': "Na", 'sas': "Na", 'afr': "Na", 'eur': "Na"}
        
        
        
        gnomad_fields = [
            'gnomadg', 'gnomadg_eas', 'gnomadg_nfe', 'gnomadg_fin',
            'gnomadg_amr', 'gnomadg_afr', 'gnomadg_asj', 'gnomadg_oth',
            'gnomadg_sas', 'gnomadg_mid', 'gnomadg_ami'
        ]
        print("frequencies", frequencies)
        ## Can probabaly add a conditional statement to check if the key is in the dictionary
        ## if not add it with a value of "Na"
        
        gnomad_values = [frequencies.get(field, "Na") for field in gnomad_fields]
        
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
            gnomadg, gnomadg_eas, gnomadg_nfe, gnomadg_fin,
            gnomadg_amr, gnomadg_afr, gnomadg_asj, gnomadg_oth,
            gnomadg_sas, gnomadg_mid, gnomadg_ami, genes, consequences,
            existing_variant
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', insert_data)
    
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
    
    db_name = os.path.join(dbs_dir, f'{assembly}_indels_variant.db')
    
    # Remove existing database if it exists
    if os.path.exists(db_name):
        os.remove(db_name)
    
    create_database(db_name)
    
    vcf_data = read_vcf_file(vcf_file)
    
    ## The APi is different based on the genome assembly used 
    if assembly == 'GRCh38':
        server = "https://rest.ensembl.org"
    elif assembly == 'GRCh37':
        server = "https://grch37.rest.ensembl.org"
        
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    # Create batches of 200 variants
    batches = [vcf_data[i:i+200] for i in range(0, len(vcf_data), 200)]
    
    # Set up multiprocessing
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    process_func = partial(process_batch, server=server, headers=headers, assembly=assembly, sleep_time=sleep_time)
    
    # Process batches in parallel
    results = pool.map(process_func, batches)
    
    # Flatten results
    all_results = [item for sublist in results for item in sublist]
    
    # Process results and insert into database
    process_results_and_insert(all_results, vcf_data, db_name)
    
    print(f'Data from {vcf_file} has been processed and inserted into {db_name}')

if __name__ == "__main__":
    main()