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
            chrom, pos, rsid, ref, alt, qual, filter_, info = fields[:8]
            variants.append({
                'chrom': chrom,
                'pos': pos,
                'rsid': rsid,
                'ref': ref,
                'alt': alt,
                'qual': qual,
                'filter': filter_,
                'info': info
            })
    return variants

def process_variant(variant, server, headers, assembly, sleep_time=1):
    data = {
        "variants": [f"{variant['chrom']} {variant['pos']} . {variant['ref']} {variant['alt']} . . ."],
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

def fetch_frequencies(rsid):
    try:
        handle = Entrez.efetch(db="snp", id=rsid.replace('rs', ''), rettype="docsum", retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        return record
    except Exception as e:
        print(f"Error fetching frequencies for rsID {rsid}: {e}")
        return None

def parse_frequencies(record, variant):
    frequencies = {}
    if record and 'DocumentSummarySet' in record:
        for docsum in record['DocumentSummarySet']['DocumentSummary']:
            chrpos = docsum['CHRPOS'].split(':')[1] if 'CHRPOS' in docsum else ''
            if ('CHR' in docsum and docsum['CHR'] == variant['chrom'] and
                chrpos.isdigit() and int(chrpos) == int(variant['pos']) and
                'REF' in docsum and docsum['REF'] == variant['ref'] and
                'ALT' in docsum and docsum['ALT'] == variant['alt']):
                if 'GLOBAL_MAFS' in docsum:
                    for maf in docsum['GLOBAL_MAFS']:
                        if 'STUDY' in maf and 'FREQ' in maf:
                            study = maf['STUDY']
                            freq = maf['FREQ']
                            allele_freq = freq.split('=')[1].split('/')[0]
                            frequencies[study] = float(allele_freq)
    return frequencies

def fetch_existing_variant_from_ncbi(chrom, pos, ref, alt):
    try:
        query = f"{chrom}[CHR] AND {pos}[CHRPOS]"
        handle = Entrez.esearch(db="snp", term=query)
        record = Entrez.read(handle)
        handle.close()
        if record['IdList']:
            for rsid in record['IdList']:
                handle = Entrez.efetch(db="snp", id=rsid, rettype="docsum", retmode="xml")
                details = Entrez.read(handle)
                handle.close()
                for docsum in details['DocumentSummarySet']['DocumentSummary']:
                    print("docsum", docsum)
                    chrpos = docsum['CHRPOS'].split(':')[1] if 'CHRPOS' in docsum else ''
                    ref_allele = docsum['DOCSUM'].split('|SEQ=[')[-1].split('/')[0]
                    alt_allele = docsum['DOCSUM'].split('|SEQ=[')[-1].split('/')[1].split(']')[0]
                    print("ref_allele", ref_allele, "alt_allele", alt_allele)
                    print("ref", ref, "alt", alt)
                    if (ref_allele == ref and alt_allele == alt):
                        return f"rs{rsid}"
    except Exception as e:
        print(f"Error fetching existing variant from NCBI for {chrom}:{pos}: {e}")
    return "Na"


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
            AF REAL,
            AF_info TEXT,
            genes TEXT,
            consequences TEXT,
            existing_variant TEXT
        )
    ''')
    conn.commit()
    conn.close()

def process_and_insert_results(result, original_variant, cursor):
    if 'error' in result:
        print(f"Failed to process variant: {original_variant}")
        return
    
    variant_data = result[0]
    frequencies, genes, consequences, existing_variant = parse_vep_data(variant_data)
    
    if existing_variant == "Na":
        existing_variant = fetch_existing_variant_from_ncbi(
            original_variant['chrom'], original_variant['pos'],
            original_variant['ref'], original_variant['alt']
        )
        print(f"Existing variant not found in VEP data. Fetching from NCBI: {existing_variant}")
    
    if not frequencies and existing_variant != "Na":
        ncbi_response = fetch_frequencies(existing_variant)
        if ncbi_response:
            frequencies = parse_frequencies(ncbi_response, original_variant)

    # Determine the source of the frequency
    af_info = "Na"
    AF = 0.0

    frequency_sources = ['gnomadg', 'AF', 'ALFA', 'TOPMED']

    for source in frequency_sources:
        if source in frequencies and frequencies[source] != 0.0:
            AF = frequencies[source]
            af_info = source
            break
    else:
        # If no known source is found, take the first available frequency
        AF = next(iter(frequencies.values()), "Na")
        af_info = next(iter(frequencies.keys()), "Na")

    info_dict = parse_info(original_variant['info'])
    variant_type = info_dict.get('TYPE', '')
    variant_len = info_dict.get('LEN', '')
    operation = f"{variant_type}:{variant_len}" if variant_type and variant_len else None
    transcript_ref = info_dict.get('TRANSCRIPT', '')
    transcript_pos = info_dict.get('TRANSCRIPT_POS', '')
    genomic_ref = info_dict.get('GENOME_REF', '')  
    
    cursor.execute('''
        INSERT INTO variants (
            chrom, pos, ref, alt, qual, filter, genomic_ref,
            operation, transcript_ref, transcript_pos, AF, AF_info,
            genes, consequences, existing_variant
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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
        AF,
        af_info,
        ','.join(genes),
        ','.join(consequences),
        existing_variant
    ))

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
    
    db_name = os.path.join(dbs_dir, f'{vcf_name}_indels_variant.db')
    
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
    
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    
    for variant in vcf_data:
        result = process_variant(variant, server, headers, assembly, sleep_time)
        process_and_insert_results(result, variant, cursor)
    
    conn.commit()
    conn.close()
    
    print(f'Data from {vcf_file} has been processed and inserted into {db_name}')

if __name__ == "__main__":
    main()
