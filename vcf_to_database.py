import sys
import requests
import sqlite3
import argparse
import os

def get_vep_data(chrom, pos, ref, alt):
    """Fetch data from Ensembl VEP REST API using VCF format"""
    server = "https://rest.ensembl.org"
    ext = "/vep/human/region"
    chrom = chrom.replace('chr', '')
    variant = f"{chrom}:{pos}:{ref}:{alt}"
    
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
        "hgvs": True,
        "transcripts": True
    }
    response = requests.post(server + ext, headers=headers, json=data)
    if response.status_code != 200:
        print(f"Failed to retrieve data for {variant}: {response.status_code}")
        return None
    return response.json()

def parse_vep_data(data):
    """Parse VEP API response for allele frequencies, SIFT score, gene, and transcript"""
    if not data or len(data) == 0:
        return {}, None, None, None, None
    
    variant_data = data[0]
    frequencies = {}
    sift_score = None
    gene = None
    transcript = None
    
    if 'colocated_variants' in variant_data:
        for cv in variant_data['colocated_variants']:
            if 'frequencies' in cv:
                for pop, freq in cv['frequencies'].items():
                    if pop.startswith('gnomad'):
                        frequencies[pop] = freq
    
    if 'transcript_consequences' in variant_data:
        for tc in variant_data['transcript_consequences']:
            if 'sift_score' in tc:
                sift_score = tc['sift_score']
            if 'gene_symbol' in tc:
                gene = tc['gene_symbol']
            if 'transcript_id' in tc:
                transcript = tc['transcript_id']
            if sift_score and gene and transcript:
                break
    
    hgvs = variant_data.get('id', '')
    return frequencies, sift_score, hgvs, gene, transcript

def create_database(db_name):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS variants (
            chrom TEXT,
            pos INTEGER,
            id TEXT,
            ref TEXT,
            alt TEXT,
            qual REAL,
            filter TEXT,
            genomic_ref TEXT,
            hgvs TEXT,
            operation TEXT,
            transcript_ref TEXT,
            transcript_pos TEXT,
            gene TEXT,
            transcript TEXT,
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
            sift_score REAL
        )
    ''')
    conn.commit()
    conn.close()

def parse_info(info_str):
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
    return info_dict

def process_vcf_and_insert(input_file, db_name):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    
    with open(input_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue  # Skip header lines
            
            fields = line.strip().split('\t')
            chrom, pos, id_, ref, alt, qual, filter_, info = fields[:8]
            
            print(f"Processing variant: {chrom}:{pos}:{ref}>{alt}")
            
            vep_data = get_vep_data(chrom, pos, ref, alt)
            frequencies, sift_score, hgvs, gene, transcript = parse_vep_data(vep_data)
            
            info_dict = parse_info(info)
            
            # Extract operation (TYPE and LEN)
            variant_type = info_dict.get('TYPE', '')
            variant_len = info_dict.get('LEN', '')
            operation = f"{variant_type}:{variant_len}" if variant_type and variant_len else None
            
            # Extract transcript reference and position
            transcript_ref = info_dict.get('TRANSCRIPT', '')
            transcript_pos = info_dict.get('TRANSCRIPT_POS', '')
            
            # Extract genomic reference
            genomic_ref = info_dict.get('GENOME_REF', '')
            
            gnomad_fields = [
                'gnomad', 'gnomad_eas', 'gnomad_nfe', 'gnomad_fin',
                'gnomad_amr', 'gnomad_afr', 'gnomad_asj', 'gnomad_oth',
                'gnomad_sas', 'gnomad_mid', 'gnomad_ami'
            ]
            
            gnomad_values = [frequencies.get(field) for field in gnomad_fields]
            
            cursor.execute('''
                INSERT INTO variants (
                    chrom, pos, id, ref, alt, qual, filter, genomic_ref, hgvs,
                    operation, transcript_ref, transcript_pos, gene, transcript,
                    gnomadg, gnomadg_eas, gnomadg_nfe, gnomadg_fin,
                    gnomadg_amr, gnomadg_afr, gnomadg_asj, gnomadg_oth,
                    gnomadg_sas, gnomadg_mid, gnomadg_ami, sift_score
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                chrom, int(pos), id_, ref, alt, float(qual) if qual != '.' else None,
                filter_, genomic_ref, hgvs, operation, transcript_ref, transcript_pos,
                gene, transcript, *gnomad_values, sift_score
            ))
    
    conn.commit()
    conn.close()

def main():
    parser = argparse.ArgumentParser(description="Process a VCF file, fetch VEP data, and store in a SQLite database.")
    parser.add_argument('vcf', help="Path to the input VCF file")
    args = parser.parse_args()

    vcf_file = args.vcf
    vcf_name = os.path.basename(vcf_file).split('.')[0]
    
    output_dir = "output"
    dbs_dir = os.path.join(output_dir, "database")
    os.makedirs(dbs_dir, exist_ok=True)
    
    db_name = os.path.join(dbs_dir, f'{vcf_name}_variant.db')
    
    # Remove existing database if it exists
    if os.path.exists(db_name):
        os.remove(db_name)
    
    create_database(db_name)
    process_vcf_and_insert(vcf_file, db_name)
    print(f'Data from {vcf_file} has been processed and inserted into {db_name}')

if __name__ == '__main__':
    main()