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
        "canonical": True
    }
    response = requests.post(server + ext, headers=headers, json=data)
    if response.status_code != 200:
        print(f"Failed to retrieve data for {variant}: {response.status_code}")
        return None
    return response.json()

# return frequencies, genes, consequences
def parse_vep_data(data):
    """Parse VEP API response for allele frequencies, genes, and consequences"""
    if not data or len(data) == 0:
        return {}, [], []
    
    variant_data = data[0]
    frequencies = {}
    genes = set()
    consequences = set()
    
    if 'colocated_variants' in variant_data.keys():
        if 'frequencies' in variant_data['colocated_variants'][0].keys():
            gnomad_fre = next(iter(variant_data["colocated_variants"][0]["frequencies"].values()))
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
 
    return frequencies, list(genes), list(consequences)

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
            consequences TEXT
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
            chrom, pos, _, ref, alt, qual, filter_, info = fields[:8]
            
            print(f"Processing variant: {chrom}:{pos}:{ref}>{alt}")
            
            vep_data = get_vep_data(chrom, pos, ref, alt)
            frequencies, genes, consequences = parse_vep_data(vep_data)
            
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
                'gnomadg', 'gnomadg_eas', 'gnomadg_nfe', 'gnomadg_fin',
                'gnomadg_amr', 'gnomadg_afr', 'gnomadg_asj', 'gnomadg_oth',
                'gnomadg_sas', 'gnomadg_mid', 'gnomadg_ami'
            ]
            
            gnomad_values = [frequencies.get(field) for field in gnomad_fields]
            print("gnomad_values",gnomad_values)
            
            cursor.execute('''
                INSERT INTO variants (
                    chrom, pos, ref, alt, qual, filter, genomic_ref,
                    operation, transcript_ref, transcript_pos,
                    gnomadg, gnomadg_eas, gnomadg_nfe, gnomadg_fin,
                    gnomadg_amr, gnomadg_afr, gnomadg_asj, gnomadg_oth,
                    gnomadg_sas, gnomadg_mid, gnomadg_ami, genes, consequences
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                chrom, int(pos), ref, alt, float(qual) if qual != '.' else None,
                filter_, genomic_ref, operation, transcript_ref, transcript_pos,
                *gnomad_values, ','.join(genes), ','.join(consequences)
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