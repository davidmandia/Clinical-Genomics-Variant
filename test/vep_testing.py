import requests
import json

def get_vep_data_for_variant(chrom, pos, ref, alt, assembly='GRCh37'):
    """
    Fetches variant data from the Ensembl VEP API.

    Args:
        chrom (str): Chromosome of the variant.
        pos (int): Position of the variant.
        ref (str): Reference allele.
        alt (str): Alternate allele.
        assembly (str): Genome assembly (default: 'GRCh37').

    Returns:
        dict: The JSON response from the VEP API.
    """
    if assembly == 'GRCh38':
        server = "https://rest.ensembl.org"
    elif assembly == 'GRCh37':
        server = "https://grch37.rest.ensembl.org"
        
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    variant = f"{chrom} {pos} . {ref} {alt} . . ."
    data = {
        "variants": [variant],
        "variant_class": True,
        "allele_number": True,
        "regulatory": True,
        "canonical": True,
        "check_existing": True,
        "assembly": assembly
    }
    
    response = requests.post(f"{server}/vep/human/region", headers=headers, json=data)
    
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Error {response.status_code}: {response.text}")
        return None

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

def main():
    # Define the variant
    chrom = "1"
    pos = 145500921
    ref = "C"
    alt = "CTA"
    assembly = "GRCh37"
    
    # Get the VEP data for the variant
    vep_data = get_vep_data_for_variant(chrom, pos, ref, alt, assembly)
    
    # Print the response and choose the best rsID
    if vep_data:
        print(json.dumps(vep_data, indent=4))
        for entry in vep_data:
            best_rsid = choose_best_rsid(entry.get('colocated_variants', []))
            if best_rsid:
                print(f"Chosen rsID: {best_rsid}")

if __name__ == "__main__":
    main()
