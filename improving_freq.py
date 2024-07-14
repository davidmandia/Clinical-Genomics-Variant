import requests
import json

def fetch_vep_data():
    server = "https://rest.ensembl.org"
    endpoint = "/vep/human/region"
    headers = {"Content-Type": "text/x-vcf", "Accept": "application/json"}

    # VCF content for the variant
    vcf_content = """##fileformat=VCFv4.2
##source=EnsemblVEP
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t8355038\t.\tTC\tT\t.\t.\t.\n"""

    response = requests.post(f"{server}{endpoint}", headers=headers, data=vcf_content)
    
    if response.status_code != 200:
        print(f"Error: API request failed with status code {response.status_code}")
        return None
    
    return response.json()

def main():
    vep_response = fetch_vep_data()
    
    if vep_response:
        print(json.dumps(vep_response, indent=2))
    else:
        print("Failed to retrieve data.")

if __name__ == "__main__":
    main()
