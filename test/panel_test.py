import requests

def get_color_from_evidence(evidence_list):
    for evidence in evidence_list:
        if 'Green' in evidence:
            return 'Green'
        elif 'Amber' in evidence:
            return 'Amber'
        elif 'Red' in evidence:
            return 'Red'
    return 'No color label found'

def test_api_for_genes(gene_list, api_url_base):
    for gene in gene_list:
        api_url = f"{api_url_base}{gene}"
        print(f"Testing API for gene: {gene}")
        print("api_url:", api_url)
        response = requests.get(api_url)
        if response.status_code == 200:
            data = response.json()
            if data["count"] > 0:
                gene_data = data["results"][0]
                evidence = gene_data.get("evidence", [])
                color_label = get_color_from_evidence(evidence)
                print(f"API is working correctly for gene: {gene}. Here's the color label found in the evidence:")
                print("Color label:", color_label)
            else:
                print(f"No results found for gene: {gene}")
        else:
            print(f"API request failed for gene: {gene} with status code: {response.status_code}")
            print(response.text)

# Example list of genes
gene_list = ["MAPT", "BRCA1", "BRCA2", "TP53", "EGFR"]

# Base API URL
api_url_base = "https://panelapp.genomicsengland.co.uk/api/v1/genes/"

test_api_for_genes(gene_list, api_url_base)
