import requests

def get_ncbi_sequence(accession):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text"
    }
    response = requests.get(base_url, params=params)

    if response.status_code == 200:
        return response.text
    else:
        response.raise_for_status()

if __name__ == "__main__":
    accession = "NC_000001.10"
    try:
        sequence = get_ncbi_sequence(accession)
        print(f"Reference sequence for {accession}:\n{sequence}")
    except Exception as e:
        print(f"Error fetching reference sequence: {e}")
