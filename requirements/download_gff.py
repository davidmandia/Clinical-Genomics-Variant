import requests
import gzip
import shutil
import os


def download_and_extract_gff(url, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    gff_gz_path = os.path.join(output_dir, "GCF_000001405.40-RS_2023_03_genomic.gff.gz")
    gff_path = os.path.join(output_dir, "GCF_000001405.40-RS_2023_03_genomic.gff")

    # Download the file
    response = requests.get(url, stream=True)
    with open(gff_gz_path, 'wb') as file:
        shutil.copyfileobj(response.raw, file)

    # Extract the file
    with gzip.open(gff_gz_path, 'rb') as f_in:
        with open(gff_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Print the first few lines
    with open(gff_path, 'r') as file:
        for _ in range(10):  # Adjust the number of lines to read as needed
            print(file.readline().strip())

# Example usage
url = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/historical/GRCh38/GCF_000001405.40-RS_2023_03_historical/GCF_000001405.40-RS_2023_03_genomic.gff.gz"
url2 = "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/105.20201022/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz"
output_dir = "gff_data"
download_and_extract_gff(url, output_dir)
download_and_extract_gff(url2, output_dir)
