import csv
import re
import requests

# Install requests if not already installed
!pip install requests

# Define the public URL of GRch 38 S3 file
public_url = 'https://rp2clinicalgenomicsdavidmandia.s3.amazonaws.com/output_GCF_000001405.40-RS_2023_03_knownrefseq_alns.sam'
local_file_name = 'GRCh38.sam'

# Download the file
response = requests.get(public_url)

# Check if the request was successful
if response.status_code == 200:
    with open(local_file_name, 'wb') as f:
        f.write(response.content)
    print(f"File downloaded successfully and saved as {local_file_name}")
else:
    print(f"Failed to download file. Status code: {response.status_code}")

# Define the public URL of your GRCh37 S3 file
public_url_2 = 'https://rp2clinicalgenomicsdavidmandia.s3.amazonaws.com/outputGCF_000001405.25_GRCh37.p13_knownrefseq_alns.sam'
local_file_name_2 = 'GRCh37.sam'

# Download the file
response_2 = requests.get(public_url_2)

# Check if the request was successful
if response_2.status_code == 200:
    with open(local_file_name_2, 'wb') as f:
        f.write(response_2.content)
    print(f"File downloaded successfully and saved as {local_file_name_2}")
else:
    print(f"Failed to download file. Status code: {response_2.status_code}")