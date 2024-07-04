import requests

# Define the genomic coordinates and alleles (replace with your specific values)
chromosome = "1"
position = "1000000"
reference_allele = "A"
alternative_allele = "G"

# Construct the URL for the Ensembl API request
url = f"https://rest.ensembl.org/vep/human/region/{chromosome}:{position}/{reference_allele}/{alternative_allele}?" \
      f"content-type=application/json;include=population_frequencies"

# Send the GET request to the Ensembl API
response = requests.get(url)

print(url)

# Check if the response is successful
if response.status_code == 200:
    # Parse the JSON response
    variant_data = response.json()
    print(variant_data)

    # Extract allele frequency information
    if 'Variation' in variant_data:
        for allele in variant_data['Variation']:
            if 'PopulationAlleleFrequencies' in allele:
                frequencies = allele['PopulationAlleleFrequencies']
                for freq in frequencies:
                    print(f"Population: {freq['population']} - Allele Frequency: {freq['frequency']}")

else:
    print(f"Failed to retrieve data. HTTP Status Code: {response.status_code}")
