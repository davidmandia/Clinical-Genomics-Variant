import requests
import time

# GraphQL endpoint
url = "https://gnomad.broadinstitute.org/api"

# The GraphQL query template
query = """
query GnomadVariant($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variant_id
    chrom
    pos
    ref
    alt
    exome {
      ac
      an
      ac_hom
      populations {
        id
        ac
        an
        ac_hom
      }
    }
    genome {
      ac
      an
      ac_hom
      populations {
        id
        ac
        an
        ac_hom
      }
    }
  }
}
"""

# Function to query for a single variant
def query_variant(variant_id, dataset_id="gnomad_r4", retries=5):
    variables = {"variantId": variant_id, "datasetId": dataset_id}
    for attempt in range(retries):
        try:
            response = requests.post(url, json={"query": query, "variables": variables})
            if response.status_code == 200:
                return response.json()
            else:
                raise Exception(f"Query failed with status code {response.status_code}: {response.text}")
        except Exception as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            time.sleep(5)

    raise Exception("Max retries exceeded")

# Function to parse VCF file and extract variant IDs
def parse_vcf(vcf_file):
    variant_ids = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue  # Skip header lines
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            variant_id = f"{chrom}-{pos}-{ref}-{alt}"
            variant_ids.append(variant_id)
    return variant_ids

# Example usage
if __name__ == "__main__":
    vcf_file = "output/VCFs/sample_False_indels.vcf"
    variants = parse_vcf(vcf_file)
    results = []

    for variant in variants:
        try:
            result = query_variant(variant)
            results.append(result)
            time.sleep(6)  # Respect rate limit
        except Exception as e:
            print(f"Failed to retrieve data for variant {variant}: {e}")

    # Save results to a text or CSV file
    with open("gnomad_results.txt", "w") as f:
        for result in results:
            f.write(str(result) + "\n")

    print("Query results saved to gnomad_results.txt")
