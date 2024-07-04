import requests
import time

# GraphQL endpoint
url = "https://gnomad.broadinstitute.org/api"

# The GraphQL query template (unchanged)
query = """
query GnomadVariant($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variant_id
    reference_genome
    chrom
    pos
    ref
    alt
    colocated_variants
    coverage {
      exome {
        mean
        over_20
      }
      genome {
        mean
        over_20
      }
    }
    exome {
      ac
      an
      ac_hemi
      ac_hom
      filters
      populations {
        id
        ac
        an
        ac_hemi
        ac_hom
      }
    }
    genome {
      ac
      an
      ac_hemi
      ac_hom
      filters
      populations {
        id
        ac
        an
        ac_hemi
        ac_hom
      }
    }
    flags
    lof_curations {
      gene_id
      gene_symbol
      verdict
      flags
      project
    }
    rsids
    transcript_consequences {
      domains
      gene_id
      gene_version
      gene_symbol
      hgvs
      hgvsc
      hgvsp
      is_canonical
      is_mane_select
      is_mane_select_version
      lof
      lof_flags
      lof_filter
      major_consequence
      polyphen_prediction
      sift_prediction
      transcript_id
      transcript_version
    }
    in_silico_predictors {
      id
      value
      flags
    }
  }
}
"""

# Function to query for a single variant (unchanged)
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

# Function to parse VCF file and extract variant IDs (plain Python)
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
    vcf_file = "/workspaces/Clinical-Genomics-Variant/output/VCFs/sample_False_indels.vcf"
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
