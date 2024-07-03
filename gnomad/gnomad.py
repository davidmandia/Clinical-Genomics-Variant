# pip install requests
import requests
import time

# GraphQL endpoint
url = "https://gnomad.broadinstitute.org/api"

# The GraphQL query template
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

# Function to query for a single variant
def query_variant(variant_id, dataset_id="gnomad_r4", retries=5):
    # Query variables
    variables = {"variantId": variant_id, "datasetId": dataset_id}

    # Retry mechanism
    for attempt in range(retries):
        try:
            # HTTP POST request
            response = requests.post(url, json={"query": query, "variables": variables})
            if response.status_code == 200:
                return response.json()  # Returns the JSON response
            else:
                raise Exception(
                    f"Query failed to run by returning code of {response.status_code}. {response.text}"
                )
        except Exception as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            time.sleep(5)  # Wait for 5 seconds before retrying

    raise Exception("Max retries exceeded")

# List of variants to query
variants = ["1-55052746-GT-G", "1-55058620-TG-T"]  # Add your variants here

# Results list
results = []

# Loop through each variant and query
for variant in variants:
    try:
        result = query_variant(variant)
        results.append(result)
        time.sleep(6)  # Sleep to respect the 10 queries per minute limit
    except Exception as e:
        print(f"Failed to retrieve data for variant {variant}: {e}")

# results now contains the response for each variant
print(results)