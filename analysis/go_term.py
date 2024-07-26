import pandas as pd
import requests
import json
import matplotlib.pyplot as plt

def enrichr_query(genes, description="Gene list"):
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr'
    ADDLIST_URL = f'{ENRICHR_URL}/addList'
    ENRICH_URL = f'{ENRICHR_URL}/enrich'
    headers = {'Content-Type': 'application/json'}

    genes_str = '\n'.join(genes)
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }
    response = requests.post(ADDLIST_URL, files=payload)
    if not response.ok:
        raise Exception(f"Error posting gene list: {response.text}")
    
    user_list_id = json.loads(response.text)['userListId']

    payload = {
        'userListId': user_list_id,
        'backgroundType': 'GO_Biological_Process_2021'
    }
    response = requests.get(ENRICH_URL, params=payload, headers=headers)
    if not response.ok:
        raise Exception(f"Error fetching enrichment results: {response.text}")

    return response.json()

def main():
    gene_counts = pd.read_csv("green_clinically_relevant_gene_counts.csv")
    gene_list = gene_counts['Gene'].tolist()
    
    print(f"Submitting {len(gene_list)} genes to Enrichr for GO term enrichment analysis.")
    
    enrichment_results = enrichr_query(gene_list)
    
    enriched_terms = enrichment_results['GO_Biological_Process_2021']
    results_df = pd.DataFrame(enriched_terms, columns=[
        'Rank', 'Term', 'P-value', 'Z-score', 'Combined Score', 'Genes', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value'
    ])
    
    results_df.to_csv("go_enrichment_results.csv", index=False)
    print("Enrichment results saved to go_enrichment_results.csv")
    
    #Combined Score=ln(P-value)×Z
    
    # Visualization: Plot top 10 GO terms by combined score
    top_terms = results_df.head(10)
    plt.figure(figsize=(10, 6))
    plt.barh(top_terms['Term'], top_terms['Combined Score'], color='skyblue')
    plt.xlabel('Combined Score')
    plt.ylabel('GO Term')
    plt.title('Top 10 Enriched GO Terms')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig("top_enriched_go_terms.png")
    plt.show()

if __name__ == "__main__":
    main()
