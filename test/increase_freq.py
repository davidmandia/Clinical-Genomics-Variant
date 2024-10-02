import xml.etree.ElementTree as ET
from Bio import Entrez

def parse_xml_response(xml_response):
    print("xml_response", xml_response)
    
    # Parse the XML response
    root = ET.fromstring(xml_response)
    
    # Define the namespace for the XML
    ns = {'ns': 'https://www.ncbi.nlm.nih.gov/SNP/docsum'}

    # Find the DocumentSummary element
    document = root.find('ns:DocumentSummary', ns)
    
    # Check if document is None (i.e., the response doesn't contain relevant data)
    if document is None:
        print("No document found in response. Likely empty or malformed response.")
        return
    
    print("document", document)
    
    # Extract SNP ID
    snp_id = document.find('ns:SNP_ID', ns).text
    print(f"SNP ID: {snp_id}")
    
    # Extract allele frequencies (GLOBAL_MAFS)
    mafs = document.find('ns:GLOBAL_MAFS', ns)
    if mafs is not None:
        for maf in mafs.findall('ns:MAF', ns):
            study = maf.find('ns:STUDY', ns).text
            freq = maf.find('ns:FREQ', ns).text
            print(f"Study: {study}, Frequency: {freq}")
    else:
        print("No allele frequencies (GLOBAL_MAFS) found.")
    
    # Extract gene information
    genes = document.find('ns:GENES', ns)
    if genes is not None:
        for gene in genes.findall('ns:GENE_E', ns):
            gene_name = gene.find('ns:NAME', ns).text
            gene_id = gene.find('ns:GENE_ID', ns).text
            print(f"Gene: {gene_name}, Gene ID: {gene_id}")
    else:
        print("No gene information found.")

def get_variant_info_dbSNP(rsid):
    Entrez.email = "david.mandia@yahoo.com"
    
    # Query dbSNP for the given RSID
    handle = Entrez.efetch(db="snp", id=rsid, rettype="xml")
    response = handle.read()
    handle.close()
    
    # Parse the XML response
    parse_xml_response(response)

if __name__ == "__main__":
    # Use the RSID for your variant
    rsid = "rs1569899370"
    
    # Query dbSNP for allele frequency information
    get_variant_info_dbSNP(rsid)
