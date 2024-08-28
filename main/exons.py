import pandas as pd

def extract_exons_from_gff(gff_file):
    exons = []

    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith('#'):
                continue
            
            columns = line.strip().split('\t')
            if len(columns) < 9:
                continue
            
            seqid = columns[0]
            source = columns[1]
            feature_type = columns[2]
            start = int(columns[3])
            end = int(columns[4])
            score = columns[5]
            strand = columns[6]
            phase = columns[7]
            attributes = columns[8]

            # Filter for exons
            if feature_type.lower() == 'exon':
                # Parse attributes for additional info
                attr_dict = {}
                for attribute in attributes.split(';'):
                    key_value = attribute.strip().split('=')
                    if len(key_value) == 2:
                        attr_dict[key_value[0]] = key_value[1]
                
                gene_id = attr_dict.get('gene_id', 'N/A')
                transcript_id = attr_dict.get('transcript_id', 'N/A')
                exon_id = attr_dict.get('ID', 'N/A')  # Capturing exon ID
                
                exons.append([seqid, start, end, gene_id, transcript_id, exon_id])

    # Create a pandas DataFrame
    exons_df = pd.DataFrame(exons, columns=['Chromosome', 'Start', 'End', 'Gene ID', 'Transcript ID', 'Exon ID'])
    return exons_df

# Example 
exons_df = extract_exons_from_gff("gff_data/GCF_000001405.40-RS_2023_03_genomic.gff")
exons_df.to_csv("main/exons_GRCh38.csv")
print(exons_df.head())
