import csv

def extract_exons_from_gff(gff_file, output_csv):
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

    # Write to CSV
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Chromosome', 'Start', 'End', 'Gene ID', 'Transcript ID', 'Exon ID'])
        writer.writerows(exons)

# Example usage:
extract_exons_from_gff('gff_data/GCF_000001405.40-RS_2023_03_genomic.gff', 'exons.csv')
