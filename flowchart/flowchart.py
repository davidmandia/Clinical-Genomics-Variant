from graphviz import Digraph

# Create a new directed graph
flowchart = Digraph("Methods_Flowchart", format="png")
flowchart.attr(rankdir='TB', size='8,10', fontname='Helvetica', fontsize='12')

# Define nodes with simple shapes and no excessive colors
flowchart.node('A', 'Retrieval of RefSeq transcript Alignments (GRCh38, GRCh37)', shape='box')
flowchart.node('B', 'Identify Indels (< 3bp)', shape='box')
flowchart.node('C', 'Generate VCF Files using NCBI-BLAST for reference sequences', shape='box')
flowchart.node('G', 'Parse VCF and Query Ensembl VEP API', shape='box')
flowchart.node('H', 'Store Allele Frequencies and Variant Data in SQLite Database', shape='box')
flowchart.node('I', 'Extract Genes and Query NHS PanelApp', shape='box')
flowchart.node('J', 'Update Database with Clinical Relevance', shape='box')
flowchart.node('K', 'VCF Analyzer Tool', shape='box')
flowchart.node('V', 'VCFs', shape='diamond')
flowchart.node('L1', 'Check if VCFs are in Database', shape='box')
flowchart.node('L2', 'Return Expected Variants', shape='box')
flowchart.node('L3', 'Return Unexpected Variants', shape='box')
flowchart.node('M', 'Generate Reports', shape='box')

# Define edges with logic for expected/unexpected variants
flowchart.edge('A', 'B')
flowchart.edge('B', 'C')
flowchart.edge('C', 'G')
flowchart.edge('G', 'H')
flowchart.edge('H', 'I')
flowchart.edge('I', 'J')
flowchart.edge('J', 'K')
flowchart.edge('V', 'K')
flowchart.edge('K', 'L1')
flowchart.edge('L1', 'L2', label='Present')
flowchart.edge('L1', 'L3', label='Not Present')
flowchart.edge('L2', 'M')
flowchart.edge('L3', 'M')

# Render and save the flowchart
flowchart.render("Methods_Flowchart_VCF_Variant_Calls_Logic", view=True)
