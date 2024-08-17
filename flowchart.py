from graphviz import Digraph

# Create a new directed graph
flowchart = Digraph("Methods_Flowchart", format="png")
flowchart.attr(rankdir='TB', size='8,10', fontname='Helvetica', fontsize='12')

# Define nodes with simple shapes and no excessive colors
flowchart.node('A', 'Download Genome Alignments (GRCh38, GRCh37)', shape='box')
flowchart.node('B', 'Identify Indels (< 3bp)', shape='box')
flowchart.node('C', 'Generate VCF Files using NCBI-BLAST for reference sequences', shape='box')
flowchart.node('G', 'Parse VCF and Query Ensembl VEP API', shape='box')
flowchart.node('H', 'Store Allele Frequencies and Variant Data in SQLite Database', shape='box')
flowchart.node('I', 'Extract Genes and Query NHS PanelApp', shape='box')
flowchart.node('J', 'Update Database with Clinical Relevance', shape='box')
flowchart.node('K', 'Deploy VCF Reader Tool', shape='box')
flowchart.node('L', 'Classify Variants (Expected/Unexpected) along with Supplementary Data', shape='box')
flowchart.node('M', 'Generate Reports', shape='box')

# Define edges
flowchart.edges(['AB', 'BC', 'CG', 'GH', 'HI', 'IJ', 'JK', 'KL', 'LM'])

# Render and save the flowchart
flowchart.render("Methods_Flowchart_Simple_Scientific", view=True)
