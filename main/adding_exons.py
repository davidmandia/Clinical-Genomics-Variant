import pandas as pd

# Load the exon data from the CSV file
exon_data = pd.read_csv('exons.csv')

# Filter rows where 'Exon ID' is not NaN
exon_with_ids = exon_data.dropna(subset=['Exon ID'])

# Display the first few rows of this filtered data
print(exon_with_ids.head())
