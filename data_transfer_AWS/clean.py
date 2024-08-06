import pandas as pd

# Load the CSV file
csv_file_path = "/workspaces/Clinical-Genomics-Variant/output/database/GRCh38_indels_variant.csv"
df = pd.read_csv(csv_file_path)

# Replace "Na" with None (which will be interpreted as NULL by PostgreSQL)
df.replace("Na", None, inplace=True)

# Save the cleaned CSV
clean_csv_file_path = "/workspaces/Clinical-Genomics-Variant/output/database/GRCh38_indels_variant_clean.csv"
df.to_csv(clean_csv_file_path, index=False)
