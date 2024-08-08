sudo apt update
pip install requests pandas matplotlib biopython bcbio-gff
# Should be already included with Python 
pip install sqlite3

## For testing 
pip install pytest


# For AWs deployment 
sqlite3 output/database/GRCh38_indels_variant.db .dump > data_transfer_AWS/GRCh38_indels_variant.sql
pip install psycopg2

sqlite3 output/database/GRCh37_indels_variant.db .dump > data_transfer_AWS/GRCh37_indels_variant.sql
## However, the conversion to csv , amend to Postgres, then deploy through CSV is way faster 

# sqlite3 /workspaces/Clinical-Genomics-Variant/output/database/GRCh38_indels_variant.db <<EOF
# .headers on
# .mode csv
# .output /workspaces/Clinical-Genomics-Variant/output/database/GRCh38_indels_variant.csv
# SELECT * FROM variants;
# EOF



