sudo apt update
pip install requests pandas matplotlib biopython bcbio-gff
# Should be already included with Python 
pip install sqlite3

## For testing 
pip install pytest


# For AWs deployment 
sqlite3 output/database/GRCh38_indels_variant.db .dump > data_transfer_AWS/GRCh38_indels_variant.sql
pip install psycopg2


