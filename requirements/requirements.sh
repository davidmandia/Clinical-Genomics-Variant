pip install requests pandas matplotlib
# Should be already included with Python 
pip install sqlite3

## For testing 
pip install pytest
sudo apt update


pip install biopython
pip install bcbio-gff

sqlite3 output/database/GRCh38_indels_variant.db .dump > GRCh38_indels_variant.sql
