sudo apt-get install -y ncbi-blast+

# Verify the installation
makeblastdb -version

# Create necessary directories
mkdir -p input/blastdb
mkdir -p data

# Download and process the GRCh38 genome file
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz -P data
gzip -d data/GRCh38_latest_genomic.fna.gz
makeblastdb -in data/GRCh38_latest_genomic.fna -out input/blastdb/GRCh38 -parse_seqids -dbtype nucl
rm data/GRCh38_latest_genomic.fna.gz
rm data/GRCh38_latest_genomic.fna