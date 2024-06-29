sudo apt-get update
sudo apt-get install ncbi-blast+
makeblastdb -version


wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
gzip -d GRCh38_latest_genomic.fna.gz
makeblastdb -in GRCh38_latest_genomic.fna -out GRCh38 -parse_seqids -dbtype nucl