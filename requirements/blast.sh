#!/bin/bash

# Update the package list and install ncbi-blast+
sudo apt-get update
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

# Download and process the GRCh37 genome file
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz -P data
gzip -d data/GRCh37_latest_genomic.fna.gz
makeblastdb -in data/GRCh37_latest_genomic.fna -out input/blastdb/GRCh37 -parse_seqids -dbtype nucl
rm data/GRCh37_latest_genomic.fna.gz
rm data/GRCh37_latest_genomic.fna

# Remove the data directory
rmdir data

echo "BLAST database creation completed successfully."
