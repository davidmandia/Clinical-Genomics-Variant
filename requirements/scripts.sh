#Bash to run multiple script at the same time 

##GRCh37 database build
python main/sam_to_vcf.py input/GRCh37.sam --pseudo False --db input/blastdb/GRCh37
python main/vcf_to_database.py output/VCFs/GRCh37_False_indels.vcf --assembly GRCh37
python main/adding_gene.py output/database/GRCh37_indels_variant.db gff_data/GCF_000001405.25_GRCh37.p13_genomic.gff
python main/check_clinical_relevance.py output/database/GRCh37_indels_variant.db

## GRCh38 database build

python main/sam_to_vcf.py input/GRCh38.sam --pseudo False --db input/blastdb/GRCh38
python main/vcf_to_database.py output/VCFs/GRCh38_False_indels.vcf --assembly GRCh38
python main/adding_gene.py output/database/GRCh38_indels_variant.db gff_data/GCF_000001405.40-RS_2023_03_genomic.gff
python main/check_clinical_relevance.py output/database/GRCh38_indels_variant.db




