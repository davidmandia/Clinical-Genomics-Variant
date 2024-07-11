#Bash to run multiple script at the same time 


python main.py input/GRCh37.sam --pseudo False --db input/blastdb/GRCh37
python vcf_to_database.py output/VCFs/GRCh37_False_indels.vcf --assembly GRCh37



python main.py input/GRCh38.sam --pseudo False --db input/blastdb/GRCh38
python vcf_to_database.py output/VCFs/GRCh37_False_indels.vcf --assembly GRCh37


