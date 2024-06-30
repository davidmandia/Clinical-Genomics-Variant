Bash to run multiple script at the same time 

python main.py input/GRCh37.sam --pseudo True --db input/blastdb/GRCh37
python database_operations_variantions.py output/VCFs/GRCh37_True_indels.vcf 

python main.py input/GRCh37.sam --pseudo False --db input/blastdb/GRCh37
python database_operations_variantions.py output/VCFs/GRCh37_False_indels.vcf 

python main.py input/GRCh38.sam --pseudo True --db input/blastdb/GRCh38
python database_operations_variantions.py output/VCFs/GRCh38_True_indels.vcf 

python main.py input/GRCh38.sam --pseudo False --db input/blastdb/GRCh38
python database_operations_variantions.py output/VCFs/GRCh38_False_indels.vcf  


