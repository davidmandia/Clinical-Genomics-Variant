# Bash to run multiple script at the same time 

# python main.py GRCh37.sam --pseudo True
# python main.py GRCh37.sam --pseudo False
# python main.py GRCh38.sam --pseudo True
# python main.py GRCh38.sam --pseudo False
python database_operations_variantions.py outputs/GRCh37_False_indels.vcf 
python database_operations_variantions.py outputs/GRCh37_True_indels.vcf 
python database_operations_variantions.py outputs/GRCh38_False_indels.vcf 
python database_operations_variantions.py outputs/GRCh38_True_indels.vcf 