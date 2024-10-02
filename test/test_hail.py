import hail as hl

# Access gnomAD data from AWS S3
gnomad = hl.read_table('s3://gnomad-public-us-east-1/release/4.0/ht/genomes/gnomad.genomes.r4.0.sites.ht')

# Load your local VCF file
mt = hl.import_vcf('sorted_sample_GRCh38_variants.vcf')

# Annotate with gnomAD allele frequencies
mt = mt.annotate_rows(gnomad_AF=gnomad[mt.locus, mt.alleles].info.AF)

# Write the annotated VCF
mt.rows().select('gnomad_AF').export('annotated_with_gnomad_AF.tsv')
