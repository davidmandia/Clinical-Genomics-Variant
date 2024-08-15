# read the VCF files of GRCh37 and GRCh38 from output
Extract unique variant calls from the data. Each variant call is constructed using source,the chromosome, position, reference allele, and alternate allele.

# query data from VariantValidator API
The API request URL is constructed dynamically based on the genome build and the variant description. The URL follows a specific pattern(https://rest.variantvalidator.org/LOVD/lovd/{genome_build}/{variant_description}/{transcript_model}/{select_transcripts}/{checkonly}/{liftover}?content-type=application/json) 

# Using filters to select variants
HGVS strings are categorized based on their prefixes: those starting with 'E' are classified as Ensembl variants, while those starting with 'N' are classified as RefSeq variants, focusing specifically on entries starting with "c" (coding sequence)

Then add filters to select 5 prime UTR variants, 3 prime UTR variants, gap in refseq and ensemble, gap in refseq not in ensemble, gap not in refseq but in ensemble and no gap. After that, count the variants in different files

