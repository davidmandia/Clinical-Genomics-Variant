import pysam

# Input and output file paths
input_vcf = "/workspaces/Clinical-Genomics-Variant/NA12878.vcf"
output_vcf = "/workspaces/Clinical-Genomics-Variant/NA12878_indels.vcf"

# Open input VCF file
vcf_in = pysam.VariantFile(input_vcf, "r")

# Create an output VCF file with the same header
vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)

# Loop through each record in the input VCF file
for record in vcf_in:
    # Check if the length of REF and ALT alleles are different
    if any(len(record.ref) != len(alt) for alt in record.alts):
        # Write the record to the output VCF file
        vcf_out.write(record)

# Close the VCF files
vcf_in.close()
vcf_out.close()

print(f"Filtered VCF with indels saved to: {output_vcf}")
