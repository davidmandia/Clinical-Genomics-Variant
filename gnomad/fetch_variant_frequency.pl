use strict;
use warnings;
use lib '/workspaces/Clinical-Genomics-Variant/ensembl-vep/ensembl/modules';
use Bio::EnsEMBL::Registry;

# Connect to Ensembl API
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',  # Replace with your Ensembl database host
    -user => 'anonymous',              # Anonymous access is common for public databases
);

# Get VariationAdaptor
my $va = $registry->get_adaptor('human', 'variation', 'variation');

# Define variant details (replace with your variant's chromosome, position, ref, alt)
my $chromosome = '1';
my $position = 1000000;
my $reference_allele = 'A';
my $alternative_allele = 'G';

# Fetch variation by location
my $variation = $va->fetch_by_location(
    $chromosome,
    $position,
    $position,
    $reference_allele,
    $alternative_allele
);

# Check if variation was found
if ($variation) {
    # Iterate over alleles and fetch frequency
    foreach my $allele (@{$variation->get_all_Alleles}) {
        my $allele_frequency = $allele->frequency;
        print "Allele Frequency: $allele_frequency\n";
    }
} else {
    print "Variant not found at specified location.\n";
}
