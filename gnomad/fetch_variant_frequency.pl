use strict;
use warnings;
use lib '/workspaces/Clinical-Genomics-Variant/ensembl-vep/ensembl/modules';
use Bio::EnsEMBL::Registry;

# Connect to Ensembl API
my $registry = 'Bio::EnsEMBL::Registry';
eval {
    $registry->load_registry_from_db(
        -host => 'ensembldb.ensembl.org',  # Ensembl public database host
        -user => 'anonymous',              # Anonymous access is common for public databases
    );
};
if ($@) {
    die "Error connecting to Ensembl database: $@\n";
}

# Get VariationAdaptor
my $va = $registry->get_adaptor('human', 'variation', 'variation');

# Check if adaptor was successfully retrieved
if (!defined $va) {
    die "Error: Unable to get variation adaptor.\n";
}

# Define variant details (replace with your variant's chromosome, position, ref, alt)
my $chromosome = '1';
my $position = 1000000;
my $reference_allele = 'A';
my $alternative_allele = 'G';

# Fetch variation by location
my $variation = $va->fetch_by_location(
    $chromosome,
    $position,
    $reference_allele,
    $alternative_allele
);

# Check if variation was found
if (defined $variation) {
    # Fetch all population-genotype frequencies
    my $allele = $variation->get_all_Alleles();
    foreach my $pf (@{$variation->get_all_PopulationGenotypeFeatures}) {
        my $population = $pf->population->name;
        my $frequency = $pf->frequency;
        print "Population: $population, Frequency: $frequency\n";
    }
} else {
    print "Variant not found at specified location.\n";
}
