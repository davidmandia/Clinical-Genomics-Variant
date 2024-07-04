# pip install requests
# pip install vobject

# pip install sqlite3

# ## For testing 

# pip install pytest

# ## Initiate the database 
# #sqlite3 variants.db


## VEP and GNOMAD

git clone https://github.com/Ensembl/ensembl-vep
cd ensembl-vep
sudo apt-get install libdbi-perl


# Install cpanm using cpan
cpan App::cpanminus

# Use cpanm to install required Perl modules
cpanm JSON Try::Tiny URI::Escape LWP::UserAgent HTML::Entities DBI DBD::mysql Archive::Extract Bio::Perl

# Navigate to the ensembl-vep directory (if not already there)
cd /ensembl-vep

# Clone necessary Ensembl API repositories
git clone https://github.com/Ensembl/ensembl.git
git clone https://github.com/Ensembl/ensembl-variation.git
git clone https://github.com/Ensembl/ensembl-funcgen.git
git clone https://github.com/Ensembl/ensembl-io.git
git clone https://github.com/Ensembl/ensembl-compara.git

# Set the PERL5LIB environment variable
export PERL5LIB=/ensembl-vep/ensembl/modules:/ensembl-vep/ensembl-variation/modules:/ensembl-vep/ensembl-funcgen/modules:/ensembl-vep/ensembl-io/modules:/ensembl-vep/ensembl-compara/modules:$PERL5LIB

cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)
eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)
source ~/.bashrc  # or ~/.bash_profile or ~/.zshrc
perl -MDBD::mysql -e 'print "DBD::mysql is installed\n";'


# Run the VEP installer script
perl INSTALL.pl

# Run VEP
./vep -i examples/homo_sapiens_GRCh38.vcf --cache

## HAIL  -- REQUIRES CLOUD COMPUTING

# sudo apt-get install -y openjdk-11-jre-headless g++ python3.9 python3-pip libopenblas-base liblapack3
# python -m pip install hail

# # Update package list and install dependencies
# sudo apt-get update && sudo apt-get install -y apt-transport-https ca-certificates gnupg

# # Add the Cloud SDK distribution URI as a package source
# echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | sudo tee -a /etc/apt/sources.list.d/google-cloud-sdk.list

# # Import the Google Cloud public key
# curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -

# # Update and install the Google Cloud SDK
# sudo apt-get update && sudo apt-get install -y google-cloud-sdk

cpan BioPerl
cpan DBI