### ETL Document for Clinical Genomics Variant Database

#### **1. Overview**

This ETL process describes the comprehensive steps taken to create a clinical genomics variant database. The final database contains the varaint information( chrom number, position, ref, alt) of indels varaints found in the reference transcritpme alignment.
It covers the extraction of data from transcriptome sequences, variant identification, data transformation, and integration with various datasets, and the final loading into a PostgreSQL database on AWS. The database includes genomic variant information, allele frequencies, gene annotations, and clinical relevance data.

#### **2. Data Sources**

1. **Transcriptome Sequences**: RefSeq transcriptme alignment to genome builds GRCh37 and GRCH38 
2. **Blast DB**: Blast NCBI command line tool to retrieve reference sequence for the VCF
3. **gnomAD Database**: Allele frequency data across various populations.
4. **GFF Files**: Annotation files for gene and exon structures.
5. **PanelApp**: Clinical relevance and disease association data.

#### **3. Extract Phase**

##### **3.1. Transcriptome Alignment and Indel Identification**

- **Tools Used**: BLAST, Ensembl VEP
- **Process**:
  - **Alignment**: Align raw transcriptome sequences to a reference genome using BLAST or other alignment tools.
  - **Variant Calling**: Identify indels (insertions and deletions). Verify transcriptomic position and genomic position.  Generate VCF (Variant Call Format) files as the output, containing variant information such as chromosome, position, reference, and alternate sequences.

##### **3.2. Extracting VCF Data**

- **Tools Used**: Python, Ensembl VEP
- **Process**:
  - Parse VCF files to extract essential variant fields including chromosome, position, reference allele, alternate allele and use them to retireve variant consequences from Ensembl VEP

##### **3.3. Extracting gnomAD Data**

- **Tools Used**: Python, Ensembl VEP
- **Process**:
  - Extract allele frequency data for each variant across different populations, which helps in assessing the commonality of the variant.
  - Populations include general population(af), nfe (Non-Finnish European), eas( East Asian), afr (African)

##### **3.4. Extracting GFF Data**

- **Tools Used**: Python, GFF parsing libraries
- **Process**:
  - Extract structural annotations such as gene IDs, gene symbols, and exon boundaries. This helps in mapping variants to specific genes and exons.

##### **3.5. Extracting Clinical Relevance from PanelApp**

- **Tools Used**: Python, PanelApp API or datasets
- **Process**:
  - Extract clinical annotations including disease associations and relevance of specific genes, providing insights into the potential clinical impact of the variants.

#### **4. Transform Phase**

##### **4.1. Data Cleaning**

- **Tools Used**: Python (Pandas)
- **Process**:
  - Standardize chromosome nomenclature and genomic positions.
  - Ensure uniform representation of alleles and remove any formatting inconsistencies.


##### **4.2. Data Integration**

- **Tools Used**: Python (Pandas, SQLAlchemy)
- **Process**:
  - Integrate VCF data with allele frequency information from gnomAD.
  - Annotate variants with gene and exon information based on GFF files.
  - Link clinical annotations from PanelApp to corresponding variants.

##### **4.3. Schema Alignment**

- **Tools Used**: Python
- **Process**:
  - Align the data to fit into the database schema, including mapping fields to the appropriate columns and converting data types.
  - Calculate additional attributes such as variant consequences, variant distribution across the genome, affected genes and transcripts.

#### **5. Load Phase**

##### **5.1. Database Preparation**

- **Tools Used**: PostgreSQL, AWS RDS
- **Process**:
  - Design the database schema, including tables for variants, allele frequencies, gene annotations, and clinical relevance.

##### **5.2. Data Loading**

- **Tools Used**: Python , PostgreSQL
- **Process**:
  - Load the cleaned and transformed data into the corresponding tables.

##### **5.3. Post-Load Validation**

- **Tools Used**: SQL, Python
- **Process**:
  - Conduct thorough data validation to check for completeness and accuracy.
  - Cross-reference with original datasets to ensure data consistency.

#### **6. Deployment**

- **Tools Used**: AWS Lambda, API Gateway, S3
- **Process**:
  - Develop a web interface for querying the database.
  - Set up API endpoints to allow external access and data retrieval.
  - Implement security measures and monitoring systems to maintain the database's integrity and performance.

#### **7. Maintenance and Updates**

- **Regular updates**: Continuously update the database with new genomic data and annotations.
- **Monitoring**: Use monitoring tools to ensure the database remains operational and to track data quality.

### Conclusion

This ETL document outlines the detailed process of building a comprehensive clinical genomics variant database, from initial variant identification in transcriptome data to the deployment of a searchable and accessible database. The system integrates multiple data sources, providing a valuable resource for clinical and research applications.