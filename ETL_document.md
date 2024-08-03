### ETL Document for Clinical Genomics Variant Database

### Complete Database 
A brief explanation of each column:

1. **chrom**
   - The chromosome number where the variant is located.

2. **pos**
   - The position of the variant on the chromosome.

3. **ref**
   - The reference allele at the variant position.

4. **alt**
   - The alternative allele that is different from the reference allele at the variant position.

5. **qual**
   - The quality score assigned to the variant, indicating the confidence in the variant call.

6. **filter**
   - Information about whether the variant passed certain quality control filters.

7. **genomic_ref**
   - Reference sequence identifier for the genomic reference used.

8. **operation**
   - The type of variant operation, such as insertion, deletion, or substitution.

9. **transcript_ref**
   - Reference sequence identifier for the transcript used.

10. **transcript_pos**
    - The position of the variant within the transcript.

11. **af**
    - Allele frequency of the variant in the gnomAD database.

12. **af_eas**
    - Allele frequency of the variant in the East Asian population.

13. **af_nfe**
    - Allele frequency of the variant in the Non-Finnish European population.

14. **af_fin**
    - Allele frequency of the variant in the Finnish population.

15. **af_amr**
    - Allele frequency of the variant in the Latino/Admixed American population.

16. **af_afr**
    - Allele frequency of the variant in the African/African American population.

17. **af_asj**
    - Allele frequency of the variant in the Ashkenazi Jewish population.

18. **af_oth**
    - Allele frequency of the variant in other or mixed populations.

19. **af_sas**
    - Allele frequency of the variant in the South Asian population.

20. **af_mid**
    - Allele frequency of the variant in the Middle Eastern population.

21. **af_ami**
    - Allele frequency of the variant in the Amish population.

22. **genes**
    - List of genes affected by the variant.

23. **consequences**
    - Predicted consequences of the variant on the gene(s) and transcript(s), such as missense, nonsense, or frameshift.

24. **existing_variant**
    - Indicates if the variant is already known or previously reported.

25. **gene_symbol**
    - The official symbol of the gene(s) affected by the variant.

26. **gene_id**
    - The unique identifier for the gene(s) affected by the variant.

27. **clinically_relevant**
    - Indicates whether the variant has clinical relevance.

28. **clinical_label**
    - The clinical significance or label assigned to the variant, such as pathogenic, benign, or uncertain significance.

29. **exon_id**
    - Identifier for the exon affected by the variant, if applicable.

#### **1. Overview**

This ETL process describes the comprehensive steps taken to create a clinical genomics variant database. The final database contains the varaint information( chrom number, position, ref, alt) of indels varaints found in the reference transcritpme alignment.
It covers the extraction of data from transcriptome sequences, variant identification, data transformation, and integration with various datasets, and the final loading into a PostgreSQL database on AWS. The database includes genomic variant information, allele frequencies, gene annotations, and clinical relevance data.

#### **2. Data Sources**

1. **Transcriptome Sequences**: RefSeq transcriptme alignment to genome builds GRCh37 and GRCH38 (Downloaded the bam.gz files)
     - Link to GRCH38 file  Index of https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/historical/GRCh38/GCF_000001405.40-RS_2023_03_historical/ 
     - Link To GRCh37 file  Index of https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/105.20201022/GCF_000001405.25_GRCh37.p13/RefSeq_transcripts_alignments/
2. **Blast DB**: Blast NCBI command line tool to retrieve reference sequence for the VCF
     - Link to tool instruction https://www.ncbi.nlm.nih.gov/books/NBK279690/  
3. **gnomAD Database**: Allele frequency data across various populations.
     - Link to Ensembl VEP APi for gnomad data and variant consequences https://rest.ensembl.org/ 
4. **GFF Files**: Annotation files for gene and exon structures.
5. **PanelApp**: Clinical relevance and disease association data.
     - PanelApp link for API https://panelapp.genomicsengland.co.uk/api/docs/

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
