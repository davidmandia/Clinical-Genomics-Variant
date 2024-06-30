# SAM File Short Indel Identifier

## Important:
I will update the README shortly, but I have changed the script and the arguments. In the meanwhile, you can check the requirement fold. The programm.sh file will run all scripts for all file including all necessary packages, input files download, and all the scripts for all input files

## Important

This Python script processes a SAM file to identify short indels (insertions and deletions) in the alignment data. The script identifies indels with lengths less than 3 and outputs the results to a text file, including the position of each indel.


If Mismatches are present, they will be present in the mismatch text file in output

## Requirements

- Python 3.6 or higher
- argparse module (usually included with Python)
- re module (usually included with Python)
- os module (usually included with Python)
- requests module (install it using `pip install requests`)

Additionally, you can download reference genome SAM files (GRCh37 and GRCh38) by running the `download_SAM_fromS3.py` Python script provided in the repository.

## Usage

1. Clone the repository or download the script.
2. Ensure you have a SAM file to process.
3. Run the script from the command line.

### Command Line Usage

```sh
python main.py path/to/your/samfile.sam
```

### Example

```sh
python main.py sample.sam
```

## Output

The script will output a VCF text file containing the identified indels in the SAM file, along with their positions. Additional info about the variant is present on the INFO section 

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	136889	.	NN	A	99	PASS	DP=100;LEN= 1 ;TYPE=DEL;TRANSCRIPT= NR_039983.1 ;TRANSCRIPT_POS= 45 ;CIGAR: 3=1X11=1X2=1X16=1X9=1D6=3X9=1X2=1X2785=93N58=227N492= ;GENOME_REF= NC_000001.11

```


### Example Output
#### Indels
```sh
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	136889	.	GG	G	99	PASS	DP=100;LEN= 1 ;TYPE=DEL;TRANSCRIPT= NR_039983.1 ;TRANSCRIPT_POS= 45 ;CIGAR: 3=1X11=1X2=1X16=1X9=1D6=3X9=1X2=1X2785=93N58=227N492= ;GENOME_REF= NC_000001.11
1	826578	.	C	CT	99	PASS	DP=100;LEN= 1 ;TYPE=INS;TRANSCRIPT= NM_024796.1 ;TRANSCRIPT_POS= 405 ;CIGAR: 33S146=1X19=1X47=1X10=1X146=1I315=1X315=1X2=1X8=1X30=1X254= ;GENOME_REF= NC_000001.11
1	1014534	.	AG	A	99	PASS	DP=100;LEN= 1 ;TYPE=DEL;TRANSCRIPT= NM_005101.2 ;TRANSCRIPT_POS= 628 ;CIGAR: 78=407N104=1X92=1X92=1X259=1D5= ;GENOME_REF= NC_000001.11
1	1014534	.	AG	A	99	PASS	DP=100;LEN= 1 ;TYPE=DEL;TRANSCRIPT= NM_005101.1 ;TRANSCRIPT_POS= 628 ;CIGAR: 78=407N104=1X92=1X92=1X259=1D5=1S ;GENOME_REF= NC_000001.11
```

```

## Converting BAM to SAM
You can convert a BAM file to SAM using the provided Bash script bam_to_sam_converter.sh. This script takes two arguments: the input BAM file and the output SAM file. Make sure you have samtools installed on your system before running this script.

### Requirements for conversion from Bam to Sam

``` sh
sudo apt-get update
sudo apt-get install samtools
```
Before using the converter you have to run to make it executable 

``` sh
chmod +x bam_to_sam_converter.sh
```
### Usage of converter

``` sh
bash convert_bam_to_sam.sh input.bam output.sam
```


