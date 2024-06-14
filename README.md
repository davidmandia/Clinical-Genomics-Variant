# SAM File Short Indel Identifier

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
python identify_indels.py path/to/your/samfile.sam
```

### Example

```sh
python identify_indels.py example.sam
```

## Output

The script will output a text file containing the identified indels in the SAM file, along with their positions. Each indel will be represented by a line with the following format:

```
{read_id}: [{genome_reference}, ({operation}, {length}, {position})]
```


### Example Output
#### Indels
```sh
NM_001291345.2: [('NC_000001.10', (('I', 1, 335),))]
NM_001290264.2: ['NC_000001.10', [('I', 1, 1261)]]
NR_146199.1: ['NC_000001.10', [('I', 1, 50)]]
NM_015378.4: ['NC_000001.10', [('I', 2, 280521)]]
`
```
#### Mismatches
```sh
NM_001291345.2: ['NC_000001.10', [('X', 1, 335)]]
NM_001290264.2: ['NC_000001.10', [('X', 1, 1261)]]
NR_146199.1: ['NC_000001.10', [('X', 1, 50)]]
NM_015378.4: ['NC_000001.10', [('X', 2, 280521)]]
```

