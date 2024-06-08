# SAM File Short Indel Identifier

This Python script processes a SAM file to identify short indels (insertions and deletions) in the alignment data. The script identifies indels with lengths less than 3 and outputs the results to a text file, including the position of each indel.

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
{read_id}: [{operation}, {length}, {position}]
```

For example:

```
NR_039983.1: [('D', 1, 143)]
NM_024796.1: [('I', 1, 112)]
NM_005101.2: [('D', 1, 350)]
```


