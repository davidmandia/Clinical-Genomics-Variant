import unittest
from unittest.mock import patch, mock_open,call
import json
import io
import requests
from InvertVCFtoValidator import process_vcf  

class TestProcessVCF(unittest.TestCase):
    
    def setUp(self):
        self.vcf_file = 'testvalidator/test.vcf'
        self.lovd_file = 'testvalidator/lovd_data.json'
        self.source = 'GRCh37'
        self.missingfile='testvalidator/missing.txt'
        self.processor = process_vcf(self.vcf_file, self.lovd_file, self.source,self.missingfile)

    #here we don't really open the file,just use patch to mock the content got from the file    
    @patch("builtins.open", new_callable=mock_open, read_data="""\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t1570922\t.\tT\tTG\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_001291345.2;TRANSCRIPT_POS=335;GENOME_REF=NC_000001.10
1\t1570922\t.\tT\tTG\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_001787.3;TRANSCRIPT_POS=335;GENOME_REF=NC_000001.10
1\t1570922\t.\tT\tTG\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_033486.3;TRANSCRIPT_POS=335;GENOME_REF=NC_000001.10
1\t1570922\t.\tT\tTG\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_033487.3;TRANSCRIPT_POS=335;GENOME_REF=NC_000001.10
1\t1570922\t.\tT\tTG\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_033489.3;TRANSCRIPT_POS=335;GENOME_REF=NC_000001.10
1\t1570922\t.\tT\tTG\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_033490.3;TRANSCRIPT_POS=335;GENOME_REF=NC_000001.10
1\t1594199\t.\tC\tCT\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_001110781.3;TRANSCRIPT_POS=1261;GENOME_REF=NC_000001.10
1\t1594199\t.\tC\tCT\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_001290264.2;TRANSCRIPT_POS=1261;GENOME_REF=NC_000001.10
1\t7442597\t.\tT\tTG\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NR_146199.1;TRANSCRIPT_POS=50;GENOME_REF=NC_000001.10
1\t12570607\t.\tG\tGTG\t99\tPASS\tDP=100;LEN=2;TYPE=INS;TRANSCRIPT=NM_015378.4;TRANSCRIPT_POS=14863;GENOME_REF=NC_000001.10
1\t12570607\t.\tG\tGTG\t99\tPASS\tDP=100;LEN=2;TYPE=INS;TRANSCRIPT=NM_018156.4;TRANSCRIPT_POS=14788;GENOME_REF=NC_000001.10
""")#read data are used as the content of vcf file
    def test_read_vcf_file(self, mock_file):
        with open(self.vcf_file, 'r') as vcf:
            for line in vcf:
                print(f"Line: {line.strip()}")
        variants = self.processor.read_vcf_file()
        self.assertEqual(variants[0]['chrom'], '1')#compare the results
        self.assertEqual(variants[0]['pos'], '1570922')
        self.assertEqual(variants[0]['ref'], 'T')
        self.assertEqual(variants[0]['alt'], 'TG')
        self.assertEqual(variants[0]['qual'], '99')
        self.assertEqual(variants[0]['filter'], 'PASS')
        self.assertEqual(variants[0]['info'], 'DP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_001291345.2;TRANSCRIPT_POS=335;GENOME_REF=NC_000001.10')

    @patch("builtins.open", new_callable=mock_open, read_data="""\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t1570922\t.\tT\tTG\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_001291345.2;TRANSCRIPT_POS=335;GENOME_REF=NC_000001.10
1\t1570922\t.\tT\tTG\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_001787.3;TRANSCRIPT_POS=335;GENOME_REF=NC_000001.10
1\t1570922\t.\tT\tTG\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_033486.3;TRANSCRIPT_POS=335;GENOME_REF=NC_000001.10
1\t1570922\t.\tT\tTG\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_033487.3;TRANSCRIPT_POS=335;GENOME_REF=NC_000001.10
1\t1570922\t.\tT\tTG\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_033489.3;TRANSCRIPT_POS=335;GENOME_REF=NC_000001.10
1\t1570922\t.\tT\tTG\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_033490.3;TRANSCRIPT_POS=335;GENOME_REF=NC_000001.10
1\t1594199\t.\tC\tCT\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_001110781.3;TRANSCRIPT_POS=1261;GENOME_REF=NC_000001.10
1\t1594199\t.\tC\tCT\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NM_001290264.2;TRANSCRIPT_POS=1261;GENOME_REF=NC_000001.10
1\t7442597\t.\tT\tTG\t99\tPASS\tDP=100;LEN=1;TYPE=INS;TRANSCRIPT=NR_146199.1;TRANSCRIPT_POS=50;GENOME_REF=NC_000001.10
1\t12570607\t.\tG\tGTG\t99\tPASS\tDP=100;LEN=2;TYPE=INS;TRANSCRIPT=NM_015378.4;TRANSCRIPT_POS=14863;GENOME_REF=NC_000001.10
1\t12570607\t.\tG\tGTG\t99\tPASS\tDP=100;LEN=2;TYPE=INS;TRANSCRIPT=NM_018156.4;TRANSCRIPT_POS=14788;GENOME_REF=NC_000001.10
""")
    def test_get_vcf_calls(self, mock_file):
        vcf_calls = self.processor.get_vcf_calls()
        expected_calls = {'GRCh37:1:1570922:T:TG', 'GRCh37:1:12570607:G:GTG', 'GRCh37:1:1594199:C:CT', 'GRCh37:1:7442597:T:TG'}
        # Use set comparison
        self.assertEqual(set(vcf_calls), expected_calls)

    @patch('sys.stdout', new_callable=io.StringIO)       
    @patch('requests.get')#mock response got from api
    def test_query_API_lovd(self, mock_get,mock_stdout):
        mock_response = {
        "1:1570922:T:TG": {
            "1:1570922:T:TG": {
            "g_hgvs": "NC_000001.10:g.1570924dup",
            "genomic_variant_error": None,
            "hgvs_t_and_p": {
                "NM_001291345.2": {
                "gap_statement": "NM_001291345.2 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 237 extra bases between c.225_463 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_001291345.2 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_001278274.1:p.?",
                "p_hgvs_tlc": "NP_001278274.1:p.?",
                "select_status": {},
                "t_hgvs": "NM_001291345.2:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                },
                "NM_001787.3": {
                "gap_statement": "NM_001787.3 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 261 extra bases between c.231_493 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_001787.3 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_001778.2:p.?",
                "p_hgvs_tlc": "NP_001778.2:p.?",
                "select_status": {},
                "t_hgvs": "NM_001787.3:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                },
                "NM_033486.3": {
                "gap_statement": "NM_033486.3 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 261 extra bases between c.231_493 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_033486.3 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_277021.2:p.?",
                "p_hgvs_tlc": "NP_277021.2:p.?",
                "select_status": {
                    "mane_select": None,
                    "refseq_select": None
                },
                "t_hgvs": "NM_033486.3:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                },
                "NM_033487.3": {
                "gap_statement": "NM_033487.3 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 259 extra bases between c.-276_-16 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_033487.3 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_277022.1:p.?",
                "p_hgvs_tlc": "NP_277022.1:p.?",
                "select_status": {},
                "t_hgvs": "NM_033487.3:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                },
                "NM_033489.3": {
                "gap_statement": "NM_033489.3 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 261 extra bases between c.129_391 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_033489.3 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_277024.2:p.?",
                "p_hgvs_tlc": "NP_277024.2:p.?",
                "select_status": {},
                "t_hgvs": "NM_033489.3:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                },
                "NM_033490.3": {
                "gap_statement": "NM_033490.3 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 260 extra bases between c.-282_-21 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_033490.3 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_277025.1:p.?",
                "p_hgvs_tlc": "NP_277025.1:p.?",
                "select_status": {},
                "t_hgvs": "NM_033490.3:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                }
            },
            "p_vcf": "1-1570922-T-TG",
            "selected_build": "GRCh37"
            },
            "errors": [],
            "flag": None
        },
        "metadata": {
            "variantformatter_version": "2.2.1.dev70+gbe03152",
            "variantvalidator_hgvs_version": "2.2.0",
            "variantvalidator_version": "2.2.1.dev689+g995c0b0",
            "vvdb_version": "vvdb_2024_5",
            "vvseqrepo_db": "VV_SR_2024_04/master",
            "vvta_version": "vvta_2024_01"
        }
        }
        mock_get.return_value.status_code = 200
        mock_get.return_value.json.return_value = mock_response
        
        result = self.processor.query_API_lovd('GRCh37', '1:1570922:T:TG', transcript_model="all", select_transcripts="raw", checkonly="False", liftover="False")
        self.assertEqual(result, mock_response)
        mock_get.side_effect = requests.exceptions.HTTPError(
            "500 Server Error: INTERNAL SERVER ERROR for url: https://rest.variantvalidator.org/LOVD/lovd/GRCh37/NW_003315937.1:12760:GGA:G/all/raw/False/False?content-type=application/json"
        )
        self.processor.query_API_lovd('GRCh37', 'NW_003315937.1:12760:GGA:G', transcript_model="all", select_transcripts="raw", checkonly="False", liftover="False")
        self.assertEqual(mock_stdout.getvalue().strip(), "HTTP error occurred: 500 Server Error: INTERNAL SERVER ERROR for url: https://rest.variantvalidator.org/LOVD/lovd/GRCh37/NW_003315937.1:12760:GGA:G/all/raw/False/False?content-type=application/json")
       
    def test_format_data(self):
        vcf_json = {
        "1:1570922:T:TG": {
            "1:1570922:T:TG": {
            "g_hgvs": "NC_000001.10:g.1570924dup",
            "genomic_variant_error": None,
            "hgvs_t_and_p": {
                "NM_001291345.2": {
                "gap_statement": "NM_001291345.2 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 237 extra bases between c.225_463 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_001291345.2 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_001278274.1:p.?",
                "p_hgvs_tlc": "NP_001278274.1:p.?",
                "select_status": {},
                "t_hgvs": "NM_001291345.2:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                },
                "NM_001787.3": {
                "gap_statement": "NM_001787.3 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 261 extra bases between c.231_493 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_001787.3 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_001778.2:p.?",
                "p_hgvs_tlc": "NP_001778.2:p.?",
                "select_status": {},
                "t_hgvs": "NM_001787.3:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                },
                "NM_033486.3": {
                "gap_statement": "NM_033486.3 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 261 extra bases between c.231_493 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_033486.3 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_277021.2:p.?",
                "p_hgvs_tlc": "NP_277021.2:p.?",
                "select_status": {
                    "mane_select": None,
                    "refseq_select": None
                },
                "t_hgvs": "NM_033486.3:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                },
                "NM_033487.3": {
                "gap_statement": "NM_033487.3 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 259 extra bases between c.-276_-16 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_033487.3 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_277022.1:p.?",
                "p_hgvs_tlc": "NP_277022.1:p.?",
                "select_status": {},
                "t_hgvs": "NM_033487.3:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                },
                "NM_033489.3": {
                "gap_statement": "NM_033489.3 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 261 extra bases between c.129_391 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_033489.3 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_277024.2:p.?",
                "p_hgvs_tlc": "NP_277024.2:p.?",
                "select_status": {},
                "t_hgvs": "NM_033489.3:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                },
                "NM_033490.3": {
                "gap_statement": "NM_033490.3 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 260 extra bases between c.-282_-21 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_033490.3 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_277025.1:p.?",
                "p_hgvs_tlc": "NP_277025.1:p.?",
                "select_status": {},
                "t_hgvs": "NM_033490.3:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                }
            },
            "p_vcf": "1-1570922-T-TG",
            "selected_build": "GRCh37"
            },
            "errors": [],
            "flag": None
        },
        "metadata": {
            "variantformatter_version": "2.2.1.dev70+gbe03152",
            "variantvalidator_hgvs_version": "2.2.0",
            "variantvalidator_version": "2.2.1.dev689+g995c0b0",
            "vvdb_version": "vvdb_2024_5",
            "vvseqrepo_db": "VV_SR_2024_04/master",
            "vvta_version": "vvta_2024_01"
        }
        }
        expected_output = {
        "<1:1570922:T:TG>": {
            "refseq": [
                "NM_001291345.2:c.*203_*205=",
                "NM_001787.3:c.*203_*205=",
                "NM_033486.3:c.*203_*205=",
                "NM_033487.3:c.*203_*205=",
                "NM_033489.3:c.*203_*205=",
                "NM_033490.3:c.*203_*205="
            ],
            "ensembl": []
        }
    }
        vcf_json2={
        "error": "local variable 'result' referenced before assignment"
        }
        expected_output2={
            "<NW_003315937.1:12760:GGA:G>": {
                "refseq": [None],
                "ensembl": [None]
            }
        }
        self.assertEqual(self.processor.format_data('1:1570922:T:TG', vcf_json), expected_output)
        self.assertEqual(self.processor.format_data('NW_003315937.1:12760:GGA:G', vcf_json2), expected_output2)
    
    @patch('builtins.open', new_callable=mock_open)
    def test_write_to_lovdjson(self, mock_file):
        data = {'test': 'data'}
        self.processor.write_to_lovdjson(data)
        expected_calls = [
        call('{'),
        call('\n    '),
        call('"test"'),
        call(': '),
        call('"data"'),
        call('\n'),
        call('}'),
        call('\n')
        ]   
        mock_file().write.assert_has_calls(expected_calls, any_order=False) #compare the results in the file
    
    @patch('builtins.open', new_callable=mock_open)
    def test_clear_write_file(self, mock_file):
        data = {'test': 'data'}
        self.processor.clear_write_file(data)
        expected_calls = [
        call('{'),
        call('\n    '),
        call('"test"'),
        call(': '),
        call('"data"'),
        call('\n'),
        call('}'),
        call('\n')
        ]   
        mock_file().write.assert_has_calls(expected_calls, any_order=False)   
    @patch('builtins.open', new_callable=mock_open)
    def test_write_missing_seq(self, mock_file):
        data = ['name1','name2']
        self.processor.write_missing_seq(data)
        expected_calls = [
        call('name1\n'),
        call('name2\n')
        ]   
        mock_file().write.assert_has_calls(expected_calls, any_order=False)

    @patch('requests.get')
    def test_variant_validator_lovd_api_call(self, mock_get):
        mock_response = {
        "1:1570922:T:TG": {
            "1:1570922:T:TG": {
            "g_hgvs": "NC_000001.10:g.1570924dup",
            "genomic_variant_error": None,
            "hgvs_t_and_p": {
                "NM_001291345.2": {
                "gap_statement": "NM_001291345.2 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 237 extra bases between c.225_463 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_001291345.2 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_001278274.1:p.?",
                "p_hgvs_tlc": "NP_001278274.1:p.?",
                "select_status": {},
                "t_hgvs": "NM_001291345.2:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                },
                "NM_001787.3": {
                "gap_statement": "NM_001787.3 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 261 extra bases between c.231_493 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_001787.3 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_001778.2:p.?",
                "p_hgvs_tlc": "NP_001778.2:p.?",
                "select_status": {},
                "t_hgvs": "NM_001787.3:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                },
                "NM_033486.3": {
                "gap_statement": "NM_033486.3 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 261 extra bases between c.231_493 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_033486.3 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_277021.2:p.?",
                "p_hgvs_tlc": "NP_277021.2:p.?",
                "select_status": {
                    "mane_select": None,
                    "refseq_select": None
                },
                "t_hgvs": "NM_033486.3:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                },
                "NM_033487.3": {
                "gap_statement": "NM_033487.3 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 259 extra bases between c.-276_-16 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_033487.3 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_277022.1:p.?",
                "p_hgvs_tlc": "NP_277022.1:p.?",
                "select_status": {},
                "t_hgvs": "NM_033487.3:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                },
                "NM_033489.3": {
                "gap_statement": "NM_033489.3 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 261 extra bases between c.129_391 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_033489.3 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_277024.2:p.?",
                "p_hgvs_tlc": "NP_277024.2:p.?",
                "select_status": {},
                "t_hgvs": "NM_033489.3:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                },
                "NM_033490.3": {
                "gap_statement": "NM_033490.3 contains 1 extra bases between c.*204_*205, and 13 extra bases between c.*210_*211, and 260 extra bases between c.-282_-21 than NC_000001.10",
                "gapped_alignment_warning": "Submitted description does not represent a true variant because it is an artefact of aligning NM_033490.3 with NC_000001.10 (genome build GRCh37)",
                "gene_info": {
                    "hgnc_id": "HGNC:1729",
                    "symbol": "CDK11B"
                },
                "p_hgvs_slc": "NP_277025.1:p.?",
                "p_hgvs_tlc": "NP_277025.1:p.?",
                "select_status": {},
                "t_hgvs": "NM_033490.3:c.*203_*205=",
                "transcript_variant_error": None,
                "transcript_version_warning": None
                }
            },
            "p_vcf": "1-1570922-T-TG",
            "selected_build": "GRCh37"
            },
            "errors": [],
            "flag": None
        },
        "metadata": {
            "variantformatter_version": "2.2.1.dev70+gbe03152",
            "variantvalidator_hgvs_version": "2.2.0",
            "variantvalidator_version": "2.2.1.dev689+g995c0b0",
            "vvdb_version": "vvdb_2024_5",
            "vvseqrepo_db": "VV_SR_2024_04/master",
            "vvta_version": "vvta_2024_01"
        }
        }
        expected_output = {
        "<1:1570922:T:TG>": {
            "refseq": [
                "NM_001291345.2:c.*203_*205=",
                "NM_001787.3:c.*203_*205=",
                "NM_033486.3:c.*203_*205=",
                "NM_033487.3:c.*203_*205=",
                "NM_033489.3:c.*203_*205=",
                "NM_033490.3:c.*203_*205="
            ],
            "ensembl": []
        }
        }
        mock_get.return_value.json.return_value = mock_response
        with patch('json.dumps') as mock_dumps:
            mock_dumps.return_value = json.dumps(mock_response)
            self.processor.clear_write_file = lambda: None  # Mock method to avoid file operations
            self.processor.write_to_lovdjson = lambda data: None  # Mock method to avoid file operations
            self.processor.write_missing_seq = lambda missing: None  # Mock method to avoid file operations
            data = self.processor.query_API_lovd('GRCh37', '1:1570922:T:TG', transcript_model="all", select_transcripts="raw", checkonly="False", liftover="False")
            formatted_data = self.processor.format_data('1:1570922:T:TG',data)
            self.processor.write_to_lovdjson(formatted_data)
            mock_get.assert_called_once_with('https://rest.variantvalidator.org/LOVD/lovd/GRCh37/1:1570922:T:TG/all/raw/False/False?content-type=application/json')
            self.assertEqual(expected_output, formatted_data)

if __name__ == '__main__':
    unittest.main()
