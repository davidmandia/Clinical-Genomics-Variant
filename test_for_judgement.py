import unittest
import json
from unittest.mock import patch, mock_open
from Judgement_variant import JudgementVariantList

class TestJudgementVariantList(unittest.TestCase):
    
    def setUp(self):
        self.final_path = 'final_path.json'
        self.judgement_file_path_five = 'five_prime_utr.json'
        self.judgement_file_path_three = 'three_prime_utr.json'
        self.judgement_file_path_refseq = 'gaps_in_refseq_not_in_ensembl.json'
        self.judgement_file_path_ensembl = 'gaps_not_in_refseq_in_ensembl.json'
        self.judgement_file_path_refseq_ensembl = 'gap_in_refseq_and_ensembl.json'
        self.judgement_file_path_no_gap = 'no_gap.json'
        self.lovd_file = 'lovd.json'
        self.source = 'GRCh37'
        
        self.judgement = JudgementVariantList(
            self.final_path,
            self.judgement_file_path_five,
            self.judgement_file_path_three,
            self.judgement_file_path_refseq,
            self.judgement_file_path_ensembl,
            self.judgement_file_path_refseq_ensembl,
            self.judgement_file_path_no_gap,
            self.lovd_file,
            self.source
        )
        
    @patch("builtins.open", new_callable=mock_open, read_data='{"<8:21966700:GC:G>": {"refseq": ["NM_024815.3:c.108_110=", "NM_024815.4:c.110_112="], "ensembl": ["ENST00000309188.6:c.108_116=", "ENST00000521807.2:c.113del", "ENST00000522405.1:c.-121del"]}}')
    def test_read_list(self, mock_file):
        expected_output = {
            "<8:21966700:GC:G>": {
                "refseq": ["NM_024815.3:c.108_110=", "NM_024815.4:c.110_112="],
                "ensembl": ["ENST00000309188.6:c.108_116=", "ENST00000521807.2:c.113del", "ENST00000522405.1:c.-121del"]
            }
        }
        result = self.judgement.read_list()
        self.assertEqual(result, expected_output)
        mock_file.assert_called_once_with(self.lovd_file, 'r')
        
    def test_is_five_prime_utr(self):
        variant = "ENST00000522405.1:c.-121del"
        self.assertTrue(self.judgement.is_five_prime_utr(variant))
        
    def test_is_three_prime_utr(self):
        variant = "NM_024815.4:c.110_112="
        self.assertFalse(self.judgement.is_three_prime_utr(variant))
        
    def test_has_gap(self):
        variant1 = "NM_024815.3:c.108_110="
        variant2 = "ENST00000309188.6:c.108_116="
        variant="NM_024815.4:c.110_112="
        self.assertTrue(self.judgement.has_gap(variant1))
        self.assertTrue(self.judgement.has_gap(variant2))
        self.assertTrue(self.judgement.has_gap(variant))
        
    @patch("builtins.open", new_callable=mock_open, read_data='{"<8:21966700:GC:G>": {"refseq": ["NM_024815.3:c.108_110=", "NM_024815.4:c.110_112="], "ensembl": ["ENST00000309188.6:c.108_116=", "ENST00000521807.2:c.113del", "ENST00000522405.1:c.-121del"]}}')
    def test_process_list(self, mock_file):
        expected_output = {
            '5_prime_UTR': ['<8:21966700:GC:G>'],
            '3_prime_UTR': [],
            'gap_in_refseq_and_ensembl': ['<8:21966700:GC:G>'],
            'gap_in_refseq_not_in_ensembl': [],
            'gap_not_in_refseq_but_in_ensembl': [],
            'no_gap_found': []
        }
        result = self.judgement.process_list()
        self.assertEqual(result, expected_output)
        
    @patch("builtins.open", new_callable=mock_open)
    def test_write_to_file(self, mock_file):
        # Simulate process_list and read_list return values
        with patch.object(self.judgement, 'process_list', return_value={
            '5_prime_UTR': ['<8:21966700:GC:G>'],
            '3_prime_UTR': [],
            'gap_in_refseq_and_ensembl': ['<8:21966700:GC:G>'],
            'gap_in_refseq_not_in_ensembl': [],
            'gap_not_in_refseq_but_in_ensembl': [],
            'no_gap_found': []
        }):
            with patch.object(self.judgement, 'read_list', return_value={
                '<8:21966700:GC:G>': {
                    "refseq": ["NM_024815.3:c.108_110=", "NM_024815.4:c.110_112="],
                    "ensembl": ["ENST00000309188.6:c.108_116=", "ENST00000521807.2:c.113del", "ENST00000522405.1:c.-121del"]
                }
            }):
                self.judgement.write_to_file()
                mock_file.assert_any_call(self.final_path, 'w')
                mock_file.assert_any_call(self.judgement_file_path_five, 'w')
                mock_file.assert_any_call(self.judgement_file_path_three, 'w')
                mock_file.assert_any_call(self.judgement_file_path_refseq, 'w')
                mock_file.assert_any_call(self.judgement_file_path_ensembl, 'w')
                mock_file.assert_any_call(self.judgement_file_path_refseq_ensembl, 'w')
                mock_file.assert_any_call(self.judgement_file_path_no_gap, 'w')

if __name__ == '__main__':
    unittest.main()
