import unittest
from unittest.mock import patch, mock_open, MagicMock
import os
import json
from collections import defaultdict
from calculate_nums_in_different_jsonfiles import VariantAnalyzer 

class TestVariantAnalyzer(unittest.TestCase):

    @patch('os.listdir')
    @patch('builtins.open', new_callable=mock_open)
    @patch('json.load')
    def test_analyze(self, mock_json_load, mock_open_file, mock_listdir):
        mock_listdir.return_value = ['file1.json', 'file2.json']
        mock_json_load.side_effect = [
            {'variant1': 'data1', 'variant2': 'data2'},
            {'variant1': 'data3', 'variant3': 'data4'}
        ]

        analyzer = VariantAnalyzer('/_directory', '/_output.txt')
        analyzer.analyze()

        expected_variant_files = {
            'variant1': {'file1.json', 'file2.json'},
            'variant2': {'file1.json'},
            'variant3': {'file2.json'}
        }

        self.assertEqual(analyzer.variant_files, expected_variant_files)

    def test_output(self):
        analyzer = VariantAnalyzer('/_directory', '/_output.txt')
        analyzer.variant_files = {
            'variant1': {'file1.json', 'file2.json'},
            'variant2': {'file1.json'},
            'variant3': {'file2.json'}
        }

        analyzer.output()

        expected_dict_name = {
            ('file1.json', 'file2.json'): ['variant1'],
            ('file1.json',): ['variant2'],
            ('file2.json',): ['variant3']
        }

        self.assertEqual(analyzer.dict_name, expected_dict_name)

    @patch('builtins.open', new_callable=mock_open)
    def test_write_variant_counts_to_file(self, mock_open_file):
        analyzer = VariantAnalyzer('/_directory', '/_output.txt')
        analyzer.dict_name = {
            ('file1.json', 'file2.json'): ['variant1'],
            ('file1.json',): ['variant2'],
            ('file2.json',): ['variant3']
        }

        analyzer.write_variant_counts_to_file()

        mock_open_file.assert_called_once_with('/_output.txt', 'w')
        handle = mock_open_file()
        handle.write.assert_any_call("file list,num of all the variants")
        handle.write.assert_any_call("\n")
        handle.write.assert_any_call("('file1.json', 'file2.json'), 1\n")
        handle.write.assert_any_call("\n")
        handle.write.assert_any_call("('file1.json',), 1\n")
        handle.write.assert_any_call("\n")
        handle.write.assert_any_call("('file2.json',), 1\n")
        handle.write.assert_any_call("\n")

if __name__ == '__main__':
    unittest.main()
