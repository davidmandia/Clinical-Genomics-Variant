import os
import json
from collections import defaultdict

class VariantAnalyzer:
    def __init__(self, directory, output_file):
        self.directory = directory
        self.output_file = output_file
        self.variant_files = defaultdict(set)
        self.dict_name = defaultdict(list) 

    def _load_json(self, file_path):
        with open(file_path, 'r') as file:
            try:
                data = json.load(file)
                return data
            except json.JSONDecodeError:
                print(f"Error decoding JSON in file {file_path}")
                return {}

    def analyze(self):
        """
        Analyze the variants in all JSON files in the directory and count the number of files where each variant appears.
        """
        file_names = [f for f in os.listdir(self.directory) if f.endswith('.json')]
        for file_name in file_names:
            file_path = os.path.join(self.directory, file_name)
            data = self._load_json(file_path)           
            for variant in data.keys():
                self.variant_files[variant].add(file_name)
    
    def output(self):
        for variant, files in self.variant_files.items():
            file_list = tuple(sorted(files))
            self.dict_name[file_list].append(variant)

    def write_variant_counts_to_file(self):
        """
        Write all variants and their statistics to an output file.
        """
        with open(self.output_file, 'w') as file:
            file.write("file list,num of all the variants")
            file.write('\n')
            for files,variants in self.dict_name.items():
                count = len(variants)
                file.write(f'{files}, {count}\n')
                file.write('\n')

if __name__ == "__main__":
    directory_GRCh37 = 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh37'
    output_file_GRCh37 = 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh37_variant_counts.txt'
    directory_GRCh38 = 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh38'
    output_file_GRCh38 = 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh38_variant_counts.txt'
    directory_GRCh37_change = 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/select_variant_changing_UTR/GRCh37'
    output_file_GRCh37_change = 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/select_variant_changing_UTR/GRCh37_variant_counts.txt'
    directory_GRCh38_change = 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/select_variant_changing_UTR/GRCh38'
    output_file_GRCh38_change = 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/select_variant_changing_UTR/GRCh38_variant_counts.txt'
    analyzer_GRCh37 = VariantAnalyzer(directory_GRCh37, output_file_GRCh37)
    analyzer_GRCh37.analyze()
    analyzer_GRCh37.output()   
    analyzer_GRCh37.write_variant_counts_to_file()
    analyzer_GRCh38 = VariantAnalyzer(directory_GRCh38, output_file_GRCh38)
    analyzer_GRCh38.analyze()
    analyzer_GRCh38.output()   
    analyzer_GRCh38.write_variant_counts_to_file()
    analyzer_GRCh37_change = VariantAnalyzer(directory_GRCh37_change, output_file_GRCh37_change)
    analyzer_GRCh37_change.analyze()
    analyzer_GRCh37_change.output()   
    analyzer_GRCh37_change.write_variant_counts_to_file()
    analyzer_GRCh38_change = VariantAnalyzer(directory_GRCh38_change, output_file_GRCh38_change)
    analyzer_GRCh38_change.analyze()
    analyzer_GRCh38_change.output()   
    analyzer_GRCh38_change.write_variant_counts_to_file()
