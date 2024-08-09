import json
import re

class JudgementVariantList:
    def __init__(self, final_path,judgement_file_path_five,judgement_file_path_three,judgement_file_path_refseq, judgement_file_path_ensembl,judgement_file_path_refseq_ensembl,judgement_file_path_no_gap, lovd_file, source):
        self.final_path=final_path
        self.judgement_file_path_five = judgement_file_path_five
        self.judgement_file_path_three=judgement_file_path_three
        self.judgement_file_path_refseq=judgement_file_path_refseq
        self.judgement_file_path_ensembl=judgement_file_path_ensembl
        self.judgement_file_path_refseq_ensembl=judgement_file_path_refseq_ensembl
        self.judgement_file_path_no_gap=judgement_file_path_no_gap
        self.source = source
        self.lovd_file = lovd_file
    
    def read_list(self):
        data = {}
        with open(self.lovd_file, 'r') as file:
            content = file.read()
            # Split by row
            lines = content.split('\n')           
            # Concatenate each line and parse it into a JSON object
            current_object = ''
            for line in lines:
                line = line.strip()
                if line:
                    current_object += line
                    # If there is a complete JSON object in the current object
                    try:
                        parsed_object = json.loads(current_object)
                        data.update(parsed_object)
                        current_object = ''
                    except json.JSONDecodeError:
                        continue
        return data
    
    def is_five_prime_utr(self, variant):
        return re.search(r'c\.-', variant)# 'c.-' is not a valid regex; changed to match 5' UTR
    
    def is_three_prime_utr(self, variant):
        return re.search(r'c\.\*', variant) # 'c.*' is not a valid regex; changed to match 3' UTR

    def has_gap(self, variant):
        return re.search(r'=$|ins[GATC]+$', variant)

    def process_list(self):
        data_dic = self.read_list()
        categories = {
            '5_prime_UTR': [],
            '3_prime_UTR': [],
            'gap_in_refseq_and_ensembl': [],
            'gap_in_refseq_not_in_ensembl': [],
            'gap_not_in_refseq_but_in_ensembl': [],
            'no_gap_found': []
            }
        for variant_key,ref in data_dic.items():
            refseq_list = ref.get('refseq', [])
            ensembl_list = ref.get('ensembl', [])
            # check whether has gaps
            refseq_gaps = any(self.has_gap(v) for v in refseq_list if v)
            ensembl_gaps = any(self.has_gap(v) for v in ensembl_list if v)
            
            if refseq_gaps and ensembl_gaps:
                categories['gap_in_refseq_and_ensembl'].append(variant_key)
            elif refseq_gaps and not ensembl_gaps:
                categories['gap_in_refseq_not_in_ensembl'].append(variant_key)
            elif not refseq_gaps and ensembl_gaps:
                categories['gap_not_in_refseq_but_in_ensembl'].append(variant_key)
            else:
                categories['no_gap_found'].append(variant_key)
            
            # check 5' and 3' UTR
            if any(self.is_five_prime_utr(v) for v in refseq_list + ensembl_list if v):
                categories['5_prime_UTR'].append(variant_key)
            if any(self.is_three_prime_utr(v) for v in refseq_list + ensembl_list if v):
                categories['3_prime_UTR'].append(variant_key)
        return categories
    
    def read_pre(self,dict,dict_categories,output_json,num):
        key_index=list(dict_categories.keys())
        for line in dict_categories[key_index[num]]:
            if not line.endswith('\n'):
                dict[line]=output_json[line]
        return dict        

    def write_to_file(self):
        content=self.process_list()
        output=self.read_list()
        with open(self.final_path,'w') as files:
            json.dump(output,files,indent=4)
        with open(self.judgement_file_path_five,'w') as file1:
            dict1={}
            self.read_pre(dict1,content,output,0)
            json.dump(dict1,file1,indent=4)
        with open(self.judgement_file_path_three,'w') as file2:
            dict2={}
            self.read_pre(dict2,content,output,1)
            json.dump(dict2,file2,indent=4)
        with open(self.judgement_file_path_refseq,'w') as file3:
            dict3={}
            self.read_pre(dict3,content,output,2)
            json.dump(dict3,file3,indent=4)
        with open(self.judgement_file_path_ensembl,'w') as file4:
            dict4={}
            self.read_pre(dict4,content,output,3)
            json.dump(dict4,file4,indent=4)
        with open(self.judgement_file_path_refseq_ensembl,'w') as file5:
            dict5={}
            self.read_pre(dict5,content,output,4)
            json.dump(dict5,file5,indent=4)
        with open(self.judgement_file_path_no_gap,'w') as file6:
            dict6={}
            self.read_pre(dict6,content,output,5)
            json.dump(dict6,file6,indent=4)
if __name__ == "__main__":
    lovd_json_path_GRCh37 = 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/GRCh37lovd.json'
    lovd_json_path_GRCh38 = 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/GRCh38lovd.json'
    final_path_GRch37='C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh37_lovd.json'
    final_path_GRch38='C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh38_lovd.json'
    judgement_file_path37_five = 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh37/GRCh37_five_prime_utr.json'
    judgement_file_path37_three= 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh37/GRCh37_three_prime_utr.json'
    judgement_file_path37_refseq= 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh37/GRCh37_gaps_in_refseq_not_in_ensembl.json'
    judgement_file_path37_ensembl= 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh37/GRCh37_gaps_not_in_refseq_in_ensembl.json'
    judgement_file_path37_refseq_ensembl= 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh37/GRCh37_gap_in_refseq_and_ensembl.json'
    judgement_file_path37_no_gap= 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh37/GRCh37_no_gap.json'
    judgement_file_path38_five = 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh38/GRCh38_five_prime_utr.json'
    judgement_file_path38_three= 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh38/GRCh38_three_prime_utr.json'
    judgement_file_path38_refseq= 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh38/GRCh38_gaps_in_refseq_not_in_ensembl.json'
    judgement_file_path38_ensembl= 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh38/GRCh38_gaps_not_in_refseq_in_ensembl.json'
    judgement_file_path38_refseq_ensembl= 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh38/GRCh38_gap_in_refseq_and_ensembl.json'
    judgement_file_path38_no_gap= 'C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/Select_variant/GRCh38/GRCh38_no_gap.json'
    Select_result_37 = JudgementVariantList(final_path_GRch37,judgement_file_path37_five,judgement_file_path37_three,judgement_file_path37_refseq, judgement_file_path37_ensembl,judgement_file_path37_refseq_ensembl,judgement_file_path37_no_gap,lovd_json_path_GRCh37, 'GRCh37')
    Select_result_38 = JudgementVariantList(final_path_GRch38,judgement_file_path38_five,judgement_file_path38_three,judgement_file_path38_refseq, judgement_file_path38_ensembl,judgement_file_path38_refseq_ensembl,judgement_file_path38_no_gap,lovd_json_path_GRCh38, 'GRCh38')
    result37=Select_result_37.write_to_file()
    result38=Select_result_38.write_to_file()
