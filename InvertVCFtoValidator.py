import requests
import json
import time
class process_vcf:
    def __init__(self,vcf_file,lovd_file,source,missingvariantfile):
        self.vcf_file=vcf_file
        self.source=source
        self.lovd_file=lovd_file
        self.missingvariantfile=missingvariantfile

    def read_vcf_file(self):# this function is to read vcf files
        variants = []
        with open(self.vcf_file, 'r') as vcf:
            for line in vcf:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                chrom, pos, _, ref, alt, qual, filter_, info = fields[:8]
                variants.append({
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'qual': qual,
                    'filter': filter_,
                    'info': info
                }) #here we drop id which is useless in this project
        return variants
    
    def write_to_lovdjson(self,data):
        with open(self.lovd_file,"w") as json_file:
            json.dump(data,json_file)

    def get_vcf_calls(self):
        vcf_data=self.read_vcf_file()
        vcf_calls=set()#set is used to remove the duplicate calls
        for seq in vcf_data:
            calls=f"{self.source}:{seq['chrom']}:{seq['pos']}:{seq['ref']}:{seq['alt']}" #select the chrom, position, reference and alt
            vcf_calls.add(calls)    
        return list(vcf_calls)
    
    #this function is used to pull back all the content from api
    def query_API_lovd(self,genome_build,variant_description, transcript_model="all", select_transcripts="raw", checkonly="False", liftover="False"):
        #the base_url is the link to get json of each variant description
        base_url = f"https://rest.variantvalidator.org/LOVD/lovd/{genome_build}/{variant_description}/{transcript_model}/{select_transcripts}/{checkonly}/{liftover}?content-type=application/json"
        try:
            response = requests.get(base_url)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.HTTPError as http_err:
            print(f'HTTP error occurred: {http_err}')
        except requests.exceptions.RequestException as req_err:
            print(f'Request error occurred: {req_err}')
        except json.JSONDecodeError as json_err:
            print(f'JSON decode error occurred: {json_err}')
    
    ##this is the last funtion to merge all the funtionality and write the results to files
    def variant_validator_lovd_api_call(self):
        data_load=self.get_vcf_calls()
        variant_description=[]
        missing_variant=[]
        count=0
        for seq in data_load:
            variant_description.append(seq[7:]) #here we don't need the sources like GRCh37, so remove it from the results of get_vcf_calls function
        for selected_variant_description in variant_description:
            lovd_json=self.query_API_lovd(self.source,selected_variant_description,transcript_model="all", select_transcripts="raw", checkonly="False", liftover="False")
            format_json=self.format_data(selected_variant_description,lovd_json)
            if not lovd_json:
                missing_variant.append(selected_variant_description)
                print(missing_variant)
            else:
                if count==0: #when the content firstly wrote in the file, we need to clear all the content to aviod duplication
                    self.clear_write_file(format_json)
                else:
                    self.write_to_lovdjson(format_json)
            time.sleep(1)
            count=count+1
        print(missing_variant)
        self.write_missing_seq(missing_variant)

    
    def write_to_lovdjson(self,data):
        if not data:
            print('Data is empty or None. Skipping write operation.')
            return       
        try:
            with open(self.lovd_file, "a") as json_file:
                json.dump(data, json_file,indent=4)
                json_file.write('\n')
                print(f'Successfully wrote data to {self.lovd_file}')
        except IOError as io_err:
            print(f'IO error occurred while writing JSON file: {io_err}')
    
    def clear_write_file(self,data):
        with open(self.lovd_file, "w") as json_file:
            json.dump(data, json_file,indent=4)
            json_file.write('\n')
            print(f'Successfully wrote data to {self.lovd_file}')
        
    def write_missing_seq(self,data):
        with open(self.missingvariantfile,"w",encoding='utf-8') as file:
            for item in data:
                file.write(item+ '\n')

    #this function is used to select the variant from refseq and ensembl
    def format_data(self,seq,vcf_json):
        refseq_variant=[]
        ensembl_variant=[]
        if not vcf_json or not vcf_json.get(seq,{}) or not vcf_json.get(seq,{}).get(seq,{}) or not vcf_json.get(seq,{}).get(seq,{}).get('hgvs_t_and_p', {}):
            print(f'{seq} has no related variant')
            ensembl_variant.append(None)
            refseq_variant.append(None)
        else:
            for transript in vcf_json.get(seq,{}).get(seq,{}).get('hgvs_t_and_p', {}).values():
                if not transript or not transript.get('t_hgvs',{}):
                    print(f'{seq} has no related variant')
                else:
                    variant_judgement=transript.get('t_hgvs',{}).split(":")
                    if variant_judgement[1][0]=="c":#only focus on the sequence with c
                        if transript.get('t_hgvs',{}).startswith('E'):
                            ensembl_variant.append(transript.get('t_hgvs',{}))
                        if transript.get('t_hgvs',{}).startswith('N'):
                            refseq_variant.append(transript.get('t_hgvs',{}))
        format_data={
            f'<{seq}>':{
                'refseq':refseq_variant,
                'ensembl':ensembl_variant
            }
        }
        return format_data
    

if __name__=="__main__":
    filepath_GRCh37='C:/Dissertation_project/Clinical-Genomics-Variant/output/VCFs/GRCh37_indels.vcf'
    filepath_GRCh38='C:/Dissertation_project/Clinical-Genomics-Variant/output/VCFs/GRCh38_indels.vcf'
    lovd_json_path_GRCh37='C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/GRCh37lovd.json'
    lovd_json_path_GRCh38='C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/GRCh38lovd.json'
    missingvariantfile_GRCh37='C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/GRCh37_missing.txt'
    missingvariantfile_GRCh38='C:/Dissertation_project/Clinical-Genomics-Variant/LOVD_json/GRCh38_missing.txt'
    process_result_37=process_vcf(filepath_GRCh37,lovd_json_path_GRCh37,'GRCh37',missingvariantfile_GRCh37)
    process_result_38=process_vcf(filepath_GRCh38,lovd_json_path_GRCh38,'GRCh38',missingvariantfile_GRCh38)
    results_37=process_result_37.variant_validator_lovd_api_call()
    results_38=process_result_38.variant_validator_lovd_api_call()
