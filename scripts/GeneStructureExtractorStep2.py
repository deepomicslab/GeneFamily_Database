import csv
import os
import pandas as pd
from tqdm import tqdm
import json
import re

def read_csv_to_array(csv_path, header=None):
    result = []
    with open(csv_path, mode='r', newline='', encoding='utf-8') as file:
        if header is not None:
            reader = csv.DictReader(file, fieldnames=header)
        else:
            reader = csv.DictReader(file)
        for row in reader:
            result.append(row)
    return result

def get_files_with_prefix_suffix(directory, prefix='', suffix=''):
    matching_files = []
    for filename in os.listdir(directory):
        if filename.startswith(prefix) and filename.endswith(suffix):
            file = os.path.join(directory, filename)
            matching_files.append(file)
    return matching_files

def save_to_json(data, file_path):
    with open(file_path, mode='w', encoding='utf-8') as json_file:
        json.dump(data, json_file, ensure_ascii=False, indent=4)

def load_json(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

def analysis_gene_structure(database_dir, gene_family_name, species, module_name, default_map=True):
    gene_family_dir = os.path.join(database_dir, module_name, gene_family_name)
    corresponding_files = [os.path.join(gene_family_dir, f'{specie}_gene_exon_info.csv') for specie in species]
    tmp_data_format = {'columns':[], 'data':[]}
    file_header = [ "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    ]
    return_data = []
    for corresponding_file in corresponding_files:
        tmp_data_format = {'columns':[], 'data':[]}
        tmp_data_format['columns'] = file_header
        tmp_structure_data = read_csv_to_array(corresponding_file, file_header)
        tmp_data_format['data'] = [{**x} for x in tmp_structure_data]
        return_data.append(tmp_data_format)

    return_file_list = [ os.path.basename(x).replace("_gene_exon_info.csv", "") for x in corresponding_files]

    return {
        'file_name_list': return_file_list,
        'data': return_data
    }

def get_all_gene_structure(database_dir, gene_family_name, module_name,):
    gene_family_dir = os.path.join(database_dir, module_name, gene_family_name)
    corresponding_files = get_files_with_prefix_suffix(gene_family_dir, prefix='', suffix='_gene_exon_info.csv')
    tmp_data_format = {'columns':[], 'data':[]}
    file_header = [ "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    ]

    return_data = []

    for corresponding_file in corresponding_files:
        tmp_data_format = {'columns':[], 'data':[]}
        tmp_data_format['columns'] = file_header
        tmp_data_format['data'] = read_csv_to_array(corresponding_file, file_header)
        return_data.extend(tmp_data_format['data'])

    return_data_format = file_header

    return {
        'header': return_data_format,
        'data': return_data
    }

def analysis_gene_structure_species_sections(database_dir, gene_family_name, module_name):
    gene_family_dir = os.path.join(database_dir, module_name, gene_family_name)
    files = get_files_with_prefix_suffix(gene_family_dir, prefix='', suffix='_gene_exon_info.csv')
    return [os.path.basename(x).replace("_gene_exon_info.csv", "") for x in files]

def get_transcript_id_mrna(attributes):
    match = re.search(r'ID=transcript:([^;]+)', attributes)
    if match:
        return match.group(1)
    return None

def get_parent_transcript_cds(attributes):
    if not isinstance(attributes, (str, bytes)):
        return None
    try:
        match = re.search(r'Parent=transcript:([^;]+)', attributes)
        if match:
            return match.group(1)
    except TypeError:
        return None
    return None


if __name__ == "__main__":
    # Input paths configuration
    gene_family_list_path = "../data/gene_family_list.csv"
    gene_family_module3_dir = "../data/"
    gene_family_module1_dir = "../data/"
    gene_family_results_dir = r'../results/gene_structure/'
    
    # Get module name from directory
    module_name = os.path.basename(gene_family_module3_dir)

    # Read gene family list
    gene_family_df = pd.read_csv(gene_family_list_path, index_col=None, header=0)
    gene_family_list = gene_family_df['genefamily'].tolist()

    # Process each gene family's structure
    for idx, gene_family_name in tqdm(enumerate(gene_family_list), total=len(gene_family_list)):        
        # Skip if gene family directory doesn't exist
        if not os.path.exists(os.path.join(gene_family_module3_dir, module_name, gene_family_name)):
            continue

        print('gene_family_name', gene_family_name)

        # Set up output directory
        gene_structure_dir = os.path.join(gene_family_results_dir, gene_family_name, "Module3", "gene_structure")
        os.makedirs(gene_structure_dir, exist_ok=True)

        # Analyze gene structure sections by species
        gene_structure_species_sections = analysis_gene_structure_species_sections(gene_family_module3_dir, 
                                                                                gene_family_name, 
                                                                                module_name)

        # Process feature data
        feature_file_path = os.path.join(gene_family_module1_dir, gene_family_name, 'Feature.all.txt')
        if not os.path.exists(feature_file_path):
            continue

        # Read and process feature data
        feature_df = pd.read_csv(feature_file_path, sep='\t', index_col=None, header=0)
        feature_df['mRNA_ID'] = [x.split('.')[0] for x in feature_df['mRNA_ID']]
        feature_df['Gene_symbol'].fillna('Unknown', inplace=True)
        feature_df['Gene_symbol'] = feature_df['Gene_symbol'].str.upper()
        feature_df['Gene_symbol'] = feature_df['Gene_symbol'].str.replace('/', '_', regex=False)
        feature_df['Protein_ID'] = [x.split('.')[0] for x in feature_df['Protein_ID']]

        # Create mapping dictionaries
        mrna_to_gene_symbol = dict(zip(feature_df['mRNA_ID'], feature_df['Gene_symbol']))
        mrna_to_species = dict(zip(feature_df['mRNA_ID'], feature_df['Species']))

        # Process gene structure by species
        for species_section in gene_structure_species_sections:
            response_data = analysis_gene_structure(gene_family_module3_dir, 
                                                 gene_family_name, 
                                                 [species_section], 
                                                 module_name)
            
            species_feature_df = feature_df[feature_df['Species'] == species_section]
            species_mrna_to_gene = dict(zip(feature_df['mRNA_ID'], feature_df['Gene_symbol']))
            response_data['mrnap_genesymbols'] = species_mrna_to_gene
            
            output_path = os.path.join(gene_structure_dir, 
                                     f'gene_structure_result.{gene_family_name}.{species_section}.json')
            save_to_json(response_data, output_path)

        # Process gene structure by gene symbol
        gene_symbol_list = list(feature_df['Gene_symbol'])

        # Get all gene structure data
        all_structure_data = get_all_gene_structure(gene_family_module3_dir, gene_family_name, module_name)
        structure_records = all_structure_data['data']
        structure_headers = all_structure_data['header']

        # Group results by family members
        family_member_results = {}

        # Process each structure record
        for i, record in enumerate(structure_records):
            mrna_id = None
            symbol = None
            
            if record['type'] == 'mRNA':
                mrna_id = get_transcript_id_mrna(record.get('attributes', ''))
            elif record['type'] in ['CDS', 'five_prime_UTR', 'three_prime_UTR']:
                mrna_id = get_parent_transcript_cds(record.get('attributes', ''))
                
            symbol = mrna_to_gene_symbol.get(mrna_id, "UNKNOWN")
            
            if symbol not in family_member_results:
                family_member_results[symbol] = []
            family_member_results[symbol].append(record)

        # Save results for each gene symbol
        for gene_symbol, structure_data in family_member_results.items():
            symbol_feature_df = feature_df[feature_df['Gene_symbol'] == gene_symbol]
            symbol_mrna_to_species = dict(zip(feature_df['mRNA_ID'], feature_df['Species']))
            
            symbol_result = {
                'file_name_list': [gene_symbol],
                'data': [{
                    'columns': structure_headers,
                    'data': structure_data,
                }],
                'mrnap_species': symbol_mrna_to_species
            }
            
            output_path = os.path.join(gene_structure_dir,
                                     f'gene_structure_result.{gene_family_name}.{gene_symbol}.json')
            save_to_json(symbol_result, output_path)

