import os
import json
import csv
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm

def save_to_json(data, file_path):
    with open(file_path, mode='w', encoding='utf-8') as json_file:
        json.dump(data, json_file, ensure_ascii=False, indent=4)

def load_json(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

def get_files_with_prefix_suffix(directory, prefix='', suffix=''):
    matching_files = []
    for filename in os.listdir(directory):
        if filename.startswith(prefix) and filename.endswith(suffix):
            file = os.path.join(directory, filename)
            matching_files.append(file)
    return matching_files

def read_tsv_to_array(tsv_path, header=None):
    result = []
    with open(tsv_path, mode='r', newline='', encoding='utf-8') as file:
        first_line = file.readline().strip()
        columns = first_line.split('\t')
        file.seek(0)
        if header is not None:
            reader = csv.DictReader(file, delimiter='\t', fieldnames=header)
        else:
            reader = csv.DictReader(file, delimiter='\t')
        result = list(reader)
        if len(result) == 0:
            na_row = {col: 'NA' for col in columns}
            result.append(na_row)
    return result


if __name__ == "__main__":
    # Input and output path configuration
    gene_family_list_path = "../data/gene_family_list.csv"
    gene_family_source_dir = '../data/'
    gene_family_results_dir = '../results/'

    # Read gene family list from CSV
    gene_family_df = pd.read_csv(gene_family_list_path, index_col=None, header=0)
    gene_family_list = gene_family_df['genefamily'].tolist()

    # Create base results directory
    os.makedirs(gene_family_results_dir, exist_ok=True)

    # Process each gene family
    for index, gene_family_name in tqdm(enumerate(gene_family_list), total=len(gene_family_list)):
        # Set up directories
        gene_family_dir = os.path.join(gene_family_source_dir, gene_family_name)
        if not os.path.exists(gene_family_dir):
            continue

        diversity_stats_dir = os.path.join(gene_family_results_dir, gene_family_name, 'diversity_statistics')
        os.makedirs(diversity_stats_dir, exist_ok=True)

        # Read and process feature data
        feature_file_path = os.path.join(gene_family_dir, 'Feature.all.txt')
        feature_df = pd.read_csv(feature_file_path, sep='\t', index_col=0, header=0)
        
        # Clean and standardize data
        feature_df['mRNA_ID'] = [x.split('.')[0] for x in feature_df['mRNA_ID']]
        feature_df['Gene_symbol'] = feature_df['Gene_symbol'].fillna('Unknown')
        feature_df['Gene_symbol'] = feature_df['Gene_symbol'].str.upper()

        # Get unique gene symbols and species
        unique_gene_symbols = sorted(list(set(feature_df['Gene_symbol'])))
        unique_species = sorted(list(set(feature_df.index)))

        # Create species-gene member matrix
        species_gene_matrix = pd.DataFrame(columns=['species'] + unique_gene_symbols, 
                                         index=unique_species)
        species_gene_matrix['species'] = species_gene_matrix.index
        species_gene_matrix.fillna(0, inplace=True)

        # Count gene occurrences per species
        for index, feature_row in feature_df.iterrows():
            species_gene_matrix.loc[index, feature_row['Gene_symbol']] += 1

        # Save species-gene matrix
        matrix_output_path = os.path.join(diversity_stats_dir, 
                                        f'{gene_family_name}.species.feature.all.flow.txt')
        species_gene_matrix.to_csv(matrix_output_path, sep='\t', index=False)

        # Read and format data for JSON output
        species_gene_data = read_tsv_to_array(matrix_output_path)
        matrix_columns = list(species_gene_data[0].keys())
        
        # Save JSON format data
        json_output_path = os.path.join(diversity_stats_dir, 
                                      f'{gene_family_name}_species_tree_link_flow_data.json')
        json_data = {
            'data': species_gene_data,
            'columns': matrix_columns
        }
        save_to_json(json_data, json_output_path)



