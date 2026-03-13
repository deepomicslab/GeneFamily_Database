from tqdm import tqdm
import pandas as pd
import os
import json

def save_to_json(data, file_path):
    with open(file_path, mode='w', encoding='utf-8') as json_file:
        json.dump(data, json_file, ensure_ascii=False, indent=4)

if __name__ == "__main__":
    # Input paths
    gene_family_list_path = "../data/gene_family_list.csv"
    gene_family_source_dir = r"../data"
    gene_family_results_dir = r'../results'

    # Read gene family list from CSV
    gene_family_df = pd.read_csv(gene_family_list_path, index_col=None, header=0)
    gene_family_list = gene_family_df['genefamily'].tolist()

    # Dictionary to store gene family members
    gene_family_members = {}

    # Process each gene family
    for index, gene_family_name in tqdm(enumerate(gene_family_list), total=len(gene_family_list)):
        # Check if gene family directory exists
        gene_family_dir = os.path.join(gene_family_source_dir, gene_family_name)
        if not os.path.exists(gene_family_dir):
            continue
        
        # Read feature data
        feature_file_path = os.path.join(gene_family_dir, 'Feature.all.txt')
        feature_df = pd.read_csv(feature_file_path, sep='\t', index_col=0, header=0)
        
        # Process gene symbols
        feature_df['Gene_symbol'].fillna('Unknown', inplace=True)
        feature_df['Gene_symbol'] = feature_df['Gene_symbol'].str.upper()

        # Extract unique gene symbols excluding 'UNKNOWN'
        unique_gene_symbols = feature_df['Gene_symbol'].unique().tolist()
        gene_family_members[gene_family_name] = [
            symbol for symbol in unique_gene_symbols if symbol != 'UNKNOWN'
        ]

    # Save results
    output_file_path = os.path.join(gene_family_results_dir, 'gene_family_members.json')
    save_to_json(gene_family_members, output_file_path)



