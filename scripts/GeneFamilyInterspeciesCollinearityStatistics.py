import csv
import os
import json
import shutil
import sys
import subprocess
import copy
from collections import Counter
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import difflib

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

def count_gene_occurrences(data):
    result = {}
    for species, genes in data.items():
        result[species] = dict(Counter(genes))
    return result

def count_gene_occurrences(data):
    result = {}
    for species, genes in data.items():
        result[species] = dict(Counter(genes))
    return result

def analyze_gene_family_evolution(genefamily_name, database_dir):
    """
    Analyze gene family evolution and generate connection matrices.

    Args:
        module_name (str): Name of the module
        genefamily_name (str): Name of the gene family
        database_dir (str): Path to the database directory

    Returns:
        tuple: (sp1_sp2_combination_same_genemembers, sp1_sp2_combination_diff_genemembers, sp_genemembers)
    """
    genefamily_dir = os.path.join(database_dir, genefamily_name, 'Module4', 'evolution')

    # Read and process Feature.all.txt
    feature_all_df = pd.read_csv(os.path.join(genefamily_dir, 'Feature.all.txt'), sep='\t', index_col=0, header=0)
    feature_all_df['mRNA_ID'] = [x.split('.')[0] for x in feature_all_df['mRNA_ID']]
    feature_all_df['Gene_symbol'] = feature_all_df['Gene_symbol'].fillna('UNKNOWN')
    feature_all_df['Gene_symbol'] = feature_all_df['Gene_symbol'].str.upper()
    genefamily_pep_genesymbols = dict(zip(feature_all_df['mRNA_ID'], feature_all_df['Gene_symbol']))

    # Get family members
    genefamily_members = list(set(genefamily_pep_genesymbols.values()))

    # Process link files
    links_files = get_files_with_prefix_suffix(genefamily_dir, prefix='', suffix='.link.txt')
    filter_links_files = [x for x in links_files if '.tandem.link.txt' not in os.path.basename(x)]

    sp1_combinations = []
    sp2_combinations = []
    sp1mrna_sp2mrna_combinations = []
    sp_counts = {}
    sp_genecounts = {}
    sp_list = []
    sp_genemembers = pd.DataFrame(columns=['Species', 'Species Counts'] + genefamily_members)

    # Process each link file
    for link_file in filter_links_files:
        tmp_link_data = read_tsv_to_array(link_file)
        if len(tmp_link_data) == 1 and list(tmp_link_data[0].keys())[0] == '':
            continue
        if 'Species_1' not in tmp_link_data[0]:
            continue
        sp1_combinations.append(tmp_link_data[0]['Species_1'])
        sp2_combinations.append(tmp_link_data[0]['Species_2'])
        sp_list.append(tmp_link_data[0]['Species_1'])
        sp_list.append(tmp_link_data[0]['Species_1'])

        for tmp_link_row in tmp_link_data:

            if tmp_link_row['Species_1'] not in sp_genecounts:
                sp_genecounts[tmp_link_row['Species_1']] = []
            if tmp_link_row['Species_2'] not in sp_genecounts:
                sp_genecounts[tmp_link_row['Species_2']] = []

            if (tmp_link_row['mRNAID_1'] in genefamily_pep_genesymbols and
                tmp_link_row['mRNAID_2'] in genefamily_pep_genesymbols):

                if (genefamily_pep_genesymbols[tmp_link_row['mRNAID_1']] == 'UNKNOWN' and
                    genefamily_pep_genesymbols[tmp_link_row['mRNAID_2']] == 'UNKNOWN') or (genefamily_pep_genesymbols[tmp_link_row['mRNAID_1']] != genefamily_pep_genesymbols[tmp_link_row['mRNAID_2']] and
                            (genefamily_pep_genesymbols[tmp_link_row['mRNAID_1']] != 'UNKNOWN') and (genefamily_pep_genesymbols[tmp_link_row['mRNAID_1']] != 'UNKNOWN') ) :
                    continue

                sp1mrna_sp2mrna_combinations.append({
                    f"{tmp_link_row['Species_1']}&&{tmp_link_row['Species_2']}":
                    f"{tmp_link_row['mRNAID_1']}&&{tmp_link_row['mRNAID_2']}"
                })

                sp_genecounts[tmp_link_row['Species_1']].append(
                    genefamily_pep_genesymbols[tmp_link_row['mRNAID_1']]
                )
                sp_genecounts[tmp_link_row['Species_2']].append(
                    genefamily_pep_genesymbols[tmp_link_row['mRNAID_2']]
                )

                # Count species
                for sp in [tmp_link_row['Species_1'], tmp_link_row['Species_2']]:
                    if sp not in sp_counts:
                        sp_counts[sp] = []
                sp_counts[tmp_link_row['Species_1']].append(tmp_link_row['Species_2'])
                sp_counts[tmp_link_row['Species_2']].append(tmp_link_row['Species_1'])
                sp_counts[tmp_link_row['Species_1']] = list(set(sp_counts[tmp_link_row['Species_1']]))
                sp_counts[tmp_link_row['Species_2']] = list(set(sp_counts[tmp_link_row['Species_2']]))

    sp_genecounts_result = count_gene_occurrences(sp_genecounts)

    # Generate species gene members table
    sp_list = list(set(sp_list))
    sp_list.sort()
    sp_genemembers['Species'] = sp_list
    sp_genemembers['Species_Counts'] = [len(sp_counts[x]) if x in sp_counts else 0 for x in sp_list]
    sp_genemembers = sp_genemembers.fillna(0)
    sp_genemembers.index = sp_genemembers['Species']

    for index, row in sp_genemembers.iterrows():
        for column in row.index[2:]:
            if column in sp_genecounts_result[row['Species']]:
                sp_genemembers.loc[row['Species'], column] += sp_genecounts_result[row['Species']][column]

    # Filter and clean sp_genemembers
    mask = (sp_genemembers[genefamily_members] != 0).any(axis=1)
    sp_genemembers = sp_genemembers[mask]
    sp_genemembers = sp_genemembers.loc[:, ~(sp_genemembers == 0).all()]

    # Generate combination tables
    header = ['Species_1', 'Species_2'] + genefamily_members
    sp1_sp2_combination_same_genemembers = pd.DataFrame(columns=header)
    sp1_sp2_combination_diff_unknown_genemembers = pd.DataFrame(columns=header)
    sp1_sp2_combination_diff_genesymbol_genemembers = pd.DataFrame(columns=header)

    for df in [sp1_sp2_combination_same_genemembers, sp1_sp2_combination_diff_unknown_genemembers, sp1_sp2_combination_diff_genesymbol_genemembers]:
        df['Species_1'] = sp1_combinations
        df['Species_2'] = sp2_combinations
        df.index = [f'{x}&&{sp2_combinations[i]}' for i, x in enumerate(sp1_combinations)]
        df.fillna(0, inplace=True)

    # Process combinations
    for item in sp1mrna_sp2mrna_combinations:
        values = list(item.values())
        keys = list(item.keys())
        mrna1, mrna2 = values[0].split('&&')

        if mrna1 not in genefamily_pep_genesymbols or mrna2 not in genefamily_pep_genesymbols:
            continue

        mrna1_symbol = genefamily_pep_genesymbols[mrna1]
        mrna2_symbol = genefamily_pep_genesymbols[mrna2]

        if mrna1_symbol != mrna2_symbol and (mrna1_symbol == 'UNKNOWN' or mrna2_symbol == 'UNKNOWN'):
            sp1_sp2_combination_diff_unknown_genemembers.loc[keys[0], mrna1_symbol] += 1
            sp1_sp2_combination_diff_unknown_genemembers.loc[keys[0], mrna2_symbol] += 1
        if mrna1_symbol != mrna2_symbol and mrna1_symbol != 'UNKNOWN' and mrna2_symbol != 'UNKNOWN':
            sp1_sp2_combination_diff_genesymbol_genemembers.loc[keys[0], mrna1_symbol] += 1
            sp1_sp2_combination_diff_genesymbol_genemembers.loc[keys[0], mrna2_symbol] += 1
        if mrna1_symbol == mrna2_symbol and mrna1_symbol != 'UNKNOWN' and mrna2_symbol != 'UNKNOWN':
            sp1_sp2_combination_same_genemembers.loc[keys[0], mrna1_symbol] += 1


    mask = (sp1_sp2_combination_same_genemembers[genefamily_members] != 0).any(axis=1)
    sp1_sp2_combination_same_genemembers = sp1_sp2_combination_same_genemembers[mask]
    sp1_sp2_combination_same_genemembers = sp1_sp2_combination_same_genemembers.loc[:, ~(sp1_sp2_combination_same_genemembers == 0).all()]


    mask = (sp1_sp2_combination_diff_unknown_genemembers[genefamily_members] != 0).any(axis=1)
    sp1_sp2_combination_diff_unknown_genemembers = sp1_sp2_combination_diff_unknown_genemembers[mask]
    sp1_sp2_combination_diff_unknown_genemembers = sp1_sp2_combination_diff_unknown_genemembers.loc[:, ~(sp1_sp2_combination_diff_unknown_genemembers == 0).all()]

    mask = (sp1_sp2_combination_diff_genesymbol_genemembers[genefamily_members] != 0).any(axis=1)
    sp1_sp2_combination_diff_genesymbol_genemembers = sp1_sp2_combination_diff_genesymbol_genemembers[mask]
    sp1_sp2_combination_diff_genesymbol_genemembers = sp1_sp2_combination_diff_genesymbol_genemembers.loc[:, ~(sp1_sp2_combination_diff_genesymbol_genemembers == 0).all()]

    return (sp1_sp2_combination_same_genemembers,
            sp1_sp2_combination_diff_unknown_genemembers,
            sp1_sp2_combination_diff_genesymbol_genemembers,
            sp_genemembers)


def process_and_save_gene_family_evolution(genefamily_name, database_dir, save_dir):
    """
    Process gene family evolution analysis and save results to files.

    Args:
        genefamily_name (str): Name of the gene family (e.g., 'Lys')
        module_name (str): Name of the module (e.g., 'module4')
        database_dir (str): Path to the database directory
        save_dir (str): Path to save the output files
    """
    # Run analysis
    same_genemembers, diff_unknown_genemembers, diff_genesymbol_genemembers, sp_genemembers = analyze_gene_family_evolution(
        genefamily_name=genefamily_name,
        database_dir=database_dir
    )

    # Create save directory
    os.makedirs(save_dir, exist_ok=True)

    # Save CSV files
    same_genemembers.to_csv(
        os.path.join(save_dir, f"{genefamily_name}.links.flow.txt"),
        sep='\t',
        index=False
    )

    diff_unknown_genemembers.to_csv(
        os.path.join(save_dir, f"{genefamily_name}.links.flow.unknown.txt"),
        sep='\t',
        index=False
    )


    sp_genemembers.to_csv(
        os.path.join(save_dir, f"{genefamily_name}.species.links.flow.txt"),
        sep='\t',
        index=False
    )

    # Read and process same gene members data
    sp1_sp2_combination_same_genemembers_df = read_tsv_to_array(
        os.path.join(save_dir, f"{genefamily_name}.links.flow.txt")
    )
    sp1_sp2_combination_same_genemembers_columns = list(sp1_sp2_combination_same_genemembers_df[0].keys())

    # Save JSON file
    file_name = f'{genefamily_name}_sp1_sp2_combination_same_genemembers.json'
    output_dir = os.path.join(save_dir, file_name)
    save_to_json({
        'data': sp1_sp2_combination_same_genemembers_df,
        'columns': sp1_sp2_combination_same_genemembers_columns
    }, output_dir)

def process_gene_family(database_dir, gene_family_name, base_save_dir):
    """Integrate and process three analysis modules for a single gene family
    
    Args:
        database_dir: Database directory path
        module_name: Module name
        gene_family_name: Gene family name
        base_save_dir: Base directory path for saving results
    """
    
    try:
        # Check if gene family directory exists
        gene_family_dir = os.path.join(database_dir, gene_family_name)
        if not os.path.exists(gene_family_dir):
            return False            
        module4_dir = os.path.join(base_save_dir, gene_family_name)
        save_dir_evolution_flow = os.path.join(module4_dir, "interspecies_collinearity_statistics") 
        process_and_save_gene_family_evolution(
            genefamily_name=gene_family_name,
            database_dir=database_dir,
            save_dir=save_dir_evolution_flow
        )
        
        return True
        
    except Exception as e:
        print(f"Error processing {gene_family_name}: {str(e)}")
        return False


if __name__ == '__main__':
    # Input and output path configuration
    gene_family_list_path = "../data/gene_family_list.csv"
    gene_family_source_dir = '../data'
    gene_family_output_dir = '../results'

    # Read gene family list from CSV
    gene_family_df = pd.read_csv(gene_family_list_path, index_col=None, header=0)
    gene_family_list = gene_family_df['genefamily'].tolist()

    # Process each gene family
    for index, gene_family_name in tqdm(enumerate(gene_family_list), 
                                      total=len(gene_family_list), 
                                      desc="Processing genes"):
        process_gene_family(gene_family_source_dir, 
                          gene_family_name, 
                          gene_family_output_dir)


