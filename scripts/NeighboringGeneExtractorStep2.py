import csv
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import copy
import re
import json
import sys
import pickle
from multiprocessing import Pool, cpu_count
from collections import defaultdict
import os
import shutil
import subprocess


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

def extract_transcript_id(attributes):
    """Extract transcript ID from attributes string
    
    Searches for transcript ID in attributes string using either 'ID=transcript:' 
    or 'Parent=transcript:' pattern. Returns the matched transcript ID or None if not found.
    """
    transcript_id = None
    # Try to find ID=transcript: pattern
    id_match = re.search(r'ID=transcript:([^;]+)', attributes)
    if id_match:
        transcript_id = id_match.group(1)
    else:
        # If not found, try to find Parent=transcript: pattern
        parent_match = re.search(r'Parent=transcript:([^;]+)', attributes)
        if parent_match:
            transcript_id = parent_match.group(1)
    return transcript_id

def parse_attributes(attr_str):
    """Parse attribute string
    
    Parses a GFF/GTF format attribute string into a dictionary
    containing transcriptID, geneID, and geneName
    """
    attrs = {}
    # Split attribute string by semicolon and remove empty strings
    items = [x.strip() for x in attr_str.split(';') if x.strip()]
    
    # Extract specific attributes
    for item in items:
        if 'ID=transcript:' in item:
            attrs['transcriptID'] = item.replace("ID=transcript:","")
        if 'Parent=gene:' in item:
            attrs['geneID'] = item.replace("Parent=gene:","")
        if 'Name=' in item:
            attrs['geneName'] = item.replace("Name=","")
    
    # Set default empty values for missing attributes
    for key in ['transcriptID', 'geneID', 'geneName']:
        if key not in attrs:
            attrs[key] = ''
    
    return attrs

def process_gene_exon_csv(csv_path, file_header=None):
    """Process gene exon CSV file
    
    Reads and processes a CSV file containing gene exon information,
    groups the data by transcript ID, and formats it into a structured output
    """
    # Read CSV data into array
    exon_info_data = read_csv_to_array(csv_path, file_header)
    df = pd.DataFrame(exon_info_data)
    
    # Extract transcript IDs from attributes column
    df['transcript_id'] = df['attributes'].apply(extract_transcript_id)
    
    # Group data by transcript ID
    grouped = df.groupby('transcript_id')
    
    # Convert each group to dictionary records
    complete_transcripts = []
    for transcript_id, group in grouped:
        complete_transcripts.append(group.to_dict(orient="records"))
    
    # Process and format the results
    results = []
    for item in complete_transcripts:
        results.append([])
        for block in item:
            # Parse attributes string into dictionary
            attrs = parse_attributes(block['attributes'])
            
            # Create formatted result row
            result_row = {
                'chr': block['chromosome'],
                'type': block['type'],
                'start': int(block['start']),
                'end': int(block['end']),
                'strand': block['strand'],
                'transcriptID': attrs['transcriptID'],
                'geneID': attrs['geneID'],
                'geneName': attrs['geneName'].upper().split('-')[0],
                'species': os.path.basename(csv_path).replace('_gene_exon_info.csv', ''),
            }
            results[-1].append(result_row)
    
    return results

def analysis_genomic_neighbors(database_dir, gene_family_name, correspondind_species, module_name, default_map=True):
    """Analyze genomic neighboring regions"""
    gene_family_dir = os.path.join(database_dir, module_name, gene_family_name)
    print("gene_family_dir ", gene_family_dir)
    
    corresponding_files = []
    if default_map:
        # Get all files with '_gene_exon_info.csv' suffix if using default mapping
        files = get_files_with_prefix_suffix(gene_family_dir, prefix='', suffix='_gene_exon_info.csv')
        corresponding_files = files
    else:
        # Get files only for specified species if not using default mapping
        for species in correspondind_species:
            files = get_files_with_prefix_suffix(gene_family_dir, prefix=species, suffix='_gene_exon_info.csv')
            corresponding_files.extend(files)

    # Define GFF file header columns
    file_header = [
        "chromosome", "source", "type", "start", "end",
        "score", "strand", "phase", "attributes"
    ]
    
    return_data = []
    return_file_list = []
    
    for corresponding_file in corresponding_files:
        tmp_data_format = {'columns':[], 'data':[]}
        test = read_csv_to_array(corresponding_file, file_header)
        if len(test) == 0:
            continue
            
        # Process gene exon information from CSV file
        tmp_mrna_blocks = process_gene_exon_csv(corresponding_file, file_header)
        tmp_data_format['data'] = tmp_mrna_blocks
        tmp_data_format['columns'] = file_header
        return_data.append(tmp_data_format)
        return_file_list.append(os.path.basename(corresponding_file))

    # Return processed data and list of processed files
    return {
        'file_name_list': return_file_list,
        'data': return_data
    }

def export_gene_data_to_tsv(gene_species_response_data, result_dir):
    """Export gene data to TSV files"""
    if 'data' not in gene_species_response_data:
        print("Error: Input data missing 'data' key!")
        return
    
    for gene, gene_data in gene_species_response_data['data'].items():
        for current_gene_id, current_gene_data in gene_data.items():
            if 'data' not in current_gene_data:
                print(f"Warning: gene '{gene}' with current_gene_id '{current_gene_id}' missing 'data' key, skipping!")
                continue
            
            data_list = current_gene_data['data']
            if not data_list:
                print(f"Warning: gene '{gene}' with current_gene_id '{current_gene_id}' has empty data, skipping!")
                continue
            
            fieldnames = list(data_list[0].keys())
            
            # Convert dictionary values to JSON strings
            for item in data_list:
                for key in item:
                    if isinstance(item[key], dict):
                        item[key] = json.dumps(item[key])
            
            output_filename = os.path.join(result_dir, f"{gene}-{current_gene_id}.tsv")
            
            with open(output_filename, 'w', newline='', encoding='utf-8') as tsv_file:
                writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter='\t')
                writer.writeheader()
                for item in data_list:
                    writer.writerow(item)
            
            print(f"File written: {output_filename}")

def transform_gene_dict(input_dict):
    """
    Transform gene data dictionary
    Converts a dictionary structured by gene_symbol -> species_gene_key -> gene_info
    into a dictionary structured by species -> combined_gene_key -> gene_info
    """
    result = {}
    for gene_symbol, gene_data in input_dict.items():
        for species_gene_key, gene_info in gene_data.items():
            species, gene_id = species_gene_key.split('&', 1)
            if species not in result:
                result[species] = {}
            new_key = f"{gene_symbol}&{gene_id}"
            result[species][new_key] = gene_info
    return result

def process_gene_family(args):
    """Process a single gene family"""
    (gene_family_name, module_name, results_base_dir, download_base_dir, base_dir, species_gff_dir) = args
    
    try:
        # Set save path
        save_dir = os.path.join(results_base_dir, 
                               gene_family_name, "Module3", "gene_neighbor")
        # Directory for downloaded files package
        result_dir = os.path.join(download_base_dir, 
                                 'download', gene_family_name, 'gene_neighborhood_pre')
        
        gene_family_data_dir = os.path.join(base_dir, module_name, gene_family_name)

        if not os.path.exists(gene_family_data_dir):
            return f"Skipped {gene_family_name}: Data directory not found at {gene_family_data_dir}"
        
        # Check if results already exist
        file_path = os.path.join(results_base_dir, gene_family_name, "Module3", f"{gene_family_name}_gene_neighbor.tar.gz")
        if os.path.exists(file_path) and os.path.getsize(file_path) > 390:
            print('???', os.path.join(results_base_dir, gene_family_name, "Module3"))
            return 'has results'
        
        # Create directories
        os.makedirs(save_dir, exist_ok=True)
        os.makedirs(result_dir, exist_ok=True)
        
        # Analyze genomic neighboring regions
        correspondind_species = []
        # Data folder location
        response_data = analysis_genomic_neighbors(base_dir, gene_family_name,
                                                 correspondind_species, module_name)

        # Process gene list
        all_gene_list = []
        all_block_items = []
        for seqs_data in response_data['data']:
            for seqs in seqs_data['data']:
                for gene in seqs:
                    gene['species'] = os.path.basename(gene['species'])
                    all_gene_list.append(gene['geneName'])
                    all_block_items.append(gene)

        # Get gene set and species list
        all_gene_set = [x.upper().split('-')[0] for x in list(set(all_gene_list)) if x != '']
        species_list = [x.replace('_gene_exon_info.csv','') for x in response_data['file_name_list']]

        # Initialize response data
        gene_species_response_data = {'data': {}}
        gene_species_response_data_download = {'data': {}}
        range_width = 1000000

        # Process each gene
        for gene in all_gene_set:
            
            #print('gene', gene)

            gene_species_response_data['data'][gene] = {}
            gene_species_response_data_download['data'][gene] = {}
            print(f'Processing gene: {gene} in family: {gene_family_name}')

            # Process each species
            for specie in species_list:

                print(f'Processing gene: {gene} species: {specie} in family: {gene_family_name}')

                # Read GFF data
                gff_blocks_path = os.path.join(species_gff_dir, 
                                             f'{specie}_species_gff_block_items.pkl')
                gff_download_path = os.path.join(species_gff_dir, 
                                               f'{specie}_species_gff_download.pkl')
                
                try:
                    with open(gff_blocks_path, 'rb') as f:
                        gff_blocks = pickle.load(f)
                    with open(gff_download_path, 'rb') as f:
                        gff_download = pickle.load(f)

                except Exception as e:
                    print(f"Error loading GFF data for species {specie}: {e}")
                    continue

                # Process all blocks
                for item in all_block_items:
                    #print('item ', item)
                    tmp_species_items = []
                    tmp_species_items_download = []
                    current_gene_id = ''

                    if item['geneName'].upper() == gene and item['species'] == specie:
                        current_gene_id = f"{item['species']}&{item['geneID']}&{item['transcriptID']}"

                        # Initialize data structures
                        gene_species_response_data['data'][gene][current_gene_id] = {'data':[], 'pos_range': []}
                        gene_species_response_data_download['data'][gene][current_gene_id] = {'data':[]}
                        
                        # Calculate position range
                        pos_start = max([item['start'] - range_width, 0])
                        pos_end = item['end'] + range_width
                        pos_chr = item['chr'].split('.')[0]

                        gene_species_response_data['data'][gene][current_gene_id]['pos_range'] = [pos_start, pos_end]

                        # Process GFF blocks
                        #print('chr block counts ', len(gff_blocks[pos_chr]))
                        for subitem in gff_blocks[pos_chr]:
                            subitem['chr'] = subitem['chr'].split('.')[0]
                            if (subitem['species'] == specie and 
                                subitem['chr'] == pos_chr and 
                                subitem['start'] >= pos_start and 
                                subitem['end'] <= pos_end):
                                tmp_species_items.append(subitem)

                        for subitem in all_block_items:
                            subitem['chr'] = subitem['chr'].split('.')[0]
                            if (subitem['species'] == specie and 
                                subitem['chr'] == pos_chr and 
                                subitem['start'] >= pos_start and 
                                subitem['end'] <= pos_end):
                                tmp_species_items.append(subitem)

                        for subitem_download in gff_download[pos_chr]:
                            subitem_download['chromosome'] = subitem_download['chromosome'].split('.')[0]
                            if (subitem_download['chromosome'] == pos_chr and 
                                int(subitem_download['start']) >= pos_start and 
                                int(subitem_download['end']) <= pos_end):
                                tmp_species_items_download.append(subitem_download)

                    if current_gene_id in gene_species_response_data['data'][gene]:
                        gene_species_response_data['data'][gene][current_gene_id]['data'].extend(tmp_species_items)
                    if current_gene_id in gene_species_response_data_download['data'][gene]:
                        gene_species_response_data_download['data'][gene][current_gene_id]['data'].extend(tmp_species_items_download)

        save_to_json(gene_species_response_data, 
                    os.path.join(save_dir, f'analysis_genomic_neighbors_response.{gene_family_name}.gene.species.json'))

        export_gene_data_to_tsv(gene_species_response_data_download, result_dir)

        # Save gene data
        for gene_data, species_info in gene_species_response_data['data'].items():
            print(f"Genes: {gene_data}, Number of corresponding species: {len(species_info)}")
            tmp_gene_species_data = {'data': {}}
            tmp_gene_species_data['data'][gene_data] = species_info
            save_to_json(tmp_gene_species_data, 
                        os.path.join(save_dir, f'analysis_genomic_neighbors_response.{gene_family_name}.{gene_data}.gene.species.json'))

        # Transform and save species data
        transformed_data = transform_gene_dict(gene_species_response_data['data'])
        for species_data, gene_info in transformed_data.items():
            print(f"Species: {species_data}, Number of corresponding genes: {len(gene_info)}")
            tmp_species_gene_data = {'data': {}}
            tmp_species_gene_data['data'][species_data] = gene_info
            save_to_json(tmp_species_gene_data, 
                        os.path.join(save_dir, f'analysis_genomic_neighbors_response.{gene_family_name}.{species_data}.species.gene.json'))

        save_to_json(transformed_data, 
                    os.path.join(save_dir, f'analysis_genomic_neighbors_response.{gene_family_name}.species.gene.json'))

        # Save sections
        save_sections = {
            'gene_sections': [{"label":x, "value":x} for x in sorted(list(gene_species_response_data['data'].keys()))],
            'species_sections': [{"label":x, "value":x} for x in sorted(list(transformed_data.keys()))]
        }
        save_to_json(save_sections, 
                    os.path.join(save_dir, f'analysis_genomic_neighbors_response.{gene_family_name}.sections.json'))

        return f"Successfully processed {gene_family_name}"
        
    except Exception as e:
        return f"Error processing {gene_family_name}: {str(e)}"


if __name__ == '__main__':

    # Input configuration
    gene_family_list_path = "../data/gene_family_list.csv"
    module_name = 'module3'

    # Directory paths configuration
    gene_family_results_dir = "../results/neighborhood_results"
    gene_family_download_dir = "../results/neighborhood_download"
    base_working_dir = "./"
    species_gff_data_dir = '../data/species_gff_blocks'

    # Read gene family list
    gene_family_df = pd.read_csv(gene_family_list_path, sep=',', index_col=None, header=0)
    gene_family_list = gene_family_df['genefamily'].tolist()
    
    print(f"Total gene families to process: {len(gene_family_list)}")

    # Prepare processing arguments for each gene family
    processing_args = [
        (
            gene_family,
            module_name,
            gene_family_results_dir,
            gene_family_download_dir,
            base_working_dir,
            species_gff_data_dir
        )
        for gene_family in gene_family_list
    ]

    # Configure parallel processing
    num_parallel_processes = 1

    # Process gene families in parallel
    with Pool(processes=num_parallel_processes) as pool:
        results = list(
            tqdm(
                pool.imap(process_gene_family, processing_args), 
                total=len(processing_args),
                desc="Processing gene families"
            )
        )



