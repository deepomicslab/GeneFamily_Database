import os
import re
import csv
import json
import pickle
from collections import defaultdict
from tqdm import tqdm
import numpy as np
import pandas as pd

def get_files_with_prefix_suffix(directory, prefix='', suffix=''):
    matching_files = []
    for filename in os.listdir(directory):
        if filename.startswith(prefix) and filename.endswith(suffix):
            file = os.path.join(directory, filename)
            matching_files.append(file)
    return matching_files

def save_to_pickle(data, output_path):
    with open(output_path, 'wb') as f:
        pickle.dump(data, f)

def load_from_pickle(file_path):
    with open(file_path, 'rb') as f:
        return pickle.load(f)

def extract_transcript_id(attributes):
    """
    Extract transcript ID from attributes string.
    First try to get transcript info from 'ID=transcript:xxx',
    if not found, extract from 'Parent=transcript:xxx'.
    
    Args:
        attributes (str): Attributes string containing transcript information
        
    Returns:
        str or None: Transcript ID if found, None otherwise
    """
    transcript_id = None
    id_match = re.search(r'ID=transcript:([^;]+)', attributes)
    if id_match:
        transcript_id = id_match.group(1)
    else:
        parent_match = re.search(r'Parent=transcript:([^;]+)', attributes)
        if parent_match:
            transcript_id = parent_match.group(1)
    return transcript_id

def read_gff_file(file_path):
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 8:
                record = {
                    'chromosome': fields[0],
                    'source': fields[1],
                    'type': fields[2],
                    'start': fields[3],
                    'end': fields[4],
                    'score': fields[5],
                    'strand': fields[6],
                    'phase': fields[7],
                    'attributes': fields[8] if len(fields) > 8 else ''
                }
                data.append(record)
    return data


def parse_attributes(attr_str):
    attrs = {}
    if ';' in attr_str:
        items = [x.strip() for x in attr_str.split(';') if x.strip()]
    else:
        items = [attr_str]

    for item in items:
        if 'ID=transcript:' in item:
            attrs['transcriptID'] = item.replace("ID=transcript:","")
        if 'Parent=gene:' in item:
            attrs['geneID'] = item.replace("Parent=gene:","")
        if 'Name=' in item:
            attrs['geneName'] = item.replace("Name=","")
        if 'Parent=transcript:' in item:
            attrs['transcriptID'] = item.replace("Parent=transcript:", "")
    if 'transcriptID' not in attrs:
        attrs['transcriptID'] = ''
    if 'geneID' not in attrs:
        attrs['geneID'] = ''
    if 'geneName' not in attrs:
        attrs['geneName'] = ''
    return attrs

def process_gene_exon_gff(gff_path):

    exon_info_data = read_gff_file(gff_path)
    
    df = pd.DataFrame(exon_info_data)
    df['transcript_id'] = df['attributes'].apply(extract_transcript_id)
    
    df = df[df['transcript_id'].notna()]
    
    grouped = df.groupby('transcript_id')
    
    results = []
    filename = os.path.basename(gff_path)
    
    for transcript_id, group in grouped:
        current_transcript = []
        for _, block in group.iterrows():
            attrs = parse_attributes(block['attributes'])
            result_row = {
                'chr': block['chromosome'],
                'type': block['type'],
                'start': int(block['start']),
                'end': int(block['end']),
                'strand': block['strand'],
                'transcriptID': attrs['transcriptID'],
                'geneID': attrs['geneID'],
                'geneName': attrs['geneName'],
                'species': filename.replace('.gff.fa', '')
            }
            current_transcript.append(result_row)
        if current_transcript:
            results.append(current_transcript)
    
    return results

def read_gff_for_donwload_file(file_path, file_header):
    """
    Read GFF file and parse it into a list of dictionaries using custom headers.
    
    Args:
        file_path (str): Path to the GFF file
        file_header (list): List of field names, e.g. ["chromosome", "source", "type", ...]
    
    Returns:
        list: List of dictionaries containing data from each row
    """
    data = []
    
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            # Skip comment lines (starting with '#') and empty lines
            if line.startswith('#') or not line.strip():
                continue
            
            # Split fields by tab
            fields = line.strip().split('\t')
            
            # Ensure number of fields matches header length
            if len(fields) != len(file_header):
                print(f"Warning: Field count mismatch (expected {len(file_header)}, got {len(fields)}), skipping line: {line}")
                continue
            
            # Combine fields with headers into dictionary
            row_dict = dict(zip(file_header, fields))
            data.append(row_dict)
    
    return data

if __name__ == "__main__":

    file_header = [ 
        "chromosome",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    ]

    gff_database_dir = '/data/genefamily/database/gff'
    species_gff_results_dir = '../data/species_gff_blocks'

    # Get GFF files list
    gff_files_list = get_files_with_prefix_suffix(gff_database_dir, prefix='', suffix='.gff.fa')

    for gff_file in tqdm(gff_files_list):
        species_name = os.path.basename(gff_file).replace('.gff.fa', '')
        tmp_mrna_blocks = process_gene_exon_gff(gff_file)
        tmp_mrna_blocks_items = []
        for tmp_mrna_block in tmp_mrna_blocks:
            tmp_mrna_blocks_items.extend(tmp_mrna_block)
        
        # Save results using path joining
        output_file = os.path.join(species_gff_results_dir, 
                                f'{species_name}_species_gff_block_items.pkl')
        save_to_pickle(tmp_mrna_blocks_items, output_file)


    for gff_file in tqdm(gff_files_list):
        species_name = os.path.basename(gff_file).replace('.gff.fa', '')
        tmp_mrna_blocks = read_gff_for_donwload_file(gff_file, file_header)
        # Save results using path joining
        output_file = os.path.join(species_gff_results_dir, 
                                f'{species_name}_species_gff_download.pkl')
        save_to_pickle(tmp_mrna_blocks, output_file)

    
        
