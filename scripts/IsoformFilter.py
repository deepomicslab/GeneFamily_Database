import pandas as pd
import os, json, csv
from tqdm import tqdm
import glob

def get_files_with_prefix_suffix(directory, prefix='', suffix=''):
    matching_files = []
    for filename in os.listdir(directory):
        if filename.startswith(prefix) and filename.endswith(suffix):
            file = os.path.join(directory, filename)
            matching_files.append(file)
    return matching_files

def get_seq_lengths(fasta_path):
    seq_lengths = {}
    current_seq = ''
    current_name = ''
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    seq_lengths[current_name.split(' ')[0].split('.')[0]] = len(current_seq)
                current_name = line[1:]
                current_seq = ''
            else:
                current_seq += line
        if current_name:
            seq_lengths[current_name.split(' ')[0].split('.')[0]] = len(current_seq)
            
    return seq_lengths

def parse_feature_file(file_path):
    """
    Parse feature file using pandas and group proteins by gene symbol

    Args:
        file_path: path to Feature.all.txt file
    Returns:
        dict: {gene_symbol: [protein_id1, protein_id2, ...]}
    """
    # Read file
    df = pd.read_csv(file_path, sep='\t', index_col=0, header=0)

    # Process columns
    df['Protein_ID'] = df['Protein_ID'].str.split('.').str[0]
    df['mRNA_ID'] = df['mRNA_ID'].str.split('.').str[0]
    df['Gene_symbol'] = df['Gene_symbol'].fillna('Unknown').str.upper()

    # Group by gene symbol and aggregate protein IDs into lists
    gene_proteins = df.groupby('Gene_symbol')['Protein_ID'].agg(lambda x: list(set(x))).to_dict()

    return gene_proteins

def select_longest_protein(gene_proteins_dict, fasta_lengths):
    """
    Args:
        gene_proteins_dict: {gene_symbol: [protein_id1, protein_id2, ...]}
        fasta_lengths: {protein_id: length}
    Returns:
        {gene_symbol: longest_protein_id}
    """
    longest_proteins = {}

    for gene, proteins in gene_proteins_dict.items():
        # Get lengths for all proteins of this gene
        protein_lengths = {}
        for protein in proteins:
            if protein in fasta_lengths:
                protein_lengths[protein] = fasta_lengths[protein]

        # Skip if no proteins found in fasta
        if not protein_lengths:
            continue

        # Get protein with max length
        longest_protein = max(protein_lengths.items(), key=lambda x: x[1])[0]
        longest_proteins[gene] = longest_protein

    return longest_proteins


def select_longest_protein_by_species(feature_df, fasta_lengths):
    """
    Select longest protein for each gene symbol within each species

    Args:
        feature_df: DataFrame containing Feature.all.txt data
        fasta_lengths: dict of protein lengths {protein_id_prefix: length}
    Returns:
        DataFrame with only longest proteins selected, keeping the original Protein_ID and mRNA_ID values (with decimals)
    """
    feature_df['protein_length'] = feature_df['Protein_ID'].str.split('.').str[0].map(fasta_lengths)
    known_df = feature_df[feature_df['Gene_ID'] != 'UNKNOWN']

    longest_proteins = (
        known_df.sort_values('protein_length', ascending=False)
                .groupby(['Gene_ID', 'Species'])
                .first()
                .reset_index()
    )

    result_df = longest_proteins.drop('protein_length', axis=1)

    return result_df

def process_genefamily(feature_file_path, fasta_lengths):
    """
    Process a gene family file to select longest proteins

    Args:
        feature_file_path: path to Feature.all.txt
        fasta_lengths: dict of protein lengths {protein_id_prefix: length}
    Returns:
        Filtered DataFrame with original Protein_ID and mRNA_ID (including decimals)
    """

    df = pd.read_csv(feature_file_path, sep='\t', index_col=None, header=0)
    df['Gene_symbol'] = df['Gene_symbol'].fillna('Unknown').str.upper()
    filtered_df = select_longest_protein_by_species(df, fasta_lengths)

    return filtered_df

def load_json(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data

def save_to_json(data, file_path):
    with open(file_path, mode='w', encoding='utf-8') as json_file:
        json.dump(data, json_file, ensure_ascii=False, indent=4)
    
if __name__ == "__main__":
    print('Start ----- ')
    
    # Base directory configurations
    species_peptide_dir = '/path/to/database/pep'  # Directory containing species peptide fasta files
    peptide_length_output_dir = '../data/species_pep_length'  # Directory to save peptide length results
    gene_family_list_path = "../data/gene_family_list.csv"  # Path to gene family list CSV
    gene_family_source_dir = '../data'  # Source directory containing gene family data
    gene_family_filtered_dir = '../results'  # Output directory for filtered gene families
    all_species_peptide_length_path = '../data/species_pep_length/all.species.pep.length.json'  # Path to consolidated peptide length data

    # Read gene family list from CSV
    gene_family_df = pd.read_csv(gene_family_list_path, index_col=None, header=0)
    gene_family_list = gene_family_df['genefamily'].tolist()

    # Process peptide sequences and calculate lengths
    peptide_fasta_files = get_files_with_prefix_suffix(species_peptide_dir, prefix='', suffix='.fa')
    all_peptide_lengths = {}
    
    # Calculate peptide lengths for each species
    for idx, peptide_file in tqdm(enumerate(peptide_fasta_files), total=len(peptide_fasta_files)):
        species_name = os.path.basename(peptide_file).replace('.pep.fa', '')
        species_peptide_lengths = get_seq_lengths(peptide_file)
        all_peptide_lengths.update(species_peptide_lengths)
    
    # Save consolidated peptide lengths
    peptide_length_output_path = os.path.join(peptide_length_output_dir, "all.species.pep.length.json")
    save_to_json(all_peptide_lengths, peptide_length_output_path)

    # Load peptide length data for filtering
    peptide_length_data = load_json(all_species_peptide_length_path)

    # Process each gene family
    for idx, gene_family_name in tqdm(enumerate(gene_family_list), total=len(gene_family_list)):
        # Define input and output paths for current gene family
        feature_input_path = os.path.join(gene_family_source_dir, gene_family_name, 'Feature.all.txt')
        
        # Skip if source file doesn't exist
        if not os.path.exists(feature_input_path):
            continue
            
        # Read and process feature data
        feature_df = pd.read_csv(feature_input_path, sep='\t')
        original_columns = feature_df.columns.tolist()
        
        # Filter gene family data
        filtered_feature_df = process_genefamily(feature_input_path, peptide_length_data)
        filtered_feature_df = filtered_feature_df.reindex(columns=original_columns)
        
        # Save filtered results
        feature_output_dir = os.path.join(gene_family_filtered_dir, gene_family_name)
        os.makedirs(feature_output_dir, exist_ok=True)
        
        feature_output_path = os.path.join(feature_output_dir, 'Feature.all.txt')
        filtered_feature_df.to_csv(feature_output_path, sep='\t', index=False)



