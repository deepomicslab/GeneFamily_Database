import os
import shutil
import glob

def organize_gene_family_data(source_root, target_root, gene_family_name):
    """
    Organize gene family analysis results by archiving scattered files into categorical directories.
    
    :param source_root: Root directory of raw data (contains folders like gene_location, motif, etc.)
    :param target_root: Root directory for the organized output
    :param gene_family_name: Name of the gene family (e.g., 'Trypsin')
    """
    
    # Define common files required for every target directory
    common_files = [
        f"{gene_family_name}.merge.pep.fa",
        "Feature.all.txt"
    ]

    # Map target directories to source subfolders and filename patterns
    # Format: "Target Directory Name": [("Source Subfolder", "Filename Pattern"), ...]
    mapping = {
        "Phylogenetic Tree": [
            ("phylogenetic_trees", "filter*.merged.pep.aln.newick"),
            ("phylogenetic_trees", f"{gene_family_name}.tree.newick")
        ],
        "Sequence Identity": [
            ("sequence_alignment", f"{gene_family_name}_*.identity_matrix.txt")
        ],
        "Conserved Motif": [
            ("motif", "meme.txt")
        ],
        "Gene Collinearities": [
            ("interspecies_collinearity", "*.link.txt")
        ],
        "Duplication Events": [
            ("intraspecies_expansion", "*.tandem.link.txt")
        ],
        "Ka/Ks Ratio": [
            ("ka_ks", "*.kaks.result.txt")
        ]
    }

    print(f"Starting organization for gene family: {gene_family_name}")

    for category, sources in mapping.items():
        # 1. Create target subdirectories
        dest_dir = os.path.join(target_root, category)
        os.makedirs(dest_dir, exist_ok=True)

        # 2. Copy common files (located in the source_root directory)
        for cf in common_files:
            src_cf_path = os.path.join(source_root, cf)
            if os.path.exists(src_cf_path):
                shutil.copy2(src_cf_path, dest_dir)
            else:
                print(f"  [Warning] Common file missing: {cf}")

        # 3. Copy category-specific files
        for sub_dir, pattern in sources:
            search_path = os.path.join(source_root, sub_dir, pattern)
            matched_files = glob.glob(search_path)
            
            if not matched_files:
                print(f"  [Notice] No files matching '{pattern}' found in '{sub_dir}'")
                continue

            for fpath in matched_files:
                shutil.copy2(fpath, dest_dir)
                print(f"  Copied {os.path.basename(fpath)} -> {category}/")

    print(f"\nOrganization complete! Results saved to: {target_root}")

# --- Usage Example ---
gene_family_name = "Trypsin"
source = "/path_to/website_download"
target = "/path_to/llm_task_upload"
organize_gene_family_data(source, target, gene_family_name)


