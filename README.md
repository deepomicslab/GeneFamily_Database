# GeneFamily: A Comprehensive Mammalian Gene Family Database with Extensive Annotation and Interactive Visualization

## Abstract
Gene families underlie genetic innovation and phenotypic diversification. However, the identification, analysis, and interpretation of gene families depend on custom bioinformatics and visualization workflows that are mainly unattainable for non-expert users. To overcome the learning curve associated with analytical tools and the time costs involved in data collection and processing, we propose GeneFamily (https://gh.deepomics.org/), a comprehensive mammalian gene family database that incorporates extensive annotations and interactive visualizations. Leveraging whole-genome data from 138 mammalian species and Pfam-A Hidden Markov models, GeneFamily encompasses 2,036 gene families. Users can explore gene family distributions, phylogenetic trees, member diversities, sequence identities, structural features, conserved motifs, chromosomal localizations, neighboring regions, gene collinearities, duplication events, and adaptive signatures via an intuitive web interface. GeneFamily reduces technical barriers and accelerates research into gene family evolution, functional diversification, and comparative genomics across mammals.

## Environment
- Anaconda
- Python 3.9

## Dependencies
- pandas >= 2.0.3
- tqdm >= 4.66.5
- scipy >= 1.10.1
- numpy >= 1.24.3
- perl >= 5.10.0

## 1. IsoformFilter.py
### Description
IsoformFilter.py is used to filter protein isoforms to ensure that only one representative protein sequence is retained for each gene. During evolutionary analysis, a gene may produce multiple protein isoforms due to alternative splicing or post-translational modifications. These isoforms do not represent independent evolutionary units, so a representative sequence needs to be selected for analysis.

### Input Files
- Protein sequence files (.pep format) from Ensembl mammalian database
- Initial gene family search results: Feature1.all.txt

### Output Files
- Filtered results file: Feature.all.txt

### Processing Logic
- Group protein sequences by gene ID and species information
- Sort protein sequences within each group by length in descending order
- Select the longest protein sequence in each group as the gene representative
- Filter out remaining isoform sequences

## 2. Gene Family Member Counter
### Function Description
This module is used to count and organize gene family member information, classifying protein sequences that have been filtered by IsoformFilter.py into corresponding Gene Symbols (as gene members of the gene family), and establishing mapping relationships for gene family members.

### Input Files
- Gene family data file processed by IsoformFilter.py: Feature1.all.txt

### Output Files
```json
{
    "gene_family_name": ["member1", "member2", ...]
}
```

### Processing Logic
- Read filtered protein sequence data
- Match corresponding gene symbols for each sequence
- Mark sequences with missing gene symbols as "UNKNOWN"
- Remove duplicates from each gene family member list
- Generate mapping dictionary from gene family to member list

### Purpose
- Establish clear gene family member relationship diagrams
- Facilitate subsequent gene family analysis and evolutionary research
- Ensure uniqueness and traceability of each gene family member

## 3. GeneFamilyDiversityExtractor.py
### Function Description
This script analyzes the distribution of gene families across different species, generates gene family diversity data, and provides foundational data support for evolutionary analysis and inter-species comparison.

### Input Files
- Gene family data file processed by IsoformFilter.py: Feature1.all.txt

### Output Files
Gene family distribution statistics across species, including:
- Number of members of specific gene families in each species
- Species distribution patterns of gene families
- Differences in gene family member numbers between species

## 4. Gene Structure Extraction Pipeline
### Overview
This pipeline contains two scripts for extracting and organizing gene structure information of gene family members, providing data support for comparative analysis of gene structures.

### 4.1 GeneStructureExtractorStep1.sh
#### Input Files
- Filtered gene family data: Feature1.all.txt
- Ensembl mammalian database GFF files

#### Output Files
- Gene structure information categorized by species for the gene family

#### Functions
- Extract structure information of target genes from GFF files
- Parse exon and intron location information
- Extract transcript structure features
- Collect gene structure-related annotation information

### 4.2 GeneStructureExtractorStep2.py
#### Input Files
- Gene structure information output from Step1
- Gene family member information and species distribution information

#### Output Files
- Gene structure information organized by gene family members
- Gene structure information organized by gene family species distribution

#### Functions
- Group gene structure information by species
- Classify by gene family members
- Generate structured data files

## 5. GeneLocationExtractor.sh
### Function Description
This script extracts genomic location information of gene family members and combines chromosome gene density data to provide data support for gene distribution pattern analysis.

### Input Files
- Filtered gene family data: Feature1.all.txt
- Ensembl mammalian database GFF files
- Chromosome length information files
- Chromosome gene density data

### Output Files
- Chromosome position information of gene transcripts
- Chromosome gene density distribution data

### Processing Flow
#### Parameter Processing
- Process input parameters (project ID, species list, etc.)
- Validate working directory

#### Data Preprocessing
- Create necessary working directories
- Extract species list
- Generate transcript ID list

#### Location Information Extraction
- Parse gene location information from GFF files
- Extract chromosome coordinates of transcripts

#### Chromosome Data Processing
- Copy chromosome length information
- Integrate gene density data


## 6. Gene Neighborhood Analysis Pipeline
### Script Components
#### 1. NeighboringGeneExtractorStep1.py
##### Function Description
- Process GFF format genome annotation files
- Extract gene and transcript structure information
- Establish gene position index for each species

##### Input Files
- Ensembl mammalian database GFF files
- Filtered gene family transcript structure information

##### Output Files
- Species-specific gene block data (.pkl format)
- Gene structure information cache files

#### 2. NeighboringGeneExtractorStep2.py
##### Function Description
- Analyze genomic neighboring regions (±1Mbp range)
- Extract genomic structure around target genes
- Organize data by species and gene family members

##### Input Files
- Gene structure information data corresponding to gene family
- Species-specific gene block data (.pkl format)
- Gene structure information cache files

##### Output Files
- Gene neighborhood structure information organized by gene family members
- Gene neighborhood structure information organized by gene family species distribution

### Processing Flow
#### Data Preprocessing
- Read gene family information
- Parse GFF annotation files
- Establish gene position index

#### Neighboring Region Analysis
- Determine 1Mbp range search window
- Extract gene information in target regions
- Collect structure features of neighboring genes

#### Data Organization
- Group by species
- Group by gene family members
- Generate JSON format output

## 7. GeneFamilyInterspeciesCollinearityStatistics.py
### Function Description
This script performs statistical analysis based on gene family collinearity analysis results between different species, revealing gene family evolutionary dynamics and functional differentiation patterns through systematic classification and quantitative assessment of gene symbol conservation and variation.

### Input Files
- Gene family data after isoform filtering: Feature1.all.txt
- Inter-species collinearity analysis results (from MCScanX analysis)

### Output Files
#### Conserved Gene Symbol Collinearity Relationships
- {gene_family_name}.links.flow.txt
- Records collinear gene pairs with identical gene symbols

#### Variable Gene Symbol Collinearity Relationships
- {gene_family_name}.links.flow.unknown.txt
- Records collinear gene pairs with unknown or different gene symbols

#### Species-Level Collinearity Statistics
- {gene_family_name}.species.links.flow.txt
- Statistics of gene family member collinearity distribution in each species

### Analysis Flow
#### 1. Data Loading and Preprocessing
##### 1.1 Gene Family Data Reading
- Read Feature.all.txt file
- Standardize mRNA_ID (remove version numbers)
- Process Gene_symbol (convert to uppercase, fill NA with "UNKNOWN")
- Establish mRNA_ID to Gene_symbol mapping dictionary

##### 1.2 Collinearity Data Processing
- Filter .link.txt files
- Filter tandem repeat (tandem.link.txt) files
- Extract inter-species collinearity relationships

#### 2. Gene Symbol Relationship Analysis
##### 2.1 Gene Pair Classification
- Same gene symbol pairs (same_genemembers)
- Unknown gene symbol pairs (diff_unknown_genemembers)
- Different gene symbol pairs (diff_genesymbol_genemembers)

##### 2.2 Species Relationship Statistics
- Establish inter-species connection counts (sp_counts)
- Calculate gene symbol frequency in each species (sp_genecounts)
- Generate species gene member matrix (sp_genemembers)

### Statistical Analysis Focus
#### Gene Symbol Conservation Analysis
- Calculate percentage of completely conserved gene symbols between species
- Identify highly conserved gene family members
- Evaluate degree of gene function conservation

#### Collinearity Pattern Statistics
- Quantify collinearity relationship strength between different species
- Analyze distribution characteristics of collinear blocks
- Study conservation level of gene arrangement