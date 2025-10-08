#!/bin/bash
# init variables
genome=/data/zyz/genefamily/database/genome
cds=/data/zyz/genefamily/database/cds
pep=/data/zyz/genefamily/database/pep
script=/data/zyz/genefamily/database/script
gff=/data/zyz/genefamily/database/gff
blastdb_pep=/data/zyz/genefamily/database/blastdb_pep
chromosome=/data/zyz/genefamily/database/gene_density
while getopts ":q:s:e:p:i:a:n:t:m:" o; do
    case "${o}" in
        q)
            queryID="$OPTARG"
            ;;
        s)
            species_list="$OPTARG"
            ;;
        e)
            evalue="$OPTARG"
            ;;
        p)
            projectID="$OPTARG"
            ;;
        i)
            minw="$OPTARG"
            ;;
        a)
            maxw="$OPTARG"
            ;;
        n)
            nmotifs="$OPTARG"
            ;;
        t)
            type="$OPTARG"
            ;;
        m)
            speciesname="$OPTARG"
            ;;
        
    esac
done


echo  "queryID = ${queryID}"

workspace=./$projectID

if [ ! -d "${workspace}" ]; then
    echo "no exsist  $workspace}"
    exit 1
fi

tmp=$workspace/tmp
mkdir $workspace/proteinfasta
proteinfasta=$workspace/proteinfasta

feature=$workspace/Feature.all.txt

mkdir $workspace/gene_structure

awk '{print $1}' $feature | sed '1d' | uniq > $workspace/species_list.txt 
cat $workspace/species_list.txt | while read id
do
  grep ${id} $feature | awk '{print $3}' \
  > $workspace/gene_structure/${id}_transIDlist.txt #convert proteinID to transcriptID
  perl $script/get_gene_exon_from_gff.pl \
  -in1 $workspace/gene_structure/${id}_transIDlist.txt \
  -in2 $gff/${id}.gff.fa \
  -out $workspace/gene_structure/${id}_gene_exon_info.gff
  python $script/gff2csv.py \
  -i $workspace/gene_structure/ \
  -o $workspace/
done




