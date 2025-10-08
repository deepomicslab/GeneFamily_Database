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

workspace=/mnt/genefamily/genefamily/module3/$projectID

if [ ! -d "${workspace}" ]; then
    echo "no exsist  $workspace}"
    exit 1
fi

tmp=$workspace/tmp
mkdir $workspace/proteinfasta
proteinfasta=$workspace/proteinfasta

feature=$workspace/Feature.all.txt

mkdir $workspace/gene_structure


awk '{print $1}' $feature | sed '1d' | sort -u > $workspace/species_list.txt 
sed -i 's/\r$//' $workspace/species_list.txt

cat $workspace/species_list.txt | while read id
do
  grep ${id} $feature | awk '{print $3}' \
  > $workspace/gene_structure/${id}_transIDlist.txt #convert proteinID to transcriptID
done

cat $workspace/species_list.txt | while read id
do
  sed -i "s/\..*//g" $workspace/gene_structure/${id}_transIDlist.txt
  perl $script/get_gene_location.pl \
  -in1 $workspace/gene_structure/${id}_transIDlist.txt \
  -in2 $gff/${id}.gff.fa \
  -out $workspace/${id}_mrna_location.txt;
  cp $chromosome/${id}.txt $workspace/${id}_chrlength.txt; #各个物种染色体长度信息
  cp $chromosome/${id}gene_density.txt $workspace/${id}.gene_density.txt; #各个物种染色体密度
done

rm -rf $tmp
rm -rf $proteinfasta

