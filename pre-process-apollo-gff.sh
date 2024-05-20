#!/bin/bash

# USAGE:

# bash pre-process-apollo-gff.sh SOURCE_GFF.gff 

# 1 = source gff

## ERRORS

 if [ ! $# -eq 1 ]
    then
    echo " One arguments needed:  "
    echo " SOURCE_GFF.gff" 
    exit 1
 fi

# set source for conda so environment can be changed
source ~/miniconda3/etc/profile.d/conda.sh
conda activate agatenv

# set output prefix
gff="$1"
outname="${gff%.*}"

## AGAT to produce a gff with unique level 3 ID and drop any without "status"
tput setaf 11; echo "Processing unique IDs for level 3 IDs..."; tput sgr0
agat_sp_manage_IDs.pl --gff $gff -p level3 -o $outname-1.gff 

grep gene $outname-1.gff  | grep status=Finished | cut -f9 | cut -d ';' -f1 | sed -s 's\ID=\\g' >> keep_list.txt
agat_sp_filter_feature_from_keep_list.pl --gff $outname-1.gff --keep_list keep_list.txt -o $outname-preprocessed.gff
rm $outname-1*

# get stats
tput setaf 11; echo "Getting stats..."; tput sgr0
agat_sq_stat_basic.pl --gff $gff -o $outname-STATS.txt

tput setaf 2; echo "Stats";printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -;tput sgr0
echo "Results are rounding to two decimal places"
echo ""
awk -F '\t' '{printf "%-40-s%-10s%-20s%-20s\n", $1,$2,$3,$4}' $outname-STATS.txt
echo ""
echo ""


# count the number of genes and how many are finished/incomplete
n_gene=$(grep gene $gff | wc -l)
n_finished_gene=$(grep gene $gff | grep -o status=Finished | wc -l)

# count the number of mrna and how many are finished/incomplete
n_mrna=$(grep mRNA $gff | wc -l)
n_finished_mrna=$(grep mRNA $gff | grep -o status=Finished | wc -l)

#print what is finished
tput setaf 2; echo "Number of finished things:" | tr ' ' -;tput sgr0
echo "$n_finished_gene genes out of $n_gene are 'Finished'"
echo "$n_finished_mrna mRNA out of $n_mrna are 'Finished'"
echo ""

# say what has been saved where
tput setaf 2; echo "New files generated";printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -;tput sgr0
echo ""
echo "processed gff:    ${outname}-preprocessed.gff"
echo "stats summary:    ${outname}-STATS.txt"







