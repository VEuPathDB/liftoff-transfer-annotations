#!/bin/bash

## USAGE:
# script genome.fasta [corresponding-gff.gff]

## DESCRIPTION
# This script provides stats for genomes: 
# Number of contigs, assembly length, longest contig, N50, N90, GC%

# It can also optionally count gff features, but you might need to 
# change lines 67 onwards to match the contents of your gff 
# (e.g. "ncRNA_gene" might need to be "ncRNA")
# This was made to deal with VEuPathDB gff3 syntax

## INPUTS
genome=$1
gff=$2

name=$(basename $genome | rev | cut -d "." -f2- | rev)

mkdir temp;
#get contig lengths, ordered from large to small; 1. remove newlines in sequences, 2. remove fasta headers
# 3. get contig sizes (line lengths) and order them from large to small.
cat $genome | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }'| sed '/^>/ d'| awk '{ print length($0) }' | sort -gr > temp/contig_lengths.txt;

#number of contigs
Y=$(cat temp/contig_lengths.txt | wc -l);

#sum of contig_lengths
X=$(paste -sd+ temp/contig_lengths.txt | bc);

# cumulative contig lengths

awk 'BEGIN {sum=0} {sum= sum+$0; print sum}' temp/contig_lengths.txt > temp/contig_lengths_cum.txt;

# get cumulative contig contributions (%) to the entire assembly

awk -v var=$X 'BEGIN {FS=OFS=","} {for (i=1; i<=NF; i++) $i/=var;}1' temp/contig_lengths_cum.txt > temp/cum_perc.txt;

# join results

paste temp/contig_lengths.txt temp/cum_perc.txt > temp/matrix.txt;

# get N50, N90, largest scaffold/contig

N50=$(awk '$2 >= 0.50' temp/matrix.txt |head -1| awk '{ print $1}');
N90=$(awk '$2 >= 0.90' temp/matrix.txt |head -1| awk '{ print $1}');
large_contig=$(head -1 temp/contig_lengths.txt);
rm -r temp;

# calculate GC
# Find all Gs and Cs in sequence array and store to another array
lengths=`grep -v ">" $genome | wc | awk '{print $3-$1}'`
C=`grep -v ">" $genome | awk -F"[Cc]" '{print NF-1}'| awk '{sum += $1} END {print sum}'`
G=`grep -v ">" $genome| awk -F"[Gg]" '{print NF-1}'| awk '{sum += $1} END {print sum}'`
GC=`expr $C + $G`
percent=`echo "$GC / $lengths * 100" | bc -l`

## Summarise genome results:

echo "assembly  $name" > ${name}-summary.txt
echo "number of contigs/scaffolds   $Y" >> ${name}-summary.txt
echo "assembly size $X" >> ${name}-summary.txt
echo "largest contig/scaffold   $large_contig" >> ${name}-summary.txt
echo "N50   $N50" >> ${name}-summary.txt
echo "N90   $N90" >> ${name}-summary.txt
echo "GC%   $percent" >> ${name}-summary.txt


## if there is a gff specified also count gff features
if [ ${gff:+1} ]
then
    protein_coding_gene=$(cut -f3 $gff | grep -x protein_coding_gene | wc -l)
    if (( ${protein_coding_gene} == 0 ))
    then
        protein_coding_gene=$( grep protein_coding $gff | wc -l)
    fi
    CDS=$(cut -f3 $gff | grep -x CDS | wc -l)
    ncRNA=$(cut -f3 $gff | grep -x ncRNA_gene | wc -l)
    if (( ${ncRNA} == 0 ))
    then
        ncRNA=$(cut -f3 $gff | grep -x ncRNA | wc -l)
    fi
    pseudogene=$(cut -f3 $gff | grep -x pseudogene | wc -l)

    echo "protein coding gene count   $protein_coding_gene" >> ${name}-summary.txt
    echo "CDS count   $CDS" >> ${name}-summary.txt
    echo "ncRNA gene count   $ncRNA" >> ${name}-summary.txt
    echo "pseudogene count   $pseudogene" >> ${name}-summary.txt
fi

cat ${name}-summary.txt
echo
echo "Output saved to ${name}-summary.txt"