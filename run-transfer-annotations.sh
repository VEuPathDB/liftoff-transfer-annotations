#!/bin/bash

# USAGE:

# bash run-transfer-annotations.sh SOURCE_GFF.gff SOURCE_GENOME.fasta NEW_GENOME.fasta output-prefix

# 1 = source gff
# 2 = source genome
# 3 = new genome
# 4 = prefix for protein check output

## ERRORS

 if [ ! $# -eq 4 ]
    then
    echo " Three arguments needed:  "
    echo " SOURCE_GFF.gff" 
    echo " SOURCE_GENOME.fasta" 
    echo " NEW_GENOME.fasta"
    echo " prefix name for the protein stats step "
    exit 1
 fi

tput setaf 6; echo "------Transfering ${1} to ${3}------"; tput sgr0

# set source for conda so environment can be changed
source ~/miniconda3/etc/profile.d/conda.sh

# run and optimise liftoff
python ~/scripts/run-and-optimise-liftoff.py -gff ${1} -old ${2} -new ${3}

tput setaf 2; echo "------Transfer step complete------"; tput sgr0
echo ""
tput setaf 6; echo "------Cleaning up transfer------"; tput sgr0

## AGAT
conda activate agatenv

lifted_gff=$(ls *lifted_annotations.gff3_polished)
outname=$(basename -s .gff3_polished $lifted_gff)

# fix phase information
agat_sp_fix_cds_phases.pl --gff ${lifted_gff} -f ${3} -o ${outname}-AGAT-phase-fixed.gff

tput setaf 2; echo "------Clean up step complete------"; tput sgr0
echo ""
tput setaf 6; echo "------Generating stats------"; tput sgr0

# get cds of source gff
agat_sp_extract_sequences.pl -g ${1} -f ${2} -t cds -p -o ${outname}-SOURCE-cds-aa.faa
# get cds of lifted gff
agat_sp_extract_sequences.pl -g ${outname}-AGAT-phase-fixed.gff -f ${3} -t cds -p -o ${outname}-LIFTED-cds-aa.faa

conda deactivate

# get missing cds
liftedgff=$(ls *-AGAT-phase-fixed.gff)
outname=$(basename -s .gff ${lifted_gff})
missingcds=${outname}-missing-CDS-list.csv
liftedaa=$(ls *-LIFTED-cds-aa.faa)
sourceaa=$(ls *-SOURCE-cds-aa.faa)
name=${4}

python ~/scripts/gff_missing_cds_finder.py --source_gff ${1} --lifted_gff ${lifted_gff} -o ${missingcds}

# run the protein change finder
python ~/scripts/gff_protein_change_finder.py -lcds ${liftedaa} -scds ${sourceaa} -nf ${3} -sf ${2} -lg ${liftedgff} -sg ${1} -m ${missingcds} -o ${name}

tput setaf 2; echo "------Annotations transferred!------"; tput sgr0

