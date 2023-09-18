#!/bin/bash

## INPUTS
target_genome=$1
reference_genome=$2
gff=$3
target_name=$(basename $target_genome | rev | cut -d "." -f2- | rev)
reference_name=$(basename $reference_genome |  rev | cut -d "." -f2- | rev)

## ERRORS

 if [ ! $# -eq 3 ]
    then
    echo " Three arguments needed, a reference genome fasta, a  "
    echo " target genome fasta, anda reference genome gff. "
    echo " Both can be relative or absolute addresses "
    exit 1
 fi

## SCRIPT
cut -f3 $gff | sort | uniq | grep -v "#" > ${reference_name}_features
liftoff -g $gff -o ${target_name}_lifted_annotations.gff3 -f ${reference_name}_features -u ${reference_name}_unmapped-f>-polish $target_genome $reference_genome