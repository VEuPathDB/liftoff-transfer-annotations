import csv
from Bio import SeqIO
import argparse
import re
import pandas as pd

parser = argparse.ArgumentParser(
     formatter_class=argparse.RawDescriptionHelpFormatter,
     description='''\
Convert liftoff GFF3 outputs into Apollo friendly GFF3 files.
-------------------------------------------------------------
** This is a brute force method for one specific purpose, please be careful **

Liftoff adds several fields to the gff attributes column that are unused by Apollo.

This script will remove added attributes and convert the header to Apollo format which uses this:
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

The removed attributes are as follows:
['coverage', 'sequence_ID', 'valid_ORFs', 'extra_copy_number', 'copy_num_ID']
     ''',
     epilog="written by Helen Rebecca Davison")
parser.add_argument('-l', '--liftoff_gff', \
                    help="A gff file produced by liftoff",
                    required=True)
parser.add_argument('-f','--genome_fasta', \
                    help="A genome fasta file that corresponds to the liftoff gff",
                    required=True)
args = parser.parse_args()

genome_fasta = args.genome_fasta
liftoff_gff = args.liftoff_gff

try:
    with open(genome_fasta, "r") as f:
        alphabet = {'dna':re.compile('^[actgn]*$', re.I)}
        first = str(f.readline())
        second = str(f.readline())
        if ">" not in first:
            raise Exception("Your genome is not a fasta file")
        elif alphabet['dna'].search(second) is None:
            raise Exception("Your genome is not nucleotide sequence")
except FileNotFoundError:
    print("Genome fasta is not in the current directory")
    
try:
    with open(liftoff_gff, "r") as f:
        if "# Liftoff" not in f.read():
            print("\033[97;41m {}\033[0;0m".format("Be careful! This is not a gff produced by Liftoff and may not work as expected. Continuing, but please check your outputs!\n")) 
except FileNotFoundError:
    print("The gff is not in the current directory")

print("\033[32m {}\033[0;0m".format("Running liftoff2apollo for "+str(liftoff_gff)+"..."))



def get_annotated_contig_names(liftoff_gff):
    annotated_contigs = list(pd.read_csv(liftoff_gff, sep='\t',comment='#', header=None)[0].unique())
    return annotated_contigs

def get_contig_length(genome_fasta, liftoff_gff):
    len_dict = {}
    record_dict = SeqIO.to_dict(SeqIO.parse(genome_fasta,'fasta'))
    contig_list = get_annotated_contig_names(liftoff_gff)
    record_dict = { k: v for (k,v) in record_dict.items() if k in contig_list}
    
    for key, value in record_dict.items():
        length = len(value)
        len_dict[key] = length
    return  len_dict

def write_apollo_gff(liftoff_gff, genome_fasta):
    '''
    simple file iterator that rewrites attributes to match Apollo format.
   '''
    print("\033[33m {}\033[0;0m".format("Looking at the gff: "+str(liftoff_gff)))

    # count the number of non-comment rows
    with open(liftoff_gff, "r", encoding="utf8") as fin:
        counter = enumerate(fin)
        comment_lines = 0
        for count, value in counter:
            if any("#" in s for s in value):
                        comment_lines+=1
            num_lines = count-comment_lines
    print("\033[33m {}\033[0;0m".format("The gff has " + str(count+1) + " rows with " +str(comment_lines)+" comment lines."))
    fin.close()

    # begin rewriting the gff
    gff_name = liftoff_gff.replace(".gff","")
    with open(liftoff_gff, "r", encoding="utf8") as fin:
        with open(str(gff_name)+"_APOLLO.gff", "w", encoding="utf8", newline='') as fout:
            gff_reader = csv.reader(fin, delimiter="\t")
            gff_writer = csv.writer(fout, delimiter="\t")

            gff_writer.writerow(['##gff-version 3'])
            
            print("\033[32m {}\033[0;0m".format("Writing sequence length comments for annotated contigs in "+str(genome_fasta)+"..."))
            lengths = get_contig_length(genome_fasta, liftoff_gff)
            for key, value in lengths.items():
                gff_writer.writerow(['## sequence-region '+str(key)+' 1 '+str(value)])
            print("\033[32m {}\033[0;0m".format("Lengths obtained!"))

            print("\033[32m {}\033[0;0m".format("Fixing gff attributes..."))
            for row in gff_reader:
                if any("#" in s for s in row):
                    continue
                else:
                    seqname, source, feature, start, end, score, strand, frame, attribute = row

                    # make a dictionary for the attributes add the index number
                    attr_dict = attribute.split(';')
                    pairs = [tuple(x.split('=')) for x in attr_dict]      
                    attr_dict = dict(pairs)

                    # remove unnecessary attributes. 
                    liftoff_keys = ['coverage', 'sequence_ID', 'valid_ORFs', 'extra_copy_number', 'copy_num_ID', \
                                    'valid_ORF', 'matches_ref_protein', 'missing_start_codon', 'missing_stop_codon'\
                                    'inframe_stop_codon']
                    for key in liftoff_keys:
                        attr_dict.pop(key, None)

                    # convert the dictionary back into a gff attribute format
                    new_attribute=""
                    dict_len = len(attr_dict)
                    num_attr_added = 1
                    for key, value in attr_dict.items():
                        if num_attr_added < dict_len:
                            new_attribute += str(key) +"="+ str(value)+";"
                            num_attr_added+=1
                        elif num_attr_added >= dict_len-1:
                            new_attribute += str(key) +"="+ str(value)

                    # write the new row
                    new_row = [str(seqname), str(source), str(feature), str(start), str(end), str(score), str(strand), str(frame), str(new_attribute)]
                    gff_writer.writerow(new_row)
    fin.close()
    fout.close()
    print("\033[32m {}\033[0;0m".format("GFF attributes fixed for "+str(liftoff_gff)+". The updated gff is called: "+ str(gff_name)+"_APOLLO.gff"))
 
    return

write_apollo_gff(liftoff_gff,genome_fasta)
        
