from Bio import SeqIO
from Bio import Align
import pandas as pd
import numpy as np
import argparse
import os
from fuzzywuzzy import fuzz
import re
from time import time
from venny4py.venny4py import *
  
  
def timer_func(func):
    # This function shows the execution time of 
    # the function object passed
    def wrap_func(*args, **kwargs):
        t1 = time()
        result = func(*args, **kwargs)
        t2 = time()
        print(f'Time taken: {(t2-t1):.4f}s')
        return result
    return wrap_func

def progress_bar(progress, total):
    percent = (100*progress/float(total))
    bar = "█"*int(percent) + "-" * (100 - int(percent))
    print(f"\r|{bar}| {percent:.2f}%", end="\r")

# function to load the gff into a pandas dataframe and unpack the attributes column to access the unique IDs
@timer_func
def read_in_gff(gff):
    df = pd.read_csv(gff, sep='\t', comment='#', header=None)
    # the columns are ordered as follows: seqname, source, feature, start, end, score, strand, phase, attribute
    # gff phase IS NOT the same as frame see: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    # the attributes comlumn must be split up so into a dictionary so it is slicable
    df[9] = df[8].str.split(';')
    attr_dicts = []
    print("\033[32m {}\033[0;0m".format("Reading "+str(gff)+"..."))
    for index, row in df.iterrows():
        pairs = [tuple(x.split('=')) for x in row[9]]
        attr_dict = dict(pairs)
        attr_dicts += [attr_dict]
        progress_bar(index+1,len(df))
    print("\033[32m {}\033[0;0m".format("\n Doing some tidying..."))
    df[9]=attr_dicts
    df = pd.concat([df.drop([9], axis=1), df[9].apply(pd.Series)], axis=1)
    df = df[[0,1,2,3,4,5,6,7,8,'ID','Parent','gene_id']]
    df['CDS_number'] = pd.to_numeric(df['ID'].str.split('-CDS', expand=True)[1])
    df = df.fillna(value=np.nan)
    print("\033[32m {}\033[0;0m".format("Done!\n"))
    return df

@timer_func
def get_matching_ids(old_gff_df,new_gff_df):
    # get a list of matching IDs
    matching = []
    absent = []
    old_cds = old_gff_df.loc[old_gff_df[2]=='CDS']
    new_cds = new_gff_df.loc[new_gff_df[2]=='CDS']
    print("\033[32m {}\033[0;0m".format("Finding matching gene CDS exon IDs. This might take a minute..."))
    for index, id in old_cds['ID'].items():
        if id in list(new_cds['ID']):
            matching.append(id)
        elif id not in list(new_cds['ID']):
            absent.append(id)
        progress_bar(index+1,len(old_gff_df))
    print("\n")
    return matching, absent # absent will be used for metrics

def get_gene_sequence(df, ID, record_dict):
    # gets sequence for a gff record
    start = df.loc[df['ID']==ID][3].item()
    end = df.loc[df['ID']==ID][4].item()
    seqname = df.loc[df['ID']==ID][0].item()
    if df.loc[df['ID']==ID][6].item() == "+":
        sequence = record_dict[seqname].seq[int(start)-1:int(end)]
        strand = "+"
    elif df.loc[df['ID']==ID][6].item() == "-":
        sequence = record_dict[seqname].seq[int(start)-1:int(end)].reverse_complement()
        strand = "-"
    return sequence, strand

def nucleotide_match(old_sequence, new_sequence):
    aligner = Align.PairwiseAligner()
    old_sequence = old_sequence.upper()
    new_sequence = new_sequence.upper()
    old_length = len(old_sequence)
    new_length = len(new_sequence)
    length_diff = old_length - new_length
    levenshtein_ratio = fuzz.ratio(old_sequence,new_sequence)
    score = aligner.score(old_sequence, new_sequence)
    nucleotide_match_stats = [old_length, new_length, length_diff, levenshtein_ratio, score]
    return nucleotide_match_stats

def protein_match(old_sequence, new_sequence, old_phase, new_phase):
    aligner = Align.PairwiseAligner()
    matrix = Align.substitution_matrices.load('BLOSUM62')
    aligner.substitution_matrix = matrix
    old_protein = old_sequence[int(old_phase):].translate(to_stop=True)
    new_protein = new_sequence[int(new_phase):].translate(to_stop=True)
    old_length = len(old_protein)
    new_length = len(new_protein)
    length_diff = old_length - new_length
    levenshtein_ratio = fuzz.ratio(old_protein,new_protein)
    if (old_length>0) and (new_length>0):
        score = aligner.score(old_protein, new_protein)
    else:
        score = np.nan
    protein_match_stats = [old_length, new_length, length_diff, levenshtein_ratio, score]
    return protein_match_stats

@timer_func
def get_matching_genes(matching, old_gff_df, new_gff_df, old_fasta, new_fasta):
    # Find percent match between nucleotide sequences with matching gene ids
    # between new and old gffs
    print("\033[32m {}\033[0;0m".format("Calculating gene similarity metrics for "+str(len(matching))+" IDs"))
    old_record_dict = SeqIO.to_dict(SeqIO.parse(old_fasta,'fasta'))
    new_record_dict = SeqIO.to_dict(SeqIO.parse(new_fasta,'fasta'))
    match_stats = [["ID", "strand_match","nucl_old_length", "nucl_new_length", "nucl_length_diff", "nucl_levenshtein_ratio", "nucl_score",\
                   "aa_old_length", "aa_new_length", "aa_length_diff", "aa_levenshtein_ratio", "aa_score"]]
    for index, ID in enumerate(matching):
        old_sequence, old_strand = get_gene_sequence(old_gff_df, ID, old_record_dict)
        new_sequence, new_strand = get_gene_sequence(new_gff_df, ID, new_record_dict)
        strand_match = (old_strand==new_strand)
        nucleotide_match_stats = nucleotide_match(old_sequence, new_sequence)
        old_phase = old_gff_df.loc[old_gff_df['ID']==ID][7].item()
        new_phase = new_gff_df.loc[new_gff_df['ID']==ID][7].item()
        protein_match_stats = protein_match(old_sequence, new_sequence, old_phase, new_phase)
        match_stats.append([ID, strand_match] + nucleotide_match_stats + protein_match_stats)
        progress_bar(index+1,len(matching))
    print()
    return match_stats

def find_missing_middle(match_stats):
    missing_middle = []
    print('Finding missing middle cds')
    for index, item in match_stats['gene_id'].items():
        cds = match_stats.loc[match_stats['gene_id']==item].sort_values(by='CDS_number')['CDS_number'].reset_index(drop=True)  
        if cds[0] == 1:
            expected=1
        elif cds[0] != 1:
            expected = int(cds[0])
        for i in list(cds):
            if expected == i: 
                expected+=1
            elif (expected != i) & (item not in missing_middle):
                missing_middle.append(item)
        progress_bar(index+1,len(match_stats))
    print()
    return missing_middle

def find_missing_end(match_stats, old_gff_df):
    print("Finding missing ends cds")
    old_max_cds = old_gff_df.loc[(old_gff_df[2]=="CDS"), ["Parent","ID","CDS_number"]].groupby('Parent').max()[['ID','CDS_number']].reset_index()
    old_max_cds['gene_id'] = old_max_cds.loc[:,'ID'].str.split('-CDS', expand=True)[0]
    match_stats_max = match_stats.groupby('gene_id').max()['CDS_number'].reset_index()
    merged_df = old_max_cds.merge(match_stats_max, on='gene_id', how='outer').fillna(np.NaN).dropna()
    missing_end = list(merged_df.loc[merged_df['CDS_number_x'] != merged_df['CDS_number_y']].dropna()['gene_id'])
    return missing_end

@timer_func
def get_stats(old_gff_df, match_stats, name):
    df = pd.DataFrame(match_stats[1:], columns=match_stats[0])
    df['CDS_number'] = pd.to_numeric(df.loc[:,'ID'].str.split('-CDS', expand=True)[1])
    df['gene_id'] = df.loc[:,'ID'].str.split('-CDS', expand=True)[0]
    df.to_csv(str(name)+"-STATS.csv", sep=',',index=None, na_rep=np.NaN)
    ## calculations
    first = (df.sort_values(['gene_id','CDS_number']).groupby('gene_id').first()['CDS_number'] != 1).reset_index()
    ## the stats to keep
    print("Finding changes in sequence length")
    changed_nucl = set(list(df.loc[df['nucl_length_diff']>0, 'gene_id']))
    print("Finding missing start cds")
    missing_start = set(list(first.loc[first['CDS_number']==True, 'gene_id']))
    missing_middle = set(find_missing_middle(df))
    missing_end = set(find_missing_end(df, old_gff_df))
    print(str(len(list(df.loc[df['nucl_length_diff']>0, 'ID'])))+" out of "+str(len(df)) + " CDS have changed sequence length")
    print(str(len(missing_start))+ " out of "+str(len(df['gene_id'].unique())) + " protein coding sequences have a missing start codon")
    print(str(len(missing_middle))+" out of "+str(len(df['gene_id'].unique())) + " protein coding sequences have a missing middle codon")
    print(str(len(missing_end))+" out of "+str(len(df['gene_id'].unique())) + " protein coding sequences have a missing end codon")
    differences = {"changed_nucl":changed_nucl, "missing_start":missing_start, "missing_middle":missing_middle, "missing_end":missing_end}
    np.save(name+'-missing_cds_lists.npy',differences)
    venny4py(sets=differences, out = name+"_venn", ext='svg')
    return differences, df

@timer_func
def main():
    parser = argparse.ArgumentParser(
        description="examine the similarity of an old and new gff for one genome",
        epilog="written by Helen Davison")
    parser.add_argument('-gff1', '--old_gff', \
                        help="The original gff for the genome",
                        required=True)
    parser.add_argument('-gff2', '--new_gff', \
                        help="The new gff for the genome to compare against the original one",
                        required=True)
    parser.add_argument('-f1','--old_genome', \
                        help="A genome fasta file that corresponds to the liftoff gff",
                        required=True)
    parser.add_argument('-f2','--new_genome', \
                        help="A genome fasta file that corresponds to the liftoff gff",
                        required=True)
    parser.add_argument('-o','--output', \
                        help="The name prefix to give your outputs",
                        required=True)
    args = parser.parse_args()


    old_genome = args.old_genome
    new_genome = args.new_genome
    old_gff = args.old_gff
    new_gff = args.new_gff
    name = args.output

    # ERROR HANDLING
    try:
        with open(old_genome, "r") as f:
            alphabet = {'dna':re.compile('^[actgn]*$', re.I)}
            first = str(f.readline())
            second = str(f.readline())
            if ">" not in first:
                raise Exception("Your old genome is not a fasta file")
            elif alphabet['dna'].search(second) is None:
                raise Exception("Your old genome is not nucleotide sequence")
    except FileNotFoundError:
        print("Old genome fasta is not in the current directory")

    try:
        with open(new_genome, "r") as f:
            alphabet = {'dna':re.compile('^[actgn]*$', re.I)}
            first = str(f.readline())
            second = str(f.readline())
            if ">" not in first:
                raise Exception("Your new genome is not a fasta file")
            elif alphabet['dna'].search(second) is None:
                raise Exception("Your new genome is not nucleotide sequence")
    except FileNotFoundError:
        print("New genome fasta is not in the current directory")

    try:
        with open(old_gff, "r") as f:
            if "#" not in str(f.readline()):
                raise Exception("No header lines found, check that your old gff file is definitely a gff")
        gff_id = [line for line in open(old_gff) if line[:1] != '#'][0].split()[0]
    except FileNotFoundError:
        print("Old gff is not in the current directory")

    try:
        with open(new_gff, "r") as f:
            if "#" not in str(f.readline()):
                raise Exception("No header lines found, check that your new gff file is definitely a gff")
        gff_id = [line for line in open(old_gff) if line[:1] != '#'][0].split()[0]
    except FileNotFoundError:
        print("New gff is not in the current directory")

    try:
        gff_id = [line for line in open(old_gff) if line[:1] != '#'][0].split()[0]
        with open(old_genome, "r") as f:
            found=False
            if any(line for line in f if gff_id in line):
                found=True
            if found == False:
                raise Exception("Contig IDs do not match between the old gff and old genome")
    except FileNotFoundError:
        print("Your old gff file is not in the current directory")
    
    try:
        gff_id = [line for line in open(new_gff) if line[:1] != '#'][0].split()[0]
        with open(new_genome, "r") as f:
            found=False
            if any(line for line in f if gff_id in line):
                found=True
            if found == False:
                raise Exception("Contig IDs do not match between the new gff and new genome")
    except FileNotFoundError:
        print("Your new gff file is not in the current directory")

    # RUN FUNCTIONS
    old_gff_df = read_in_gff(old_gff)
    new_gff_df = read_in_gff(new_gff)
    matching, absent = get_matching_ids(old_gff_df,new_gff_df)
    match_stats = get_matching_genes(matching, old_gff_df, new_gff_df, old_genome, new_genome)
    print(str(len(absent))+ " ids out of a total of "+str(len(matching)+len(absent)) + " did not have a match in the new gff.")
    differences, stats_df = get_stats(old_gff_df, match_stats, name)

    #for item in list(stats_df['gene_id']):
    #    if item not in list(new_gff_df.loc[new_gff_df[2]=='CDS']['Parent']):
    #        missing.append(item)
    #        print(item)
    return
main()

