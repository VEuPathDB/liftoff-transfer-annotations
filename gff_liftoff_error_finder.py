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
from tqdm import tqdm
tqdm.pandas()
  
  
def timer_func(func):
    # This function shows the execution time of 
    # the function object passed
    def wrap_func(*args, **kwargs):
        t1 = time()
        result = func(*args, **kwargs)
        t2 = time()
        print(f'Time taken: {(t2-t1):.4f}s')
        print()
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
    print("\033[33m {}\033[0;0m".format("Doing some tidying..."))
    df[9]=attr_dicts
    df = pd.concat([df.drop([9], axis=1), df[9].apply(pd.Series)], axis=1)
    df = df[[0,1,2,3,4,5,6,7,8,'ID','Parent','gene_id']]
    df['CDS_number'] = pd.to_numeric(df['ID'].str.split('-CDS', expand=True)[1])
    df = df.fillna(value=np.nan)
    print("\033[33m {}\033[0;0m".format("Done!\n"))
    return df

@timer_func
def get_matching_ids(old_gff_df,new_gff_df):
    # get a list of matching IDs
    old_cds = set(old_gff_df.loc[old_gff_df[2]=='CDS', 'ID'])
    new_cds = set(new_gff_df.loc[new_gff_df[2]=='CDS', 'ID'])
    print("\033[32m {}\033[0;0m".format("Finding matching gene CDS exon IDs."))
    matching = list(old_cds.intersection(new_cds))
    absent = list(old_cds.symmetric_difference(new_cds))
    print("\033[33m {}\033[0;0m".format("Done!"))
    return matching, absent

def get_gene_sequence(df, ID, record_dict):
    # gets sequence for a gff record
    start = df.loc[df['ID']==ID][3].item()
    end = df.loc[df['ID']==ID][4].item()
    seqname = df.loc[df['ID']==ID][0].item()
    if df.loc[df['ID']==ID][6].item() == "+":
        sequence = record_dict[seqname].seq[int(start)-1:int(end)]
    elif df.loc[df['ID']==ID][6].item() == "-":
        sequence = record_dict[seqname].seq[int(start)-1:int(end)].reverse_complement()
    return sequence

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

def get_matching_genes(matching, old_gff_df, new_gff_df, old_fasta, new_fasta, name):
    # Find percent match between nucleotide sequences with matching gene ids
    # between new and old gffs
    print("\033[32m {}\033[0;0m".format("Calculating gene similarity metrics for "+str(len(matching))+" IDs..."))
    old_record_dict = SeqIO.to_dict(SeqIO.parse(old_fasta,'fasta'))
    new_record_dict = SeqIO.to_dict(SeqIO.parse(new_fasta,'fasta'))
    # make the dataframe
    match_stats_df = pd.DataFrame(columns=["ID", "Parent", "strand_match",'old_strand', 'old_phase', 'old_nucl_seq', \
                                           'new_strand', 'new_phase', 'new_nucl_seq', "nucl_old_length", "nucl_new_length",\
                                            "nucl_len_diff", "nucl_levenshtein_ratio", "aa_old_length", "aa_new_length",\
                                            "aa_len_diff", "aa_levenshtein_ratio"])
    print("\033[33m {}\033[0;0m".format("Getting IDs..."))
    match_stats_df['ID'] = matching
    # nucleotide sequences are stored as biopython Seq class objects
    print("\033[33m {}\033[0;0m".format("Getting sequence info..."))
    all = match_stats_df['ID'].progress_apply(lambda x: (get_gene_sequence(new_gff_df, x, new_record_dict), \
                                                         get_gene_sequence(old_gff_df, x, old_record_dict),
                                                         new_gff_df.loc[new_gff_df['ID']==x]['Parent'].item(),\
                                                         old_gff_df.loc[old_gff_df['ID']==x][7].item(),\
                                                         new_gff_df.loc[new_gff_df['ID']==x][7].item(),\
                                                         old_gff_df.loc[old_gff_df['ID']==x][6].item(),\
                                                         new_gff_df.loc[new_gff_df['ID']==x][6].item()))
    match_stats_df['new_nucl_seq'], match_stats_df['old_nucl_seq'], match_stats_df['Parent'], match_stats_df['old_phase'],match_stats_df['new_phase'], match_stats_df['old_strand'],match_stats_df['new_strand']  = zip(*all)
    match_stats_df['strand_match'] = match_stats_df['old_strand'] == match_stats_df['new_strand']
    print("\033[33m {}\033[0;0m".format("Getting nucleotide stats..."))
    match_stats_df["nucl_levenshtein_ratio"] =match_stats_df[['new_nucl_seq','old_nucl_seq']].progress_apply(lambda x: fuzz.ratio(x[0],x[1]), axis=1)
    match_stats_df["nucl_old_length"] = match_stats_df['old_nucl_seq'].str.len()
    match_stats_df["nucl_new_length"] = match_stats_df['new_nucl_seq'].str.len()
    match_stats_df["nucl_len_diff"] = match_stats_df['nucl_old_length'] - match_stats_df["nucl_new_length"]
    print("\033[33m {}\033[0;0m".format("Getting protein stats..."))
    aa_lev_ratio =[]
    old_aa_len = []
    new_aa_len = []
    for index, row in tqdm(match_stats_df.iterrows()):
        old_protein = row['old_nucl_seq'][int(row['old_phase']):].translate(to_stop=True)
        new_protein = row['new_nucl_seq'][int(row['new_phase']):].translate(to_stop=True)
        if (len(old_protein)>0) and (len(new_protein)>0):
            aa_lev_ratio.append(fuzz.ratio(old_protein, new_protein))
        else:
            aa_lev_ratio.append(np.nan)
        old_aa_len.append(len(old_protein))
        new_aa_len.append(len(new_protein))
    match_stats_df["aa_old_length"] = old_aa_len
    match_stats_df["aa_new_length"] = new_aa_len
    match_stats_df["aa_levenshtein_ratio"] = aa_lev_ratio
    match_stats_df["aa_len_diff"] = (match_stats_df['aa_old_length'] - match_stats_df["aa_new_length"])
    print("\033[33m {}\033[0;0m".format('Calculating CDS numbers...'))
    match_stats_df['CDS_number'] = pd.to_numeric(match_stats_df.loc[:,'ID'].str.split('-CDS', expand=True)[1])
    match_stats_df['gene_id'] = match_stats_df.loc[:,'ID'].str.split('-CDS', expand=True)[0]
    match_stats_df.to_csv(str(name)+"-STATS.csv", sep=',',index=None, na_rep=np.NaN)
    print("\033[33m {}\033[0;0m".format("Done! Saved as "+str(name)+"-STATS.csv"))
    return match_stats_df

def find_missing_middle(match_stats_df):
    missing_middle = []
    for index, item in tqdm(match_stats_df['Parent'].items()):
        cds = match_stats_df.loc[match_stats_df['Parent']==item].sort_values(by='CDS_number')['CDS_number'].reset_index(drop=True)  
        if cds[0] == 1:
            expected=1
        elif cds[0] != 1:
            expected = int(cds[0])
        for i in list(cds):
            if expected == i: 
                expected+=1
            elif (expected != i) & (item not in missing_middle):
                missing_middle.append(item)
    print()
    return missing_middle

def find_missing_end(match_stats_df, old_gff_df):
    old_max_cds = old_gff_df.loc[(old_gff_df[2]=="CDS"), ["Parent","ID","CDS_number"]].groupby('Parent').max()[['ID','CDS_number']].reset_index()
    old_max_cds['gene_id'] = old_max_cds.loc[:,'ID'].str.split('-CDS', expand=True)[0]
    match_stats_max = match_stats_df.groupby('gene_id').max()['CDS_number'].reset_index()
    merged_df = old_max_cds.merge(match_stats_max, on='gene_id', how='outer').fillna(np.NaN).dropna()
    missing_end = list(merged_df.loc[merged_df['CDS_number_x'] != merged_df['CDS_number_y']].dropna()['Parent'])
    return missing_end

def get_stats(old_gff_df, match_stats_df, name):
    ## finding differences
    first = (match_stats_df.sort_values(['Parent','CDS_number']).groupby('Parent').first()['CDS_number'] != 1).reset_index()
    print("\033[32m {}\033[0;0m".format("Finding changes in sequence length"))
    changed_aa = set(list(match_stats_df.loc[match_stats_df['aa_levenshtein_ratio']!=100, 'Parent']))
    print("\033[32m {}\033[0;0m".format("Finding missing start cds"))
    missing_start = set(list(first.loc[first['CDS_number']==True, 'Parent']))
    print("\033[32m {}\033[0;0m".format('Finding missing middle cds'))
    missing_middle = set(list(find_missing_middle(match_stats_df)))
    print("\033[32m {}\033[0;0m".format("Finding missing end cds"))
    missing_end = set(find_missing_end(match_stats_df, old_gff_df))
    # make a venn diagram
    differences = {"changed_aa":changed_aa, "missing_start":missing_start, "missing_middle":missing_middle, "missing_end":missing_end}
    venny4py(sets=differences, ext='svg')
    a = pd.DataFrame(differences['changed_aa'], columns=['changed_aa'])
    b = pd.DataFrame(differences['missing_start'], columns=['missing_start'])
    c = pd.DataFrame(differences['missing_middle'], columns=['missing_middle'])
    d = pd.DataFrame(differences['missing_end'], columns=['missing_end'])
    differences = a.join(b, how='outer').join(c, how='outer').join(d, how='outer')
    differences.to_csv(name+'-missing-cds-list.csv', index=None, header=True)
    # summary to print to the console
    print()
    print("\033[32m {}\033[0;0m".format("Summary of differences between gffs:"))
    print(str(len(list(match_stats_df.loc[match_stats_df['nucl_len_diff']>0, 'ID'])))+" out of "+str(len(match_stats_df)) + " CDS have changed sequence length")
    print(str(len(list(match_stats_df.loc[match_stats_df['aa_levenshtein_ratio']!=100, 'ID'])))+" out of "+str(len(match_stats_df)) + " CDS have changed amino acid sequence")
    print(str(len(missing_start))+ " out of "+str(len(match_stats_df['Parent'].unique())) + " transcripts have a missing start codon")
    print(str(len(missing_middle))+" out of "+str(len(match_stats_df['Parent'].unique())) + " transcripts have a missing middle codon")
    print(str(len(missing_end))+" out of "+str(len(match_stats_df['Parent'].unique())) + " transcripts have a missing end codon")
    return differences

@timer_func
def main():
    parser = argparse.ArgumentParser(
        description="examine the similarity of an old and new gff for one genome",
        epilog="written by Helen Davison")
    parser.add_argument('-sg', '--source_gff', \
                        help="The source gff for the genome",
                        required=True)
    parser.add_argument('-lg', '--lifted_gff', \
                        help="The lifted gff for the genome to compare against the original one",
                        required=True)
    parser.add_argument('-sf','--source_genome', \
                        help="The source genome fasta file that corresponds to the liftoff and source gff",
                        required=True)
    parser.add_argument('-nf','--lifted_genome', \
                        help="The new genome fasta file that corresponds to the liftoff gff",
                        required=True)
    parser.add_argument('-o','--output', \
                        help="The name prefix to give your outputs",
                        required=True)
    args = parser.parse_args()


    old_genome = args.source_genome
    new_genome = args.lifted_genome
    old_gff = args.source_gff
    new_gff = args.lifted_gff
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
        gff_id = [line for line in open(new_gff) if line[:1] != '#'][0].split()[0]
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
    match_stats_df = get_matching_genes(matching, old_gff_df, new_gff_df, old_genome, new_genome, name)
    print(str(len(absent))+ " ids out of a total of "+str(len(matching)+len(absent)) + " did not have a match in the new gff.")
    differences = get_stats(old_gff_df, match_stats_df, name)

    #for item in list(stats_df['gene_id']):
    #    if item not in list(new_gff_df.loc[new_gff_df[2]=='CDS']['Parent']):
    #        missing.append(item)
    #        print(item)
    return
main()

