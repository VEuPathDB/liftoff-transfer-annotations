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
    print("\033[33m {}\033[0;0m".format("\n Making attributes dataframe..."))
    df[9]=attr_dicts
    df = pd.concat([df.drop([9], axis=1), df[9].apply(pd.Series)], axis=1)
    df = df[[0,1,2,3,4,5,6,7,8,'ID','Parent','gene_id']]
    print("\033[33m {}\033[0;0m".format('Calculating CDS numbers...'))
    df['CDS_number'] = pd.to_numeric(df['ID'].str.split('-CDS', expand=True)[1])
    df = df.fillna(value=np.nan)
    print("\033[33m {}\033[0;0m".format("Done!"))
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

def get_matching_genes(matching, new_gff_df, name):
    # make the dataframe
    match_stats_df = pd.DataFrame(columns=["ID"])
    print("\033[33m {}\032[0;0m".format("Cross referencing IDs..."))
    match_stats_df['ID'] = matching
    # nucleotide sequences are stored as biopython Seq class objects
    match_stats_df = match_stats_df.merge(new_gff_df.loc[new_gff_df[2]=='CDS',['ID','Parent','gene_id','CDS_number']], how='left', on = 'ID').sort_values(by='ID')
    match_stats_df.to_csv(str(name)+"-all-CDS-ids.csv", sep=',',index=None, na_rep=np.NaN)
    print("\033[33m {}\033[0;0m".format("Done! Saved as "+str(name)+"-STATS.csv"))
    return match_stats_df

def find_missing_middle(match_stats_df):
    missing_middle = []
    def checkConsecutive(l):
        return sorted(l) == list(range(min(l), max(l)+1))
    def find_missing(lst):
        return sorted(set(range(lst[0], lst[-1])) - set(lst))
    transcripts = list(match_stats_df['Parent'].unique())
    for item in tqdm(transcripts,total=len(transcripts)):
        numbers = list(map(int, match_stats_df.loc[match_stats_df['Parent']==item, 'CDS_number']))
        if checkConsecutive(numbers) == False:
            missing_middle.append(item)
    print()
    return missing_middle

def find_missing_end(match_stats_df, old_gff_df):
    old_max_cds = old_gff_df.loc[(old_gff_df[2]=="CDS"), ["Parent","ID","CDS_number"]].groupby('Parent').max()[['ID','CDS_number']].reset_index()
    match_stats_max = match_stats_df.groupby('Parent').max()['CDS_number'].reset_index()
    merged_df = old_max_cds.merge(match_stats_max, on='Parent', how='outer').fillna(np.NaN).dropna()
    missing_end = list(merged_df.loc[merged_df['CDS_number_x'] != merged_df['CDS_number_y']].dropna()['Parent'])
    return missing_end

def get_stats(old_gff_df, match_stats_df, name):
    ## finding differences
    first = (match_stats_df.sort_values(['Parent','gene_id','CDS_number']).groupby('Parent').first()['CDS_number'] != 1).reset_index()
    print("\033[32m {}\033[0;0m".format("Finding missing start cds"))
    missing_start = set(list(first.loc[first['CDS_number']==True, 'Parent']))
    print("\033[32m {}\033[0;0m".format('Finding missing middle cds'))
    missing_middle = set(list(find_missing_middle(match_stats_df)))
    print("\033[32m {}\033[0;0m".format("Finding missing end cds"))
    missing_end = set(find_missing_end(match_stats_df, old_gff_df))
    # make a venn diagram
    a = pd.DataFrame(list(missing_start), columns=['missing_start'])
    b = pd.DataFrame(list(missing_middle), columns=['missing_middle'])
    c = pd.DataFrame(list(missing_end), columns=['missing_end'])
    differences = a.join(b, how='outer').join(c, how='outer')
    differences.to_csv(name+'-missing-cds-list.csv', index=None, header=True)
    # summary to print to the console
    print()
    print("\033[32m {}\033[0;0m".format("Summary of differences between gffs:"))
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
    parser.add_argument('-o','--output', \
                        help="The name prefix to give your outputs",
                        required=True)
    args = parser.parse_args()


    old_gff = args.source_gff
    new_gff = args.lifted_gff
    name = args.output

    # ERROR HANDLING

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

    # RUN FUNCTIONS
    old_gff_df = read_in_gff(old_gff)
    new_gff_df = read_in_gff(new_gff)
    matching, absent = get_matching_ids(old_gff_df,new_gff_df)
    match_stats_df = get_matching_genes(matching, new_gff_df, name)
    print(str(len(absent))+ " ids out of a total of "+str(len(matching)+len(absent)) + " did not have a match in the new gff.")
    get_stats(old_gff_df, match_stats_df, name)
    return

main()

