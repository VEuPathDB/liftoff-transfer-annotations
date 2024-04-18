import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import csv
import os
import re
from tabulate import tabulate
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import substitution_matrices
from fuzzywuzzy import fuzz
import argparse
from tqdm import tqdm
tqdm.pandas()

'''
A script to annotate a gff with transfer vailidty information and protein identity
takes:
- missing CDS csv from gff_missing_cds_finder.py
- the source CDS protein fasta
- the new CDS protein fasta
- the original gff for the source genome
- the lifted gff for the new genome

Also produces a database of what the changes are and figures to illustrate changes.

'''

def read_in_gff(gff):
    df = pd.read_csv(gff, sep='\t', comment='#', header=None)
    # the columns are ordered as follows: seqname, source, feature, start, end, score, strand, phase, attribute
    # gff phase IS NOT the same as frame see: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    # the attributes comlumn must be split up so into a dictionary so it is slicable
    df[9] = df[8].str.split(';')
    attr_dicts = []
    print("\033[32m {}\033[0;0m".format("Reading "+str(gff)+"..."))
    for index, row in tqdm(df.iterrows(),total=len(df)):
        pairs = [tuple(x.split('=')) for x in row[9]]
        attr_dict = dict(pairs)
        attr_dicts += [attr_dict]
    print("\033[33m {}\033[0;0m".format("Doing some tidying..."))
    df[9]=attr_dicts
    df = pd.concat([df.drop([9], axis=1), df[9].apply(pd.Series)], axis=1)
    df = df.fillna(value=np.nan)
    print("\033[33m {}\033[0;0m".format("Done!"))
    return df

def check_valid_CDS_transfer(source_cds_aa,lifted_cds_aa,missing_cds):
    # load the fasta
    source_record_dict = SeqIO.to_dict(SeqIO.parse(source_cds_aa,"fasta"))
    lifted_record_dict = SeqIO.to_dict(SeqIO.parse(lifted_cds_aa,"fasta"))
    missing_cds_df = pd.read_csv(missing_cds)
    aligner = Align.PairwiseAligner(mode='global')
    matrix = substitution_matrices.load("BLOSUM62")
    aligner.substitution_matrix = matrix
    df = pd.DataFrame(columns=['transcript_ID','type','source_aa','lifted_aa','source_seq',\
                               	'lifted_seq','internal_stop_codons','has_start','has_end', \
                                'has_missing_cds','global_alignment_score',\
                                'percent_similarity_levenshtein','valid_transfer','proteins_match_source'])
    df['transcript_ID'] = source_record_dict.keys()
    df['type'] = 'CDS'
    print("\033[32m {}\033[0;0m".format("Getting CDS sequences and alignment scores"))
    for record in tqdm(source_record_dict, total=len(source_record_dict)):
        df.loc[df['transcript_ID']==record,'source_aa'] = str(source_record_dict[record].seq)
        if record in lifted_record_dict.keys():
            df.loc[df['transcript_ID']==record,'lifted_aa'] = str(lifted_record_dict[record].seq)
            df.loc[df['transcript_ID']==record,'global_alignment_score'] = aligner.score(source_record_dict[record].seq, lifted_record_dict[record].seq)
            # levenshtein similarity score
            df.loc[df['transcript_ID']==record,'percent_similarity_levenshtein'] = fuzz.ratio(source_record_dict[record].seq,lifted_record_dict[record].seq)
            # check if transcript has any CDS are missing
            if (record in list(missing_cds_df['missing_start'])) or (record in list(missing_cds_df['missing_middle'])) or (record in list(missing_cds_df['missing_end'])):
                df.loc[df['transcript_ID']==record, 'has_missing_cds'] = True
            else:
                df.loc[df['transcript_ID']==record, 'has_missing_cds'] = False
        # count internal stop codons in lifted aa
        df['internal_stop_codons'] = df['lifted_aa'].str[1:-2].str.count('\*')
        # check start codon and end codon in lifted aa
        df['has_start'] = df['lifted_aa'].str[0] == 'M'
        df['has_end'] = df['lifted_aa'].str[-1] == '*'
    print()
    print("\033[32m {}\033[0;0m".format("Marking valid transfers and protein matches"))
    for index, row in tqdm(df[~df['lifted_aa'].isna()].iterrows(), total = len(df[~df['lifted_aa'].isna()])):
        if (row['internal_stop_codons']!=0) or \
                (row['has_missing_cds']==True):
            df.loc[df['transcript_ID']==row['transcript_ID'], 'valid_transfer'] = False
        else:
            df.loc[df['transcript_ID']==row['transcript_ID'], 'valid_transfer'] = True
        if (row['source_aa']==row['lifted_aa']):
            df.loc[df['transcript_ID']==row['transcript_ID'], 'proteins_match_source'] = True
        elif (row['source_aa']!=row['lifted_aa']):
            df.loc[df['transcript_ID']==row['transcript_ID'], 'proteins_match_source'] = False
    print("\033[33m {}\033[0;0m".format("Done!"))
    df.loc[(df['internal_stop_codons']>0)&(~df['lifted_aa'].isna()), 'has_internal_stops'] = True
    df.loc[(df['internal_stop_codons']==0)&(~df['lifted_aa'].isna()), 'has_internal_stops'] = False
    return df
    
def check_ncRNA_and_pseudogene(source_fasta, new_fasta, source_gff, lifted_gff):
    print("\033[32m {}\033[0;0m".format("Checking ncRNA and Pseudogenes"))
    aligner = Align.PairwiseAligner(mode='global')
    matrix = substitution_matrices.load("BLOSUM62")
    aligner.substitution_matrix = matrix
    source_gff_df= read_in_gff(source_gff)
    lifted_gff_df = read_in_gff(lifted_gff)
    lifted_record_dict = SeqIO.to_dict(SeqIO.parse(new_fasta, "fasta"))
    source_record_dict = SeqIO.to_dict(SeqIO.parse(source_fasta, "fasta"))
    tr_df = pd.DataFrame(columns=['transcript_ID','type','source_seq','lifted_seq','global_alignment_score',\
                               'percent_similarity_levenshtein','valid_transfer','proteins_match_source'])
    # get transcript IDs
    def get_transcript_IDs(gff_df):
        Parent_IDs = list(gff_df.loc[gff_df[2].isin(['ncRNA_gene','pseudogene'])]['ID'])
        ncRNA_pseudo_gff_df = list(gff_df.loc[gff_df['Parent'].isin(Parent_IDs), 'ID'])
        return ncRNA_pseudo_gff_df
    lifted_ncRNA_pseudo_IDs = get_transcript_IDs(lifted_gff_df)
    source_ncRNA_pseudo_IDs = get_transcript_IDs(source_gff_df)
    tr_df['transcript_ID'] = source_ncRNA_pseudo_IDs
    # get the sequence information for lifted transcripts
    print("\033[32m {}\033[0;0m".format("Getting sequence information for lifted transcripts"))
    for index, record in tqdm(lifted_gff_df.loc[lifted_gff_df['ID'].isin(lifted_ncRNA_pseudo_IDs)].iterrows(),
                              total=len(lifted_gff_df.loc[lifted_gff_df['ID'].isin(lifted_ncRNA_pseudo_IDs)])):
        sequence = lifted_record_dict[record[0]].seq[record[3]:record[4]-1]
        seq_type = record[2]
        if record[6] == '-':
            sequence = sequence.reverse_complement()
        tr_df.loc[tr_df['transcript_ID']==record['ID'], 'lifted_seq'] = str(sequence)
        tr_df.loc[tr_df['transcript_ID']==record['ID'], 'type'] = str(seq_type)
    # get the sequence information for source transcripts
    print("\033[32m {}\033[0;0m".format("Getting sequence information for source transcripts"))
    for index, record in tqdm(source_gff_df.loc[source_gff_df['ID'].isin(source_ncRNA_pseudo_IDs)].iterrows(),
                              total=len(source_gff_df.loc[source_gff_df['ID'].isin(source_ncRNA_pseudo_IDs)])):
        sequence = source_record_dict[record[0]].seq[record[3]:record[4]-1]
        if record[6] == '-':
            sequence = sequence.reverse_complement()
        tr_df.loc[tr_df['transcript_ID']==record['ID'], 'source_seq'] = str(sequence)
    # get stats
    for index, record in tr_df.loc[tr_df['transcript_ID'].isin(lifted_ncRNA_pseudo_IDs)].iterrows():
            tr_df.loc[(tr_df['transcript_ID']==record['transcript_ID']),'global_alignment_score'] = aligner.score(Seq(record['lifted_seq']).translate(), Seq(record['source_seq']).translate())
            # levenshtein similarity score
            tr_df.loc[(tr_df['transcript_ID']==record['transcript_ID']),'percent_similarity_levenshtein'] = fuzz.ratio(Seq(record['lifted_seq']).translate(), Seq(record['source_seq']).translate())
            tr_df.loc[(tr_df['transcript_ID']==record['transcript_ID']),'proteins_match_source'] = Seq(record['lifted_seq']).translate() == Seq(record['source_seq']).translate()
    tr_df['valid_transfer'] =  tr_df['proteins_match_source']
    return tr_df, lifted_gff_df

def add_attribute_to_gff(lifted_gff_df, transcript_changes_df, lifted_gff, name, outdir):
    merged = lifted_gff_df.merge(transcript_changes_df[['transcript_ID','has_missing_cds','has_internal_stops','valid_transfer','proteins_match_source']], left_on='ID', right_on='transcript_ID', how='left')
    # attributes to keep
    # remove ebi biotypes from all rows
    merged = merged[[0,1,2,3,4,5,6,7,8,'ID','description','Parent','Note','gene_id', 
                     'protein_source_id','Name','has_missing_cds',
                     'has_internal_stops','valid_transfer','proteins_match_source']]
    merged.fillna(value=np.nan, inplace=True)
    attribute_names = ['ID','description','Parent','Note','gene_id', 
                     'protein_source_id','Name','has_missing_cds',
                     'has_internal_stops','valid_transfer','proteins_match_source']
    new_attributes = []
    for index, row in tqdm(merged.iterrows(),total=len(merged)):
        attr = ''
        for name in attribute_names:
            if pd.isna(row[name])==False:
                attr = attr + name + "="+str(row[name])+';'
        attr = attr[:-1]
        new_attributes.append(attr)
    merged['col_9'] = new_attributes
    # write new gff
    gff_name = str(os.path.basename(lifted_gff).replace(".gff","_changes_marked.gff"))
    out_name = os.path.join(outdir, gff_name)
    print("\033[97;46m {}\033[0;0m".format("Writing new gff to:"))
    print(tabulate([['Annotated gff', out_name]], 
                    headers=['Type', 'File name']))
    print('')
    comments=[]
    with open(lifted_gff, "r", encoding="utf8") as fin:
         gff_reader = csv.reader(fin, delimiter="\t")
         for row in gff_reader:
            if any("#" in s for s in row):
                    comments.append(row)
    f = open(str(out_name), 'a', newline='')
    for i in comments:
        f.write(str(i[0])+"\n")
    merged[[0,1,2,3,4,5,6,7,'col_9']].to_csv(f, sep='\t', index=None, header=None)
    f.close()
    return

def plot_and_save(df, outdir,name):
    # save the dataframe
    csvname = name+'_protein-change-stats.csv'
    summarypng = name+'_CDS-summary.png'
    summarypng2 = name+'_ncRNA-pseudogene-summary.png'
    failurepng = name+'_CDS-failure-summary.png'
    print("\033[97;46m {}\033[0;0m".format("Saving plots and summaries of liftoff successes and failures:"))
    print(tabulate([['Stats', csvname], 
                    ['Transfer successes and failures', summarypng],
                    ['Breakdown of failures', failurepng],
                    ['ncRNA and pseudogenes transfers', summarypng2]],
                    headers=['', 'File name']))
    print(os.path.join(outdir, csvname))
    df.to_csv(os.path.join(outdir, csvname),header=True, index=None)
    # plot overall valid CDS transfers and protein matches
    hue_order=[True, False]
    fig, ax = plt.subplots(figsize=(6,5))
    sns.countplot(df[df['type']=='CDS'], x='valid_transfer', hue='proteins_match_source',
                  hue_order=hue_order, order=[True,False], palette='bone')
    plt.title('CDS transfer rates')
    plt.xlabel('Valid transfer')
    plt.ylabel('Count')
    plt.legend(title='Proteins match')
    plt.savefig(os.path.join(outdir, summarypng),dpi=150)
    # plot a breakdown of invalid CDS transfers
    fig, ax = plt.subplots(figsize=(6,5))
    sns.countplot(df.loc[(df['valid_transfer']==False)&(df['type']=='CDS')][['has_missing_cds','has_internal_stops','proteins_match_source']].melt()\
                  .replace({'has_missing_cds':'Has missing cds','has_internal_stops':'Has internal stops','proteins_match_source':'Proteins match source'}),\
                   hue='value', hue_order=hue_order, x='variable', palette='bone')
    plt.title('CDS transfer failure breakdown for '+str(len(df.loc[(df['valid_transfer']==False)&(df['type']=='CDS')]))+' failures')
    plt.ylabel('Count')
    plt.xlabel('Failure type')
    plt.legend(title='Value')
    plt.savefig(os.path.join(outdir, failurepng),dpi=150)
    # plot valid transfers for ncRNA and pseudogenes
    fig, ax = plt.subplots(figsize=(6,6))
    sns.countplot(df[df['type']!='CDS'], x='type', hue='valid_transfer', hue_order=hue_order, palette='bone')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.title('ncRNA and pseudogene transfer rates')
    plt.xlabel('Type')
    plt.ylabel('Count')
    plt.legend(title='Valid transfer')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, summarypng2),dpi=150)
    return

def main():
    '''
    
    '''
    parser = argparse.ArgumentParser(
     formatter_class=argparse.RawDescriptionHelpFormatter,
     description='''\
Adds record of liftoff errors to gff attributes
-------------------------------------------------------------
A script to annotate a gff with transfer vailidty information and protein identity
takes:
- A list of missing CDS output from gff_missing_cds_finder.py
- the source genome nucleotide fasta
- the new genome nucleotide fasta
- the original gff for the source genome
- a gff of annotations lifted from the source genome to the new genome

Also produces a database of what the changes are and figures to illustrate changes.

     ''',
     epilog="written by Helen Rebecca Davison") 
    parser.add_argument('-lcds', '--lifted_cds_aa', \
                    help="the new CDS proteins retrieved with AGAT",
                    required=True)
    parser.add_argument('-scds','--source_cds_aa', \
                    help="the source CDS proteins retrieved with AGAT",
                    required=True)
    parser.add_argument('-nf', '--new_fasta', \
                    help="the new genome fasta associated with the lifted gff",
                    required=True)
    parser.add_argument('-sf','--source_fasta', \
                    help="the source genome fasta associated with the source gff and lifted gff",
                    required=True)
    parser.add_argument('-lg','--lifted_gff', \
                    help="the gff associated with the new genome produced by liftoff and phase fixed with AGAT to be annotated with information",
                    required=True)
    parser.add_argument('-sg','--source_gff', \
                    help="the gff associated with the source genome and lifted gff",
                    required=True)   
    parser.add_argument('-m','--missing_cds_lists', \
                    help="csv of missing cds with the headings 'missing_start', 'missing_middle', 'missing_end'",
                    required=True)
    parser.add_argument('-o','--output', \
                    help="output name",
                    required=True)
     
    args = parser.parse_args()

    # define input variables
    lifted_cds = args.lifted_cds_aa
    source_cds= args.source_cds_aa
    missing_cds = args.missing_cds_lists
    
    new_fasta = args.new_fasta
    source_fasta = args.source_fasta

    lifted_gff = args.lifted_gff
    source_gff = args.source_gff

    name = args.output

    # ERROR HANDLING
    try:
        def read(fasta):
            record_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
            sequence = record_dict[list(record_dict.keys())[0]].seq
            return str(sequence)
        def validate(seq, alphabet='dna'):
            alphabets = {'dna': re.compile('^[acgtn]*$', re.I), 
             'protein': re.compile('^[acdefghiklmnpqrstvwy*]*$', re.I)}
            if alphabets[alphabet].search(seq) is not None:
                return True
            else:
                return False
        if validate(read(source_cds),'protein') == False:
            raise Exception("Your source cds is not protein sequence")
        if validate(read(lifted_cds),'protein') == False:
            raise Exception("Your lifted cds is not protein sequence")
        if validate(read(source_fasta),'dna') == False:
            raise Exception("Your source fasta is not nucleotide sequence")
        if validate(read(new_fasta),'dna') == False:
            raise Exception("Your new fasta is not nucleotide sequence")
    except FileNotFoundError:
        print("One or more sequence files are not in the in the current directory")

    try:
        def is_fasta(filename):
            with open(filename, "r") as handle:
                fasta = SeqIO.parse(handle, "fasta")
                return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file
        if is_fasta(source_cds) == False:
            raise Exception("Your source cds is not fasta")
        if is_fasta(lifted_cds) == False:
            raise Exception("Your lifted cds is not fasta")
        if is_fasta(source_fasta) == False:
            raise Exception("Your source fasta is not fasta")
        if is_fasta(new_fasta) == False:
            raise Exception("Your new fasta is not fasta")
    except FileNotFoundError:
        print("One or more sequence files are not in the in the current directory")
        
    try:
        with open(source_gff, "r") as f:
            if "#" not in str(f.readline()):
                raise Exception("No header lines found, check that your source gff file is definitely a gff")
        gff_id = [line for line in open(source_gff) if line[:1] != '#'][0].split()[0]
    except FileNotFoundError:
        print("Source gff is not in the current directory")

    try:
        with open(lifted_gff, "r") as f:
            if "#" not in str(f.readline()):
                raise Exception("No header lines found, check that your new gff file is definitely a gff")
        gff_id = [line for line in open(lifted_gff) if line[:1] != '#'][0].split()[0]
    except FileNotFoundError:
        print("New gff is not in the current directory")

    try:
        gff_id = [line for line in open(source_gff) if line[:1] != '#'][0].split()[0]
        with open(source_fasta, "r") as f:
            found=False
            if any(line for line in f if gff_id in line):
                found=True
            if found == False:
                raise Exception("Contig IDs do not match between the source gff and source genome")
    except FileNotFoundError:
        print("Your source gff file is not in the current directory")
    
    try:
        gff_id = [line for line in open(lifted_gff) if line[:1] != '#'][0].split()[0]
        with open(new_fasta, "r") as f:
            found=False
            if any(line for line in f if gff_id in line):
                found=True
            if found == False:
                raise Exception("Contig IDs do not match between the lifted gff and new genome")
    except FileNotFoundError:
        print("Your new gff file is not in the current directory")
    
    # define output dir
    filename = name+"-cds-content-check"
    outdir = os.path.join(os.getcwd(),filename)
    os.makedirs(outdir, exist_ok=False)

    print("\033[97;46m {}\033[0;0m".format("Begining to find and annotate changes "))
    print("\033[33m {}\033[0;0m".format("Input files being used: \n"))
    print(tabulate([['Source genome', source_fasta],
                    ['New genome', new_fasta],
                    ['Source gff', source_gff],
                    ['Lifted gff', lifted_gff],
                    ['Source CDS protein', source_cds],
                    ['Lifted CDS protein', lifted_cds]
                    ], 
                    headers=['Type', 'File name']))
    print('')

    # get protein records
    transcript_changes_df = check_valid_CDS_transfer(source_cds,lifted_cds,missing_cds)
    tr_df, lifted_gff_df= check_ncRNA_and_pseudogene(source_fasta, new_fasta, source_gff, lifted_gff)
    transcript_changes_df = transcript_changes_df.merge(tr_df, how='outer')
    # output gff
    add_attribute_to_gff(lifted_gff_df, transcript_changes_df, lifted_gff, name, outdir)
    # output plots
    plot_and_save(transcript_changes_df, outdir, name)

    return()

main()