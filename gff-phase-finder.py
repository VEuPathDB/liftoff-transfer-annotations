import pandas as pd
import numpy as np
import csv
import argparse

def progress_bar(progress, total):
    percent = (100*progress/float(total))
    bar = "█"*int(percent) + "-" * (100 - int(percent))
    print(f"\r|{bar}| {percent:.2f}%", end="\r")

def read_in_gff(liftoff_gff):
    # save the commentlines from the gff
    comments=[]
    with open(liftoff_gff, "r", encoding="utf8") as fin:
         gff_reader = csv.reader(fin, delimiter="\t")
         for row in gff_reader:
            if any("#" in s for s in row):
                    comments.append(row)
        
    df = pd.read_csv(liftoff_gff, sep='\t', comment='#', header=None)
    # the columns are ordered as follows: seqname, source, feature, start, end, score, strand, phase, attribute
    # gff phase IS NOT the same as frame see: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    # the attributes comlumn must be split up so into a dictionary so it is slicable
    # any protein coding genes that are partial mappings are discarded
    df[9] = df[8].str.split(';')
    attr_dicts = []
    print("\033[32m {}\033[0;0m".format("Reading "+str(liftoff_gff)+". This might take a minute..."))
    for index, row in df.iterrows():
        pairs = [tuple(x.split('=')) for x in row[9]]
        attr_dict = dict(pairs)
        attr_dicts += [attr_dict]
        progress_bar(index+1,len(df))
    print("\033[32m {}\033[0;0m".format("Loaded! Doing some tidying..."))
    df[9]=attr_dicts
    df = pd.concat([df.drop([9], axis=1), df[9].apply(pd.Series)], axis=1)
    #valid_genes = df[(df[2]=='protein_coding_gene')&(df['partial_mapping']!='True')].ID.unique()
    parents = df[df[2]=="CDS"].Parent.unique()
    df = df[[0,1,2,3,4,5,6,7,8,'ID','Parent','gene_id']]
    df.loc[:,'CDS_number']  = pd.to_numeric(df.loc[:,'ID'].str.split('-CDS', expand=True)[1])
    df = df.fillna(value=np.nan)
    return df, parents, comments #,valid_genes


def find_phase(df, parents, comments):
    print("\033[32m {}\033[0;0m".format("Finding phase. This might take a hot minute..."))
    for index, name in enumerate(list(parents)):
        start = []
        end = []
        remainder_of_bases = 0
        if  df.loc[(df['Parent']==name) & (df[2]=='CDS')].sort_values(by='CDS_number').iloc[0,-1] == 1:
            df.loc[(df['Parent']==name) & (df[2]=='CDS') & (df['CDS_number']==1), 7] = 3
            start = df.loc[(df['Parent']==name) & (df[2]=='CDS') & (df['CDS_number']==1), 3].item()
            end = df.loc[(df['Parent']==name) & (df[2]=='CDS') & (df['CDS_number']==1), 4].item()
            remainder_of_bases = 3-((start-end-1)%3)
        else:
            first = df.loc[(df['Parent']==name) & (df[2]=='CDS')].sort_values(by='CDS_number').iloc[0,-1]
            df.loc[(df['Parent']==name) & (df[2]=='CDS') & (df['CDS_number']==first), 7] = 3
            start = df.loc[(df['Parent']==name) & (df[2]=='CDS') & (df['CDS_number']==first), 3].item()
            end = df.loc[(df['Parent']==name) & (df[2]=='CDS') & (df['CDS_number']==first), 4].item()
            remainder_of_bases = 3-((start-end-1)%3)
        cds = df.loc[(df['Parent']==name) & (df[2]=='CDS')].sort_values(by='CDS_number')[1:]['CDS_number'] 
        for number in cds:
            df.loc[(df['Parent']==name) & (df['CDS_number']==number) & (df[2]=='CDS'), 7] = remainder_of_bases
            start = df.loc[(df['Parent']==name) & (df[2]=='CDS') & (df['CDS_number']==number) & (df[2]=='CDS'), 3].item()
            end = df.loc[(df['Parent']==name) & (df[2]=='CDS') & (df['CDS_number']==number) & (df[2]=='CDS'), 4].item()
            remainder_of_bases = 3-((start-end-1-remainder_of_bases)%3)
        progress_bar(index+1,len(parents))
    for index, row in df.iterrows():
        if row[7] == 3:
                df.iloc[index, 7] = 0
        elif row[7] == 2:
                df.iloc[index, 7] = 1
        elif row[7] == 1:
                df.iloc[index, 7] = 2
    print()
    return



# need to -1 from all of column for the negative row (check the positive), so that it is phase rather than the number of needed bases

def write_gff(liftoff_gff, df, comments):
    f = open(liftoff_gff+"-PHASE-CORRECTED.gff", 'a', newline='')
    for i in comments:
        f.write(str(i[0])+"\n")
    df[[0,1,2,3,4,5,6,7,8]].to_csv(f, sep='\t', index=None, header=None)
    f.close()
    return

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
     description='''
        Calculates GFF phase information for CDS which lack it. 
        -------------------------------------------------------
        This was made for liftoff gff outputs which lack this information.
        GFF phase is not the same as reading frame! Please see the documentation:
        https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md''',
        epilog="written by Helen R. Davison")
    parser.add_argument('-gff', '--liftoff_gff', dest="liftoff_gff",\
                        type=argparse.FileType('r'), help="The gff without phase information (made for liftoff)",
                        required=True)
    args = parser.parse_args()
    # Name the inputs
    liftoff_gff = args.liftoff_gff.name


    # Run the scripts
    df, parents, comments = read_in_gff(liftoff_gff)
    print("\033[32m {}\033[0;0m".format("Done!"))
    find_phase(df, parents, comments)
    write_gff(liftoff_gff, df, comments)
    return

main()
