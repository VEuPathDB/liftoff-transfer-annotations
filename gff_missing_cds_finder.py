import pandas as pd
import numpy as np
import argparse
from time import time
from tabulate import tabulate
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
    if 'gene_id' in df.columns:
        df = df[[0,1,2,3,4,5,6,7,8,'ID','Parent','gene_id']]
    elif 'owner' in df.columns:
        df = df[[0,1,2,3,4,5,6,7,8,'ID','Parent','Name']]
        # gff_type = 'Apollo'
    print("\033[33m {}\033[0;0m".format("Done!"))
    return df


def find_missing_cds(old_gff_df, new_gff_df, name):
    print("\033[33m {}\033[0;0m".format("Cross referencing IDs and finding missing CDS..."))
    # get locations
    old_cds = old_gff_df.loc[old_gff_df[2]=='CDS', ['ID', 'Parent']].sort_values(['Parent','ID']).reset_index().drop(columns='index')
    new_cds = new_gff_df.loc[new_gff_df[2]=='CDS', ['ID','Parent']].sort_values(['Parent','ID']).reset_index().drop(columns='index')
    # find and label mismatching rows, then sort so coordinates are in series
    matching = old_cds.merge(new_cds,on=['ID', 'Parent'], how='left',indicator='Exist')
    matching['Exist'] = np.where(matching.Exist == 'both', True, False)
    # find missing start, middle, and end:
    parents_with_missing = matching.loc[matching['Exist']==False,'Parent'].unique()
    if len(parents_with_missing) == 0:
        matching['status']='present'
    else:
        for i in parents_with_missing:
            n_cds = range(0, len(matching.loc[matching['Parent']==i]))
            status = []
            if matching.loc[matching['Parent']==i,'Exist'].eq(False).all():
                status=['missing_all_cds']*(len(n_cds))
            else:
                if matching.loc[matching['Parent']==i].iloc[0,4] == False:
                    status.append('missing_start')
                else:
                    status.append('present')
                for n in matching.loc[matching['Parent']==i].iloc[1:n_cds[-1],4]:
                    if n == False:
                        status.append('missing_middle')
                    else:
                        status.append('present')
                if matching.loc[matching['Parent']==i].iloc[-1,4] == False:
                    status.append('missing_end') 
                else:
                    status.append('present')
            # add status to matching dataframe
            matching.loc[matching['Parent']==i,'status'] = status
    matching['status'] = matching['status'].fillna('present')
    matching = matching.rename(columns={3:'start',4:'end'})
    matching.to_csv(name+'-CDS-status.csv', index=None, header=True)
    
    summary = pd.DataFrame({'status':['present','missing_all_cds','missing_start','missing_middle','missing_end','total_cds']}).set_index('status')
    summary = summary.join(matching.groupby('status',as_index=False)
                 .count()[['status','Exist']]
                 .set_index('status'), on='status',how='outer').fillna(0).astype('int')
    summary.loc['total_cds'] = summary.Exist.sum()
    print("\033[97;46m {}\033[0;0m".format("CDS status summary: "))
    print()
    print(tabulate(summary,headers=['Status', 'Count']))
    print()
    return matching


@timer_func
def main():
    parser = argparse.ArgumentParser(
        description="examine the similarity of an old and new gff for one genome and produce a list of missing CDS",
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
    except FileNotFoundError:
        print("New gff is not in the current directory")

    # RUN FUNCTIONS
    old_gff_df = read_in_gff(old_gff)
    new_gff_df = read_in_gff(new_gff)
    find_missing_cds(old_gff_df, new_gff_df, name)
        
    return

main()

