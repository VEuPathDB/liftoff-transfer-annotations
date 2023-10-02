import pandas as pd
import numpy as np

def progress_bar(progress, total):
    percent = (100*progress/float(total))
    bar = "█"*int(percent) + "-" * (100 - int(percent))
    print(f"\r|{bar}| {percent:.2f}%", end="\r")

# function to load the gff into a pandas dataframe and unpack the attributes column to access the unique IDs
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
    df = df.fillna(value=np.nan)
    print("\033[32m {}\033[0;0m".format("Done!\n"))
    return df

def find_closest_CDS_coords(old_gff_df,new_gff_df, old_type, new_type):
    def overlap(start1, end1, start2, end2):
        return max(max((end2-start1), 0) - max((end2-end1), 0) - max((start2-start1), 0), 0)
    def find_match(coord_dictionary):
        for name, data in coord_dictionary.items():
            coords, cds = data
            cds.reset_index(inplace=True)
            cds['targets'] = [np.array((row['3'],row['4'])) for index, row in cds.iterrows()]   
            print('finding overlaps for '+name)     
            cds['overlaps'] = pd.Series(dtype='object')
            for index, target in cds['targets'].items():
                overlaps=[]
                for i in coords:
                    o = overlap(i[0], i[1], target[0],target[1])
                    if (o != 0) and (len(i)>0):
                        overlaps.append(i)
                if len(overlaps) >0:
                    cds.loc[cds.index==index, 'overlaps'] = [overlaps]               
                progress_bar(index+1,len(cds))
            for index, row in cds.loc[cds['overlaps'].isna()==False].iterrows():    
                distances = np.linalg.norm(row['overlaps']-target, axis=1)
                min_index = np.argmin(distances)
                closest = row['overlaps'][min_index]
                new_contig = row['0']
                old_contig = name
                old_start =  closest[0]
                old_end = closest[1]
                new_start = row['3']
                new_end = row['4']
                distance = distances[min_index]
                seq_length_change = ((new_end-new_start)-(old_end-old_start))
                annotation_change = (list(old_gff_subset.loc[(old_gff_subset['3']==old_start)&(old_gff_subset['4']==old_end)]['2'])[0], row['2'])
                dumpfile_ID = row['ID']
                veupath_ID = list(old_gff_subset.loc[(old_gff_subset['3']==closest[0]) & (old_gff_subset['4']== closest[1])]['ID'])
                overlaps = row['overlaps']
                if (new_contig==old_contig) & (len(veupath_ID)==1):
                    corresponding_coords.append([old_contig, new_contig, annotation_change, old_start, old_end, new_start, new_end, seq_length_change, distance, veupath_ID[0], overlaps, dumpfile_ID])
                elif (new_contig==old_contig) & (len(veupath_ID)>1):
                    for i in old_gff_subset.loc[(old_gff_subset['3']==closest[0]) & (old_gff_subset['2']== old_type) & (old_gff_subset['4']== closest[1])]['ID']:
                        corresponding_coords.append([old_contig, new_contig, annotation_change, old_start, old_end, new_start, new_end, seq_length_change, distance, i, overlaps, dumpfile_ID])
            for index, row in cds.loc[cds['overlaps'].isna()==True].iterrows():
                corresponding_coords.append([np.nan, row['0'], np.nan, np.nan, np.nan, row['3'], row['4'], np.nan, np.nan, np.nan, overlaps, row['ID']])
        df = pd.DataFrame(corresponding_coords[1:], columns=corresponding_coords[0])
        print("\nClosest coordinates, done")
        return df
    corresponding_coords = [["old_contig","new_contig","types_compared","old_closest_start","old_closest_end","new_start", "new_end", "seq_length_change", "distance_to_closest_coord", "closest_old_ID", "all_overlapping_coords", "new_ID"]]
    new_cds = new_gff_df[new_gff_df['2']==new_type]
    contig_names = new_cds['0'].unique()
    old_gff_subset = old_gff_df[(old_gff_df['2']==old_type)& (old_gff_df['0'].isin(contig_names))]    
    # find closest CDS
    neg_coord_dictionary = {}
    pos_coord_dictionary = {}
    for index, name in enumerate(contig_names):
        print("\nLooking for contig "+ name+ "\n")
        coords_array_pos = np.array(list(zip(old_gff_subset.loc[(old_gff_subset['0']==name)&(old_gff_subset['6']=='+')]['3'],old_gff_subset.loc[(old_gff_subset['0']==name)&(old_gff_subset['6']=='+')]['4'])))
        new_cds_subset_pos = new_cds.loc[(new_cds['0']==name) & (new_cds['6']=='+')]
        coords_array_neg = np.array(list(zip(old_gff_subset.loc[(old_gff_subset['0']==name)&(old_gff_subset['6']=='-')]['3'],old_gff_subset.loc[(old_gff_subset['0']==name)&(old_gff_subset['6']=='-')]['4'])))
        new_cds_subset_neg = new_cds.loc[(new_cds['0']==name) & (new_cds['6']=='-')]
        if len(coords_array_pos)>0:
            pos_coord_dictionary[name] = [coords_array_pos, new_cds_subset_pos]
        if len(coords_array_neg)>0:
            neg_coord_dictionary[name] = [coords_array_neg, new_cds_subset_neg]
        progress_bar(index+1,len(contig_names))
    neg_df = find_match(neg_coord_dictionary)
    pos_df = find_match(pos_coord_dictionary)
    all_df = neg_df.append(pos_df)
    all_df['all_overlapping_coords'] = all_df['all_overlapping_coords'].astype('str')
    all_df.drop_duplicates()
    return all_df

#ncbi_gff = 'Coccidioides_silveira\FungiDB-65_CposadasiiSilveira2022.gff'
#liftoff_gff = 'Coccidioides_silveira\FungiDB-65_CposadasiiSilveira2022_Genome_d1_f0.1_lifted_annotations_polished_APOLLO.gff-PHASE-CORRECTED.gff'
#old_gff_df = read_in_gff(ncbi_gff)
#new_gff_df = read_in_gff(liftoff_gff)
#old_gff_df.to_csv("Coccidioides_silveira\FungiDB-65_CposadasiiSilveira2022_gff_db.csv", index=None, header=True)
#new_gff_df.to_csv("Coccidioides_silveira\FungiDB-65_CposadasiiSilveira2022_Genome_d1_f0.1_lifted_annotations_polished_APOLLO.gff-PHASE-CORRECTED_gff_db.csv", index=None, header=True)

old_gff_df = pd.read_csv("Coccidioides_silveira\FungiDB-65_CposadasiiSilveira2022_gff_db.csv")
new_gff_df = pd.read_csv("Coccidioides_silveira\FungiDB-65_CposadasiiSilveira2022_Genome_d1_f0.1_lifted_annotations_polished_APOLLO.gff-PHASE-CORRECTED_gff_db.csv")
df = find_closest_CDS_coords(old_gff_df, new_gff_df,'CDS', 'CDS')
df2 = find_closest_CDS_coords(old_gff_df, new_gff_df,'protein_coding_gene', 'protein_coding_gene')
#df.to_csv('Coccidioides_silveira\CposadasiiSilveira2022_old-to-new_cds_comparison.csv', index=None, header=True)
#df2.to_csv('Coccidioides_silveira\CposadasiiSilveira2022_old-to-new_gene_comparison.csv', index=None, header=True)
#len(new_cds.loc[new_cds['ID'].isin(list(all_df[all_df['closest_old_ID'].isna()]['new_ID']))]['gene_id'].unique())
df3 = find_closest_CDS_coords(new_gff_df,old_gff_df,'CDS', 'CDS')
df4 = find_closest_CDS_coords(new_gff_df,old_gff_df,'protein_coding_gene', 'protein_coding_gene')


new_cds = new_gff_df[new_gff_df['2']=='CDS']
old_cds = old_gff_df[old_gff_df['2']=='CDS']

new_gene = new_gff_df[new_gff_df['2']=='protein_coding_gene']
old_gene = old_gff_df[old_gff_df['2']=='protein_coding_gene']

len(new_cds.loc[new_cds['ID'].isin(list(df[df['closest_old_ID'].isna()]['new_ID']))]['gene_id'].unique())
len(old_cds.loc[old_cds['ID'].isin(list(df2[df2['closest_old_ID'].isna()]['new_ID']))]['gene_id'].unique())
df2.to_csv('Coccidioides_silveira\CposadasiiSilveira2022_old-to-new_cds_comparison.csv', index=None, header=True)

len(new_gene.loc[new_gene['ID'].isin(list(df3[df3['closest_old_ID'].isna()]['new_ID']))]['ID'].unique())
len(old_gene.loc[old_gene['ID'].isin(list(df4[df4['closest_old_ID'].isna()]['new_ID']))]['ID'].unique())
df3.to_csv('Coccidioides_silveira\CposadasiiSilveira2022_old-to-new_gene_comparison.csv', index=None, header=True)
df4.to_csv('Coccidioides_silveira\CposadasiiSilveira2022_new-to-old_gene_comparison.csv', index=None, header=True)


names = new_gff_df[2].unique()