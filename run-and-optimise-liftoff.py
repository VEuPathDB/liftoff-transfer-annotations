import os
import sys
from pathlib import Path
from glob import glob
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import numpy as np
from tabulate import tabulate
import re


def run_liftoff(old_gff, old_genome, new_genome, new_name, old_name, distance, flank):
    '''
    runs the shell script that makes liftoff go and substitutes in the input files
    '''
    os.system("cut -f3 {} | sort | uniq | grep '#' -v > {}_features"\
              .format(str(old_gff), str(old_name)))
    os.system("liftoff -g {} -o {}_d{}_f{}_lifted_annotations.gff3 -f {}_features -u {}_d{}_f{}_unmapped-features -d {} -flank {} -polish {} {}"\
              .format(old_gff, new_name,str(distance),str(flank), old_name, old_name,str(distance),str(flank), str(distance), str(flank), new_genome, old_genome))


def optimise_distance(old_gff, old_genome, new_genome, new_name, old_name):
    '''
    Optimises the distance to use until a gff file is produced.
    '''
    distance=5
    while (any(new_name in s for s in glob("*.gff3")) == False) and (distance < 50):
        print("\033[33m {}\033[0;0m".format("Trying distance:"+str(distance)+"..."))
        run_liftoff(old_gff, old_genome, new_genome, new_name, old_name, distance, 0.0)
        distance+=10
    
    if distance > 50:
        print("Aborting.")
        raise Exception("There were more than 5 attempts to optimise distance and nothing happened.")

    return distance
 

def optimise_flank(old_gff, old_genome, new_genome, new_name, old_name, distance):
    '''
    optimise flank
    '''
    num_unmapped = []
    flank_tested = []
    for flank in [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]:
        print("\033[33m {}\033[0;0m".format("Trying flank "+str(flank)))
        unmapped_file=str(old_name)+"_d"+str(distance)+"_f"+str(flank)+"_unmapped-features"
        gff_file=str(new_name)+"_d"+str(distance)+"_f"+str(flank)+"_lifted_annotations.gff3"

        if os.path.exists(gff_file) == False:
            run_liftoff(old_gff, old_genome, new_genome, new_name, old_name, distance, flank)
            print("\033[33m {}\033[0;0m".format("liftoff for distance "+str(distance)+" and flank "+str(flank)+" now done. Moving on..."))
        elif os.path.exists(gff_file) == True:
            print("\033[33m {}\033[0;0m".format("gff file for distance "+str(distance)+" and flank "+str(flank)+" already exists. Moving on..."))
        
        with open(unmapped_file, 'r') as fp:
            for count, line in enumerate(fp):
                pass
        flank_tested.append(flank)
        num_unmapped.append(count + 1)

    flank_dict = dict(zip(flank_tested,num_unmapped))
    smallest_flank = {key:val for key,val in flank_dict.items() if val == min(flank_dict.values())}

    return smallest_flank, flank_dict


def plot_stats(df, old_name, new_name, best_distance):
    '''
    Plots unmapped features with changing flank
    '''
    fig, ax = plt.subplots(figsize=(10,7))
    plt.plot(df['Flank'], df['Unmapped Features'], label=best_distance)
    plt.title(old_name+" to "+new_name+" Unmapped Features")
    plt.legend(title="Best Distance")
    plt.tight_layout()
    plt.savefig(old_name+"-to-"+new_name+"-UNMAPPED-FEATURE-COUNTS.png")

def clean_up(df):
    best = df.loc[(df['Unmapped Features'] == df['Unmapped Features'].min())]
    best = best.loc[best['Flank']==best['Flank'].min()]
    flanks = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    flanks.remove(best['Flank'].item())
    paths=[]
    for flank in flanks:
        paths.append([Path(p) for p in  glob("*_f"+str(flank)+"*")])
    paths = [item for sublist in paths for item in sublist]
    for file in paths:
        if os.path.isfile(file):
            os.remove(file)
        else:
            # If it fails, inform the user.
            print("Error: %s file not found" % file)


def optimise_liftoff(old_gff, old_genome, new_genome, new_name, old_name, plot):
    '''
    refines what parameters are needed
    '''
    print("\033[32m {}\033[0;0m".format("Begining lift off to transfer annotations from "+old_name+" to "+new_name))
    print("Files being used: \n")
    print(tabulate([['Old gff', old_gff], ['Old genome', old_genome],['New genome', new_genome]], headers=['Type', 'File name']))
    print('')

    if (os.path.exists(new_name+"_d1_f0.0_lifted_annotations.gff3")==False):

        run_liftoff(old_gff, old_genome, new_genome, new_name, old_name, 1, 0.0)
        print("\033[46m {}\033[0;0m".format("Liftoff has run with default settings "))

    else:
        print("\033[46m {}\033[0;0m".format("Liftoff has already run with default settings"))

    if (os.path.exists(new_name+"_d1_f0.0_lifted_annotations.gff3")==True):
        if (os.path.getsize(old_name+'_d1_f0.0_unmapped-features') == 0):
            print("\033[46m {}\033[0;0m".format("Default run was a success and there are no unmapped features, hooray! "))
            sys.exit()  # Terminate the program with status code 0
        elif (os.path.getsize(old_name+'_d1_f0.0_unmapped-features') != 0):
            print("\033[33m {}\033[0;0m".format("Default setting ran but not all features were mapped :("))
            print("\033[33m {}\033[0;0m".format("Default distance of 1 is fine."))
            best_distance=1
            print("\033[33m {}\033[0;0m".format("Optimising flank..."))
            best_flank, all_flank_dict = optimise_flank(old_gff, old_genome, new_genome, new_name, old_name, 1)
            print("\033[32m {}\033[0;0m".format("The number of unmapped features for the best flank "+str(list(best_flank.keys()))+" is: "+str(list(best_flank.values()))))
            # make a dataframe to save and pass to plot_stats()
            df = pd.DataFrame(all_flank_dict, index=[0]).T.reset_index().rename(columns={'index':'Flank',0:'Unmapped Features'})
            df['Distance'] = best_distance
            df['old_genome'] = old_name
            df['new_genome'] = new_name
            df.to_csv(old_name+"-to-"+new_name+"-UNMAPPED-FEATURE-COUNTS.csv", index=None, header=True)
            print(df)

    elif (os.path.exists(new_name+"_d1_f0.0_lifted_annotations.gff3")==False):
        if (os.path.getsize(old_name+'_unmapped-features') == 0):
            print("\033[33m {}\033[0;0m".format("Default setting ran but no features were mapped :("))
            print("\033[33m {}\033[0;0m".format("Optimising distance..."))
            best_distance = optimise_distance(old_gff, old_genome, new_genome, new_name, old_name)
            print("\033[33m {}\033[0;0m".format("Best distance is "+str(best_distance) ))
            print("\033[33m {}\033[0;0m".format("Optimising flank..."))
            best_flank, all_flank_dict = optimise_flank(old_gff, old_genome, new_genome, new_name, old_name, best_distance)
            print(best_flank, all_flank_dict)
            # make a dataframe to save and pass to plot_stats()
            df = pd.DataFrame(all_flank_dict, index=[0]).T.reset_index().rename(columns={'index':'Flank',0:'Unmapped Features'})
            df['Distance'] = best_distance
            df['old_genome'] = old_name
            df['new_genome'] = new_name
            df.to_csv(old_name+"-to-"+new_name+"-UNMAPPED-FEATURE-COUNTS.csv", index=None, header=True)

    if plot is not None: 
        print("Making a plot...")
        plot_stats(df, old_name, new_name, best_distance)
        print("Plotted: " + old_name+"-to-"+new_name+"-UNMAPPED-FEATURE-COUNTS.png")
    
    #clean_up(old_name, new_name, best_distance, best_flank)

    return df


def main():
    '''
    executes the arguments and functions to run and optimise liftoff
    '''
    parser = argparse.ArgumentParser(
     formatter_class=argparse.RawDescriptionHelpFormatter,
     description='''\
Transfer old annotations to a new genome with liftoff and
optimise your parameters.
-------------------------------------------------------------

     ''',
     epilog="written by Helen Rebecca Davison") 
    parser.add_argument('-gff', '--old_gff', \
                    help="The original gff that needs to be transfered",
                    required=True)
    parser.add_argument('-old','--old_genome_fasta', \
                    help="The old genome that corresponds to the old gff",
                    required=True)
    parser.add_argument('-new','--new_genome_fasta', \
                    help="The new genome that you want to transfer annotations to",
                    required=True)
    parser.add_argument('-p','--plot',  action='store', nargs='*',\
                    help="If true, plots will be made for the number of unmapped features\
                        and the number of lines in the resulting gff files across each change\
                        in flank.")
    args = parser.parse_args()


    old_gff = args.old_gff
    old_genome = args.old_genome_fasta
    new_genome = args.new_genome_fasta
    plot = args.plot

    new_name = Path(new_genome).stem
    old_name = Path(old_genome).stem


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
            if "#" not in f.readline:
                raise Exception("No header lines found, check that your file is definitely a gff")
        gff_id = [line for line in open(old_gff) if line[:1] != '#'][0].split()[0]
    except FileNotFoundError:
        print("New genome fasta is not in the current directory")

    try:
        gff_id = [line for line in open(old_gff) if line[:1] != '#'][0].split()[0]
        with open(old_genome, "r") as f:
            found=False
            if any(line for line in f if gff_id in line):
                found=True
            if found == False:
                raise Exception("Contig IDs do not match between the old gff and old genome")
    except FileNotFoundError:
        print("Your gff file is not in the current directory")

    optimise_liftoff(old_gff, old_genome, new_genome, new_name, old_name, plot)
    df = pd.read_csv(glob("*UNMAPPED-FEATURE-COUNTS.csv"))
    clean_up(df)

main()