import csv
import argparse

parser = argparse.ArgumentParser(
     description="Forcefully adds a numbered index to the attributes of each row\
          in a gff to act as a unique identifier just in case there isn't one",
     epilog="written by Helen Davison")
parser.add_argument('-g', '--gff', \
                    help="The original gff for the genome that needs annotating",
                    required=True)
args = parser.parse_args()

gff = args.gff

def gff_indexer(gff):
    print("\033[46m {}\033[0;0m".format("Indexing the gff " + gff + " "))
    # count the number of non-comment rows to use as a stop point for the indexer
    with open(gff, "r", encoding="utf8") as fin:
        counter = enumerate(fin)
        comment_lines = 0
        for count, value in counter:
            if any("#" in s for s in value):
                    comment_lines+=1
        num_lines = count-comment_lines
    print("\033[33m {}\033[0;0m".format("The gff has " + str(count+1) + " rows with " +str(comment_lines)+" comment lines."))
    fin.close()

    gff_name = gff.replace(".gff","")
    fin = open(gff, "r", encoding="utf8")
    fout = open(str(gff_name)+"_indexed.gff", "w", encoding="utf8", newline='')
    gff_reader = csv.reader(fin, delimiter="\t")
    gff_writer = csv.writer(fout, delimiter="\t")
    i = 1
    # read each line in and if its not a comment line, get all of the values
    for row in gff_reader:
        if any("#" in s for s in row):
            gff_writer.writerow([str(line) for line in row])            
        else:
            seqname, source, feature, start, end, score, strand, frame, attribute = row
            # make a dictionary for the attributes add the index number
            attr_dict = attribute.split(';')
            pairs = [tuple(x.split('=')) for x in attr_dict]      
            attr_dict = dict(pairs)
            attr_dict['index_num'] = i

            # convert the index back into a gff attribute format
            new_attribute=""
            dict_len = len(attr_dict)
            num_attr_added = 1
            for key, value in attr_dict.items():
                if num_attr_added < dict_len:
                    new_attribute += str(key) +"="+ str(value)+";"
                    num_attr_added+=1
                elif num_attr_added >= dict_len-1:
                    new_attribute += str(key) +"="+ str(value)
                    
                    # write new attribute to the indexed gff
            new_row = [str(seqname), str(source), str(feature), str(start), str(end), str(score), str(strand), str(frame), str(new_attribute)]                                    
            gff_writer.writerow(new_row)
            i+=1
    fin.close()
    fout.close()
    print("\033[33m {}\033[0;0m".format("Index numbers 1 to " + str(num_lines+1) + " have been generated for the gff attributes"))    
    print("\033[42m {}\033[0;0m".format("Finished indexing " + gff + " "))
    
    return

gff_indexer(gff)