# Transfer annotations
These scripts will run liftoff, fix the problem with missing phase information and then evaluate the differences in protein and nucleotide sequence between the original gff and the liftoff transferred gff.

**The reccomended pipeline is:**
1. `run-and-optimise-liftoff.py`
2. Retrieve the liftoff with lowest flank number and the fewest missing features
3. `liftoff2apollo.py`
4. `gff-phase-finder.py`
5. `gff_protein_change_finder.py`

## The scripts

### `run-liftoff.sh`

Run one instance of liftoff with default parameters and `-polish`. Add it to a loop and run it as many times as you want.

**Usage:**
```
bash run-liftoff.sh target_genome.fasta reference_genome.fasta reference_gff.gff
```
**Output:**

- Standard [Liftoff output](https://github.com/agshumate/Liftoff#output)

#
### `run-and-optimise-liftoff.py`

Automate running liftoff with different `-flank` and `-distance` parameters to get the fewest unmapped features.

**Usage:**
```
run-and-optimise-liftoff.py [-h] -gff OLD_GFF -old OLD_GENOME_FASTA -new NEW_GENOME_FASTA [-p [PLOT [PLOT ...]]]

optional arguments:
  -h, --help            show this help message and exit
  -gff OLD_GFF, --old_gff OLD_GFF
                        The original gff that needs to be transfered
  -old OLD_GENOME_FASTA, --old_genome_fasta OLD_GENOME_FASTA
                        The old genome that corresponds to the old gff
  -new NEW_GENOME_FASTA, --new_genome_fasta NEW_GENOME_FASTA
                        The new genome that you want to transfer annotations to
  -p [PLOT [PLOT ...]], --plot [PLOT [PLOT ...]]
                        If true, plots will be made for the number of unmapped features and the number of lines in the resulting gff files across each change  
                        in flank.
```
**Output:**

- Standard [Liftoff output](https://github.com/agshumate/Liftoff#output)

- `[Old_Genome_name]-to-[New_Genome_name]-UNMAPPED-FEATURE-COUNTS.csv` - summary of the number of unmapped features for each level of flank tested

- `[Old_Genome_name]-to-[New_Genome_name]-UNMAPPED-FEATURE-COUNTS.png"` - a plot of the number of unmapped features

#
### `liftoff2apollo.py`

Convert liftoff GFF3 outputs into Apollo friendly GFF3 files.
Liftoff adds several fields to the gff attributes column that are unused by Apollo.

This script will remove added attributes and convert the header to Apollo format which uses this:
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

The removed attributes are as follows:
['coverage', 'sequence_ID', 'valid_ORFs', 'extra_copy_number', 'copy_num_ID']

**Usage:**
```
liftoff2apollo.py [-h] -l LIFTOFF_GFF -f GENOME_FASTA

optional arguments:
  -h, --help            show this help message and exit
  -l LIFTOFF_GFF, --liftoff_gff LIFTOFF_GFF
                        A gff file produced by liftoff
  -f GENOME_FASTA, --genome_fasta GENOME_FASTA
                        A genome fasta file that corresponds to the liftoff gff
```
**Output:**
- `liftoff-gff-name_APOLLO.gff` - A gff format annotation file


#
### `gff-phase-finder.py`
Calculates GFF phase information for CDS which lack it. This was made for liftoff gff outputs which lack this information.
NOTE! GFF phase is not the same as reading frame! Please see the [Sequence Ontology documentation](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

**Usage:**
```
gff-phase-finder.py [-h] -gff LIFTOFF_GFF

optional arguments:
  -h, --help            show this help message and exit
  -gff LIFTOFF_GFF, --liftoff_gff LIFTOFF_GFF
                        The gff without phase information (made for liftoff)
```
**Output:**

- `your-gff-name-PHASE-CORRECTED.gff` - A gff format annotation file

#
### `gff_indexer.py`

Forcefully adds a numbered index to the attributes of each row in a gff to act as a unique identifier just in case there isn't one

**Usage:**
```
gff_indexer.py [-h] -g GFF

optional arguments:
  -h, --help         show this help message and exit
  -g GFF, --gff GFF  The original gff for the genome that needs annotating
```
**Output:**

- `your-gff-name_indexed.gff` - A gff format annotation file

#
### `gff_protein_change_finder.py`
Examine the similarity of an old and new gff for one genome. 

The output will report changes in the nucleotide and protein sequence

**Usage:** 

```
gff_protein_change_finder.py [-h] -gff1 OLD_GFF -gff2 NEW_GFF -f1 OLD_GENOME -f2 NEW_GENOME -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -gff1 OLD_GFF, --old_gff OLD_GFF
                        The original gff for the genome
  -gff2 NEW_GFF, --new_gff NEW_GFF
                        The new gff for the genome to compare against the original one
  -f1 OLD_GENOME, --old_genome OLD_GENOME
                        A genome fasta file that corresponds to the liftoff gff
  -f2 NEW_GENOME, --new_genome NEW_GENOME
                        A genome fasta file that corresponds to the liftoff gff
  -o OUTPUT, --output OUTPUT
                        The name prefix to give your outputs
```
**Output:**

- `your-output-name-missing_cds_lists.npy` - a python dictionary of gene models with the following changes:
  - _'nucl_change'_, additions/deletions to their nucleotide length 
  - _'missing_start'_, missing CDS at the start
  - _'missing_middle'_, missing CDS at the middle
  - _'missing_end'_, missing CDS at the end of the gene model

- `your-output-name-STATS.csv` - a comma delimited table with the following information:
  - _ID_, The unique CDS ID
  - _strand_match_, whether both compared sequences are on the same strand
  - _nucl_old_length_, nucleotide sequence length for the original gff
  - _nucl_new_length_, nucleotide sequence length for the new gff
  - _nucl_length_diff_, difference in nucleotide sequence length
  - _nucl_levenshtein_ratio_, nucleotide sequence similarity calculated with [fuzzywuzzy](https://pypi.org/project/fuzzywuzzy/) sequence matching with levenshtein ratio
  - _nucl_score_, nucleotide sequence alignment score using BLOSUM62
  - _aa_old_length_, amino acid sequence length for the original gff
  - _aa_new_length_, amino acid sequence length for the new gff
  - _aa_length_diff_, amino acid in nucleotide sequence length
  - _aa_levenshtein_ratio_, amino acid sequence similarity calculated with [fuzzywuzzy](https://pypi.org/project/fuzzywuzzy/) sequence matching with levenshtein ratio
  - _aa_score_, amino acid sequence alignment score using BLOSUM62

#
### `find-reciprocal-best-hits.py` 
Calculate reciprocal best hits and make some nice histograms and density plots
Based on code from [here](https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/05-blast_for_rbh.html)

**Usage:**
```
find-reciprocal-best-hits.py [-h] -s1 QUERY_AA -s2 SUBJECT_AA -o OUTPUT

optional arguments:
  -h, --help    show this help message and exit
  -s1 QUERY_AA, --query_aa QUERY_AA   
                one set of aa sequence
  -s2 SUBJECT_AA, --subject_aa SUBJECT_AA    
                second set of aa sequence
  -o OUTPUT, --output OUTPUT    
                Name for outputs
```
**Output:**
Blast results for forward and reverse reads:
`outputName-fwd-results.tab`
`outputName-rev-results.tab`

Reciprocal best hits:
`Cposadasi-test-rbbh.csv`

histograms of normalised bitscores for forward, reverse and reciprocal reads:
![Forward and reverse normalised bitscores histogram](images-for-git/Cposadasi-test_norm_bitscores_histograms.png)
![RBH normalised bitscores histogram](images-for-git/Cposadasi-test_RBH_norm_bitscores_histogram.png)

2D density plots of normalised bitscores for forward, reverse and reciprocal reads
![Forward and reverse normalised bitscores 2D density](images-for-git/Cposadasi-test_norm_bitscores_2D_density_histogram.png)
![RBH normalised bitscores 2D density](images-for-git/Cposadasi-test_RBH_norm_bitscores_2D_density_histogram.png)