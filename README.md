# Transfer annotations
These scripts will run liftoff, fix the problem with missing phase information and then evaluate the differences in protein and nucleotide sequence between the original gff and the liftoff transferred gff.

**The reccomended pipeline is:**
1. Run `run-and-optimise-liftoff.py`
2. Retrieve the liftoff with lowest flank number and the fewest missing features
3. AGAT CDS phase fixer, [here](https://agat.readthedocs.io/en/latest/tools/agat_sp_fix_cds_phases.html)
   ```
   agat_sp_fix_cds_phases.pl --gff infile.gff -f fasta [ -o outfile ]
   ``` 
4. Extract CDS protein code for the source and lifted gff with [AGAT](https://agat.readthedocs.io/en/latest/tools/agat_sp_extract_sequences.html)
   ```
   agat_sp_extract_sequences.pl -g infile.gff -f infile.fasta -t cds -p [-o outfile]
   ```
6. Run `gff_missing_cds_finder.py`
7. Run `gff_protein_change_finder.py`

## The scripts

You can run everything at once with `run-transfer-annotations.sh`

## The scripts

### `run-transfer-annotations`
For this to work:
1. change "~/scripts/" infront of each of the python script lines (30, 64, and 67) to the path they are stored on your machine
2. install AGAT in a new environment called agatenv (or change the name in the script if you already have agat)
3. be in an environment with liftoff installed

**Usage:**
```
bash run-transfer-annotations.sh SOURCE_GFF.gff SOURCE_GENOME.fasta NEW_GENOME.fasta prefix

Arguments:
   SOURCE_GFF.gff - gff from original genome to be transferred
   SOURCE_GENOME.fasta - genome fasta associated with the source gff
   NEW_GENOME.fasta - the genome fasta that wat you want to lift the gff to
   prefix - a prefix for protein identity check output files
```

#
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
### `gff_missing_cds_finder.py`
Examine the similarity of an old and new gff for one genome and retireve a list of missing CDS.

**Usage:**
```
gff_missing_cds_finder.py [-h] -sg SOURCE_GFF -lg LIFTED_GFF -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -sg SOURCE_GFF, --source_gff SOURCE_GFF
                        The source gff for the genome
  -lg LIFTED_GFF, --lifted_gff LIFTED_GFF
                        The lifted gff for the genome to compare against the original one
  -o OUTPUT, --output OUTPUT
                        The name prefix to give your outputs
```
**Output:**
- `your-name-missing-cds-list.csv` - a list of transcript IDS that have a missing start, middle or end CDS

#
### `gff_protein_change_finder.py`
A script to annotate a gff with transfer vailidty information and protein identity matches 

The output will report changes in the nucleotide and protein sequence

**Usage:** 

```
gff_protein_change_finder.py [-h] -lcds LIFTED_CDS_AA -scds SOURCE_CDS_AA -nf NEW_FASTA -sf SOURCE_FASTA -lg LIFTED_GFF -sg SOURCE_GFF -m MISSING_CDS_LISTS -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -lcds LIFTED_CDS_AA, --lifted_cds_aa LIFTED_CDS_AA
                        the new CDS proteins retrieved with AGAT
  -scds SOURCE_CDS_AA, --source_cds_aa SOURCE_CDS_AA
                        the source CDS proteins retrieved with AGAT
  -nf NEW_FASTA, --new_fasta NEW_FASTA
                        the new genome fasta associated with the lifted gff
  -sf SOURCE_FASTA, --source_fasta SOURCE_FASTA
                        the source genome fasta associated with the source gff and lifted gff
  -lg LIFTED_GFF, --lifted_gff LIFTED_GFF
                        the gff associated with the new genome produced by liftoff and phase fixed with AGAT to be annotated with information
  -sg SOURCE_GFF, --source_gff SOURCE_GFF
                        the gff associated with the source genome and lifted gff
  -m MISSING_CDS_LISTS, --missing_cds_lists MISSING_CDS_LISTS
                        csv of missing cds with the headings 'missing_start', 'missing_middle', 'missing_end'
  -o OUTPUT, --output OUTPUT
                        output names
```
**Output:**
- `your-output-name_changes_marked.gff` - gff with the added attributes for 'has_missing_cds', 'has_internal_stops','valid_transfer','proteins_match_source'
- `your-output-name_protein-change-stats.csv` - a csv containing a summary of changes in CDS between source and lifted genes
- `your-output-name_CDS-summary.png` - A bar chart of the number of valid transfers and whether or not their protein sequence matches
- `your-output-name_ncRNA-pseudogene-summary.png` - A bar chart of the number and types of non-coding genes coloured by the validity of the transfer
- `your-output-name_CDS-failure-summary.png` - A bar chart of the failed CDS transfers and what is wrong with them

#
### `genome-stats-calculator.sh`
This script provides stats for genomes: 
Number of contigs, assembly length, longest contig, N50, N90, GC%

It can also optionally count gff features, but you might need to change lines 67 onwards to match the contents of your gff (e.g. "ncRNA_gene" might need to be "ncRNA")
It does not check that the gff matches the genome, so make sure you use the right one.

This was made to deal with VEuPathDB gff3 syntax

**Usage:**
```
script genome.fasta [corresponding-gff.gff]
```

**Output:**
- genome-name-summary.txt
