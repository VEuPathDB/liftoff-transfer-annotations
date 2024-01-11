
#
### `gff-phase-finder.py`
Calculates GFF phase information for CDS which lack it. This was made for liftoff gff outputs which lack this information.
NOTE! GFF phase is not the same as reading frame! Please see the [Sequence Ontology documentation](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

*This script is relatively slow and does not deal with instances where the start phase is not 0. I have since found a better piece of code to do this: [AGAT](https://agat.readthedocs.io/en/latest/tools/agat_sp_fix_cds_phases.html)*

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