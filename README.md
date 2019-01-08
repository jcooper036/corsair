# Corsair
A python module for running Phylogenetic Analysis by Maximum Likeliehood (PAML, Yang 2005) on a genomic scale.
v1.00

## Installation

Corsair comes with several of the decencies it needs to run. So that the program can be run in a distributed manner, the /bin/ of this module actually has all the depenecies.
- tBlastn
- samtools faidx
- Mafft
- Clustal Omega
- Exonerate
- Tcoffee

If you would like to use a local version instead, then these are the files that you need to alter:
- Corsair/exec_defs/run_blast.py
Look for the subprocess commands, change the executable path to whatever you like.


To install the package, run pip:
```bash
pip install foobar
```

## Pre-run requirements -- IMPORTANT!
- Control file
    - fill out a control file. To generate a sample control file:

```python
import Corsair as cor
cor.sample_ctl(file_name)
```

- Parse reference CDS Fasta file and gene list
    - the reference CDS file needs to contain ONLY the gene or transcript name in the header
    - genes in the gene list must EXACTLY match a header in the fasta file
    - for some tools to do this, see:

- Genome names need to match the clade names in the tree, followed by '_', and need to end as .fasta
    - the genome name has to fit the exact (case sensitve) name give in the tree
    - ex: Hsapiens_GRCH36.fasta

## Usage
Corsair requires a control file to be executed. build the ctl object with the following:

```python
import Corsair as cor
ctl = cor.Ctlog(sample.ctl)
```
Then, any command corsair uses takes just the ctl object as the input. ex:
```python
import Corsair as cor
ctl = cor.load_ctl('sample.ctl')
cor.blast(ctl)
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)