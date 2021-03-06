# Corsair
A python module for running Phylogenetic Analysis by Maximum Likeliehood (PAML, Yang 2005) on a genomic scale.
v1.00
jcooper036@gmail.com

## Installation

MAC users:
Corsair comes with several of the decencies it needs to run. So that the program can be run in a distributed manner, the /bin/ of this module actually has all the depenecies.
- Clustal Omega
- Exonerate
- Tcoffee
- Muscle
- paml

Linux users:
Corsair does not have the pre-compiled execs for linux. Run this command to install them:
apt-get install -y python3 python3-pip emboss exonerate clustalo muscle t-coffee paml ncbi-blast+ nfs-common curl


Some programs need to be present locally. They should be available to install with brew or apt-get
- blast (specifically tblastn)
- python3

There are also some python packages that are requied, but are easy to find with pip. Hopefully we didn't miss any - if any "cannot find module" error pop up, try installing that python module.
- subprocess
- scipy
- numpy
- pickle
- pyfaidx

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

There are several python scripts in the main directory that contain examples of how to string together the various commands.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)