# Corsair
A python module for running Phylogenetic Analysis by Maximum Likeliehood (PAML, Yang 2005) on a genomic scale.
v1.00

## Installation

Corsair comes with several of the decencies it needs to run. Namely:
- tBlastn
- Exonerate
- Mafft
- Tcoffee
- Clustal Omega


To install the package, run pip:
```bash
pip install foobar
```

## Pre-run requirements
- Control file
    - fill out a control file. A sample control file can be generated in the working directory with:

```python
import Corsair as cor
cor.sample_ctl()
```

- Parse reference CDS Fasta file
    - the reference CDS file needs to contain ONLY the gene or transcript name in the header
    - for some tools to do this, see: 
- Genome names need to match the clade names in the tree, followed by _
    - the genome name has to fit the exact (case sensitve) name give in the tree

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