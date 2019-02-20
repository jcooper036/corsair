#!/usr/bin/env python3
"""
Usage: by_gene.py <ctl_file>

Arguments:
    <isoform>     Name of the isoform to be run
    <ctl_file>    Path to ctl file

"""
#imports
import docopt
import Corsair as cor
import sys

## Initialize docopt
if __name__ == '__main__':

    try:
        arguments = docopt.docopt(__doc__)
        ctl_file = str(arguments['<ctl_file>'])
    except docopt.DocoptExit as e:
        print(e)

## parse the ctl file, initialize the control object
ctl = cor.load_ctl(ctl_file)

## loads the gene list from the control object
ctl.load_gene_list()

for iso_name in ctl.gene_list:

    # just do exonerate
    # cor.scaffolds_from_file(ctl, iso_name)
    cor.run_exonerate(ctl, iso_name)
    cor.load_species_sequences(ctl, iso_name)

    # # just do the alignment and PAML (if Blast and Exonerate are already done)
    cor.align_and_paml(ctl, iso_name)

    # ## run all the results gathering functions for the whole gene list
    cor.read_paml_output(ctl, iso_name)

    ## read get results
    cor.results_processing(ctl, iso_name)

## do the results stuff
cor.write_output(ctl)