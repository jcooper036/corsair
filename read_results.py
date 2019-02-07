#!/usr/bin/env python3
"""
Usage: by_gene.py <ctl_file>

Arguments:
    <ctl_file>    Path to ctl file

"""
#imports
import docopt
import Corsair as cor

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

## read get results
for iso_name in ctl.gene_list:

    ## align
    # cor.run_alignment(ctl, iso_name)
    
    ## trim to minimum alignment
    cor.trim_sequences(ctl, iso_name)

    ## back translate
    cor.back_translate(ctl, iso_name)



## do the results stuff
cor.write_output(ctl)