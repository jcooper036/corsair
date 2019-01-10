#!/usr/bin/env python3
import Corsair as cor

sample_ctl = '/Users/Jacob/corsair/primates/primates.ctl'

## parse the ctl file
ctl = cor.load_ctl(sample_ctl)

## for each gene in that list
for iso in ctl.gene_list['all']:
    
    # ## run all the results gathering functions for the whole gene list
    cor.read_paml_output(ctl, iso, 'clustal')