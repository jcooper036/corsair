#!/usr/bin/env python3
import Corsair as cor
import sys


restart = False

# sample_ctl = '/Users/Jacob/corsair/primates/primates_full.ctl'
sample_ctl = '/Users/Jacob/corsair/murinae/mouse.ctl'

## parse the ctl file, initialize the control object
ctl = cor.load_ctl(sample_ctl)

## loads the gene list from the control object
ctl.load_gene_list()

## for the first time only - will over-write saves otherwise
if restart:
    cor.corsair_initialize(ctl)

## for each gene in that list
for iso in ctl.gene_list:

    # # just do blast and exonerate
    if restart:
        cor.blast_and_exonerate(ctl, iso)

    # # just do the alignment and PAML (if Blast and Exonerate are already done)
    if restart:
        cor.align_and_paml(ctl, iso)

    # ## run all the results gathering functions for the whole gene list
    cor.read_paml_output(ctl, iso)

    ## read get results
    cor.results_processing(ctl, iso)

cor.write_output(ctl)
