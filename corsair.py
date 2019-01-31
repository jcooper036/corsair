#!/usr/bin/env python3
import Corsair as cor

restart = False
sample_ctl = '/Users/Jacob/corsair/primates/primates.ctl'

## parse the ctl file, initialize the control object
ctl = cor.load_ctl(sample_ctl)

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

## standard data output for all genes
for iso in ctl.gene_list['all']:
    cor.corsair_results(ctl, iso)

    ## load the isoform object
    iso = cor.load_isoform(ctl, iso)
    print(iso.paml_results)


