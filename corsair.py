#!/usr/bin/env python3
import Corsair as cor
import sys

## can take the gene list from the call
gene_list = False
if sys.argv[0]:
    if "gene_list=['" in sys.argv[0]:
        gene_list = sys.argv[0].split("[")[1].split("]")[0].replace("'", "").split(",")

restart = False
sample_ctl = '/Users/Jacob/corsair/primates/primates.ctl'

## parse the ctl file, initialize the control object
ctl = cor.load_ctl(sample_ctl)

## for the first time only - will over-write saves otherwise
if restart:
    cor.corsair_initialize(ctl)

## set the gene list to the thing we read in, if it is there
if gene_list:
    ctl.gene_list = gene_list

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

    ## print them to the terminal
    iso_ob = cor.load_isoform(ctl, iso)
    # print(iso_ob.paml_results)
