#!/usr/bin/env python3
import Corsair as cor

sample_ctl = '/Users/Jacob/corsair/primates/primates.ctl'

## parse the ctl file
ctl = cor.load_ctl(sample_ctl)

## for the first time only - will over-write saves otherwise
cor.corsair_initialize(ctl)

## loop through the gene list
for iso in ctl.gene_list:

    ## run all the execs on a gene by gene basis - allows for parallelizing this part
    cor.corsair_execs(ctl, iso)

    ## run all the results gathering functions for the whole gene list
    cor.corsair_results(ctl, iso)


