#!/usr/bin/env python3
import Corsair as cor

restart = False
sample_ctl = '/Users/Jacob/corsair/primates/primates.ctl'

## parse the ctl file
ctl = cor.load_ctl(sample_ctl)

## for the first time only - will over-write saves otherwise
if restart:
    cor.corsair_initialize(ctl)

## build the lists of different genes to run for each category
ctl = cor.build_gene_lists(ctl)

## while there is anything left to do
for i in range(len(ctl.aligners) + 1):
    
    ## for each aligner
    for aligner in ctl.aligners:
        
        ## for each gene in that list
        for iso in ctl.gene_list[aligner]:

            ## run all the execs on a gene by gene basis - allows for parallelizing this part
            # cor.align_and_m7m8(ctl, iso, aligner)

            # # just do blast and exonerate
            if restart:
                cor.blast_and_exonerate(ctl, iso, aligner)

            # # just do the alignment and PAML (if Blast and Exonerate are already done)
            if restart:
                cor.align_and_paml(ctl, iso, aligner)

            # ## run all the results gathering functions for the whole gene list
            cor.read_paml_output(ctl, iso, aligner)
    
    ## remake the aligners 'to do' lists
    ctl = cor.build_gene_lists(ctl)

## standard data output for all genes
for iso in ctl.gene_list['all']:
    cor.corsair_results(ctl, iso)

    ## load the isoform object
    iso = cor.load_isoform(ctl, iso)
    print(iso.paml_results)


