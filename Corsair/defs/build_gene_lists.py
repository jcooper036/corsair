#!/usr/bin/env python3

import Corsair as cor

def build_gene_lists(ctl, aligners):
    """
    Input: control object
    Output: returns control object
    """

    ## clear the old lists of aligners. notice how this does NOT include ctl.gene_list['all']
    for aligner in aligners:
        ctl.gene_list[aligner] = []

    ## for every gene
    for iso_name in ctl.gene_list['all']:
        
        ## load the isoform object
        iso = cor.load_isoform(ctl, iso_name)

        check = {}

        ## for each aligner, check to see if it has a p-value that is less than 0.05
        for idx, aligner in enumerate(aligners):
            check[aligner] = {'run' : False, 'pass' : False}
            if iso.paml_results[aligner]:
                if iso.paml_results[aligner]['pval'] < 0.05:
                    check[aligner]['pass'] = True
                    check[aligner]['run'] = True
                else:
                    check[aligner]['run'] = True

        ## if it does, and the next aligner does not, add it to the list of the next aligner
        for idx, aligner in enumerate(aligners):
            if idx == 0:
                if not check[aligner]['run']:
                    ctl.gene_list[aligner].append(iso.name)
            if (idx + 1) < len(aligners):
                if check[aligner]['pass'] and not check[aligners[(idx + 1)]]['run']:
                    ctl.gene_list[aligners[(idx + 1)]].append(iso.name)

        ## save the isoform object
        cor.save_isoform(ctl, iso)

    return ctl