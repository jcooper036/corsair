#!/usr/bin/env python3

import Corsair as cor

def align_and_paml(ctl, iso, aligner):
    """
    Input: control object, isoform name
    Output: results in the isoform folder in the project directory
    
    Temporary def, for practicing doing everything. Need to keep in mind that variable 
    passing needs to be such that everything can be run by itself, and not kill memory,
    and be divided onto as many CPUS as possible, and it needs to be able to quit and
    not lose everything.
    """
    
    ## align
    cor.run_alignment(ctl, iso, aligner)
    
    ## back translate
    cor.back_translate(ctl, iso, aligner)

    ## build the tree
    cor.build_tree(ctl, iso)

    ## run PAML
    cor.run_paml_M7M8(ctl, iso, aligner)

    ## load PAML results

    ## check p-value

    ## run aligners and PAML again if necessary

    ## run M8-M8a if necessary