#!/usr/bin/env python3

import Corsair as cor

def corsair_execs(ctl, iso):
    """
    Input: control object, isoform name
    Output: results in the isoform folder in the project directory
    
    Temporary def, for practicing doing everything. Need to keep in mind that variable 
    passing needs to be such that everything can be run by itself, and not kill memory,
    and be divided onto as many CPUS as possible, and it needs to be able to quit and
    not lose everything.
    """

    ## blast
    cor.run_blast(ctl, iso)

    ## parse the blast results into scaffolds

    ## exonerate on the scaffolds

    ## align

    ## build the tree

    ## run PAML

    pass