#!/usr/bin/env python3
import Corsair as cor


ctl = cor.load_ctl('/Users/Jacob/corsair/sample.ctl')
print(ctl)





def full_send():
    """
    Temporary def, for practicing doing everything. Need to keep in mind that variable 
    passing needs to be such that everything can be run by itself, and not kill memory,
    and be divided onto as many CPUS as possible, and it needs to be able to quit and
    not lose everything.
    """
        
    ## parse the ctl file

    ## load all the reference CDS sequences
    
    ## blast

    ## parse the blast results into scaffolds

    ## exonerate on the scaffolds

    ## align

    ## build the tree

    ## run PAML

    ## manage results
    pass