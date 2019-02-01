#!/usr/bin/env python3

import Corsair as cor

def results_processing(ctl, iso_name):
    """
    Input: control object, iso name
    Output: none, executes the different results processing which save the files
    """    
    
    ## load the isoform object
    iso = cor.load_isoform(ctl, iso_name)

    ## p-values should already be calculated

    ## BEB site alignment
    cor.beb_site_processing(ctl, iso_name)

    ## dNdS averages and maxes
    cor.dnds_processing(ctl, iso_name)

    ## save the isoform object
    cor.save_isoform(ctl, iso)