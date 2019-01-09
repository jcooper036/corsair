#!/usr/bin/env python3

import Corsair as cor
import os

def read_paml_output(ctl, iso_name, aligner):
    """
    Input: control object, iso name, aligner name
    Output: PAML results file for that aligner, updated isoform object
    """    
    
    ## load the isoform object
    iso = cor.load_isoform(ctl, iso_name)

    ## this dictionary contains a file path as the value, otherwise, contains false
    for aligner in iso.paml_output_files:
        if iso.paml_output_files[aligner]:
            iso.paml_results[aligner] = read_paml_file(iso.paml_output_files[aligner])

    ## save the isoform object
    cor.save_isoform(ctl, iso)

def read_paml_file(file):
    """
    Input : PAML output file
    Output : Returns a dicitonary of all variables and values
    """
    results = {}
    with open(file, 'r') as f:
        for line in f.readlines():
            line = line.strip
            


    return results