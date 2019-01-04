#!/usr/bin/env python3

import Corsair as cor
import os
import pickle

def load_isoform(ctl):
    """
    input : control object (ctl)
    output : dictionary of isoform objects
    """
    isoforms = {}

    ## load the list of isoforms
    iso_list = []
    with open(ctl.gene_list, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if len(line) > 0:
                iso_list.append(line)
    
    ## load the isoforms
    missing = []
    for iso in iso_list:
        pkl_file = "{}/genes/{}/{}.pkl".format(ctl.project_path, iso, iso)
        if os.path.exists(pkl_file):
            with open(pkl_file, 'rb') as f:
                isoforms[iso] = pickle.load(f)
        else:
            missing.append(iso)
    
    ## if some weren't found
    if len(missing) > 0:
        print("WARNING: Some or all of the genes that attempted to load have not been initialized")
        print("Make sure to initialize and save them SEPERATELY or it will erase the data for all genes")
        for gene in missing:
            print('\t' + str(gene))


    return isoforms