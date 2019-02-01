#!/usr/bin/env python3

import Corsair as cor
from copy import deepcopy

def beb_site_processing(ctl, iso_name):
    """
    Input: control object, iso name
    Output: BEB site alignment, saved isoform object
    """    
    
    ## load the isoform object
    iso = cor.load_isoform(ctl, iso_name)
    
    ## make a positiion variable that we can incriment
    for site in iso.paml_results['beb_sites_raw']:
        iso.paml_results['beb_sites_raw'][site]['pos'] = site

    ## sort the mask - it is a set before this
    iso.mask_sort = sorted(list(iso.mask))

    ## iterate over the mask to tell us what sites need re-indexing
    for m_site in iso.mask_sort:
        mask_site = m_site + 1
        for site in iso.paml_results['beb_sites_raw']:
            if iso.paml_results['beb_sites_raw'][site]['pos'] >= mask_site:
                iso.paml_results['beb_sites_raw'][site]['pos'] += 1
    
    ## make a new data set
    iso.paml_results['beb_sites'] = {}
    for site in iso.paml_results['beb_sites_raw']:
        site_key = iso.paml_results['beb_sites_raw'][site]['pos']
        site_id = iso.paml_results['beb_sites_raw'][site]['ID']
        site_beb = iso.paml_results['beb_sites_raw'][site]['BEB']
        iso.paml_results['beb_sites'][site_key] = {"ID" : site_id, "BEB" : site_beb}

    ## save the isoform object
    cor.save_isoform(ctl, iso)