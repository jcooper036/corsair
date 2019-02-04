#!/usr/bin/env python3

import Corsair as cor
from datetime import datetime

def write_gene_results(iso, ctl):
    """
    Input: Isoform object, control object
    Writes file: Isoform results file
    """
    
    results_file = iso.results_file(ctl)

    f = open(results_file, 'w')

    f.write('Gene name: ' + str(iso.name) + '\n')
    
    ## species run
    if iso.good_species:
        f.write('Number of species analyzed: {} {}\n'.format(len(iso.good_species), iso.good_species))
    else:
        f.write('ERROR - no comparable species were identified.')

    ## M7-M8 p-value
    if iso.paml_results['pval']['M7M8']:
        f.write('M7-M8 p-value: {}\n'.format(iso.paml_results['pval']['M7M8']))
    else:
        f.write('ERROR: No M7-M8 PAML p-value')

    ## M8-M8a pvalue 
    if iso.paml_results['pval']['M8M8a']:
        f.write('M8-M8a p-value: {}\n'.format(iso.paml_results['pval']['M8M8a']))

    ## BEB sites
    if iso.paml_results['beb_hit_sites']:
        f.write('BEB sites with posterior prob. > {}: {}\n'.format(ctl.beb_threshold, len(iso.paml_results)))
        for site in iso.paml_results['beb_hit_sites']:
            f.write('{}{} '.format(iso.paml_results['beb_hit_sites'][site]['ID'], site))
        f.write('\n\n')
    
    ## write the alingment
    if iso.alignment['clustal']:
        
        ## for the scale
        species_length = len(ctl.ref_species)
        sequence_length = len(iso.alignment['clustal'][ctl.ref_species])
        scale = scale_string(species_length, sequence_length)
        f.write(scale + '\n')
        
        ## write the sequence
        f.write(ctl.ref_species + '\t' + iso.alignment['clustal'][ctl.ref_species] + '\n')
        for species in iso.alignment['clustal']:
            if species != ctl.ref_species:
                f.write(species + '\t' + iso.alignment['clustal'][species] + '\n')

    print(iso.alignment['clustal'])

    f.close()


def scale_string(species_length, sequence_length):
    """Takes the length of a sequence, returns the string to print as the scale"""
    ret_string = ''
    for letter in range(species_length):
        ret_string += ' '
    ret_string += '\t'
    tens = int(sequence_length / 10)
    num = 10
    for idx in range(tens):
        adder = '        {}'.format(num)
        while len(adder) > 10:
            adder = adder[1:]
        ret_string += adder
        num += 10


    return ret_string