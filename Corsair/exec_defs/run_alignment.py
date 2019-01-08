#!/usr/bin/env python3

import Corsair as cor
import os

def run_alignment(ctl, iso_name, aligner):
    """
    Input: control object, isoform name, aligner name
    Output: Alignment file in isoform folder
    """

    ## load the isoform object
    iso = cor.load_isoform(ctl, iso_name)

    ## specifiy the input and output file names
    infile = ctl.project_path + 'temp/' + iso.name + '_' + aligner + '.fasta'
    outfile = iso.alignment_file(ctl, aligner)

    ## make the file
    cor.write_fasta(infile, iso.blast_prot)

    ## different commands based on aligner
    if aligner == 'clustal':
        command = 'clustalo -i ' + infile +' -o ' + outfile + ' --force'
    
    if aligner == 'tcoffee':
        command = 't_coffee ' + infile + ' -outfile ' + outfile + ' -multi_core -quiet -output=fasta'
    
    if aligner == 'muscle':
        command = 'muscle -quiet -in ' + infile + ' -out ' + outfile
    
    cor.shell(command)
    print('Sequences for ' + iso.name + ' aligned with ' + aligner + '.')

    ## cleans up a file that is sometimes made, but useless
    try:
        os.remove(iso + '_coff.dnd')
    except:
        pass

    ## load the alignment as a dictionary into the isoform object
    iso.alignment[aligner] = cor.read_fasta(outfile)
    
    ## save the isoform object
    cor.save_isoform(ctl, iso)
