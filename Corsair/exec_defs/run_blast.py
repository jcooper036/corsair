#!/usr/bin/env python3

import Corsair as cor
import os
import subprocess

def shell(command):
    return subprocess.check_output(command, shell=True, executable='/bin/bash')

def run_blast(ctl, iso_name):
    """
    Input: ctl file, isoform name
    Output: blast scaffold data in the isform file
    """
    
    ## load the iso object
    iso = cor.load_isoform(ctl, iso_name)

    ## check that the isoform has the appropriate data
    if not type(iso.ref_nt) == str:
        print("ERROR: Isoform object was not prepared for BLAST - no reference protein sequence")
        return None
    
    ## check if the BLAST database has been made, make it if not
    for genome in ctl.genome_paths:
        if not os.path.exists(genome + '.nhr'):
            print("UPDATE: making BLAST database for " + genome)
            command = 'makeblastdb -in ' + genome + ' -parse_seqids -dbtype nucl'
            shell(command)
        if not os.path.exists(genome + '.fai'):
            print("UPDATE: making faidx database for " + genome)
            command = 'samtools faidx ' + genome
            shell(command)

    ## run blast - output is a scaffold file in the gene folder. iso object knows them too
    for species in ctl.species:
        for genome in ctl.genome_paths:
            if species in genome:
                command = 'tblastn -outfmt "6 sseqid" -query <(echo ' + iso.ref_aa + ') -db ' + genome + ' -max_target_seqs 1| head -n 1'
                var = str(shell(command))
                var = var.lstrip("b'").rstrip("\\n'")
                iso.add_scaffold(species, var)

    ## save the iso object
    cor.save_isoform(ctl, iso)