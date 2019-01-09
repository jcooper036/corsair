#!/usr/bin/env python3

import Corsair as cor
import os

def run_paml_M7M8(ctl, iso_name, aligner):
    """
    Input: control object, isoform name, alinger name
    Output: PAML results file for that aligner, updated isoform object
    """
    ## load the isoform object
    iso = cor.load_isoform(ctl, iso_name)

    print('Running PAML(M7 vs M8) for {} aligned by {}'.format(iso.name, aligner))

    ## write the control file
    cor.m7_m8_control_file(ctl, iso, aligner)

    ## run codeml
    command = 'corsair/Corsair/bin/paml/4.9e/bin/codeml ' + iso.paml_control_file(ctl)
    cor.shell(command)

    ## remove all the files codeml leaves behind
    try:
        for f in ['2NG.dN','2NG.dS','2NG.t','4fold.nuc','lnf','rst','rst1','rub']:
            os.remove(iso.iso_files(ctl) + f)
    except:
        pass
    
    print("PAML complete for {}({})".format(iso.name, aligner))  

    ## save the isoform object
    cor.save_isoform(ctl, iso)

def run_paml_M8M8a(ctl, iso_name, aligner):
    """
    Input: control object, isoform name, alinger name
    Output: PAML results file for that aligner, updated isoform object
    """
    ## load the isoform object
    iso = cor.load_isoform(ctl, iso_name)

    print('Running PAML(M8 vs M8a) for {} aligned by {}'.format(iso.name, aligner))

    ## write the control file
    cor.m8_m8a_control_file(ctl, iso, aligner)

    ## run codeml
    command = 'corsair/Corsair/bin/paml/4.9e/bin/codeml ' + iso.paml_control_file(ctl)
    cor.shell(command)

    ## remove all the files codeml leaves behind
    try:
        for f in ['2NG.dN','2NG.dS','2NG.t','4fold.nuc','lnf','rst','rst1','rub']:
            os.remove(iso.iso_files(ctl) + f)
    except:
        pass
    
    print("PAML complete for {}({})".format(iso.name, aligner))  

    ## save the isoform object
    cor.save_isoform(ctl, iso)