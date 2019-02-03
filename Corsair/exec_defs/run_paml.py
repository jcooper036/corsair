#!/usr/bin/env python3

import Corsair as cor
import os

def run_paml_M7M8(ctl, iso_name):
    """
    Input: control object, isoform name, alinger name
    Output: PAML results file, updated isoform object
    """
    ## load the isoform object
    iso = cor.load_isoform(ctl, iso_name)

    print('Running PAML(M7 vs M8) for {}'.format(iso.name))

    ## write the control file
    cor.m7_m8_control_file(ctl, iso)

    ## run codeml
    command = ctl.mod_path + 'Corsair/bin/paml/4.9e/bin/codeml ' + iso.paml_control_file(ctl)
    cor.shell(command)

    ## remove all the files codeml leaves behind
    try:
        for f in ['2NG.dN','2NG.dS','2NG.t','4fold.nuc','lnf','rst','rst1','rub']:
            os.remove(iso.iso_files(ctl) + f)
            os.remove(os.getcwd() + f)
    except:
        pass
    
    print("PAML complete for {}".format(iso.name))  

    ## save the isoform object
    cor.save_isoform(ctl, iso)

def run_paml_M8M8a(ctl, iso_name):
    """
    Input: control object, isoform name, alinger name
    Output: PAML results file, updated isoform object
    """
    ## load the isoform object
    iso = cor.load_isoform(ctl, iso_name)

    print(iso.paml_results)

    ## check the previous p-value:
    if not iso.paml_results['pval']['M7M8'] < 0.05:
        return None

    print('Running PAML(M8 vs M8a) for {}'.format(iso.name))

    ## write the control file
    cor.m8_m8a_control_file(ctl, iso)

    ## run codeml
    command = ctl.mod_path + 'Corsair/bin/paml/4.9e/bin/codeml ' + iso.paml_control_file(ctl)
    cor.shell(command)

    ## remove all the files codeml leaves behind
    try:
        for f in ['2NG.dN','2NG.dS','2NG.t','4fold.nuc','lnf','rst','rst1','rub']:
            os.remove(iso.iso_files(ctl) + f)
    except:
        pass
    
    print("PAML complete for {}".format(iso.name))  

    ## save the isoform object
    cor.save_isoform(ctl, iso)