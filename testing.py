#!/usr/bin/env python3

import os
import pickle
from scipy.stats import chi2

file='/Users/Jacob/corsair/primates/genes/EIF2AK2/EIF2AK2_files/EIF2AK2_PAML_M7M8_output.txt'

paml_results = {
    'logLvals' : {'M7' : False, 'M8' : False, 'M8a' : False},
    '2delta' : {'M7M8' : False, 'M8M8a' : False},
    'pval' : {'M7M8' : False, 'M8M8a' : False},
    'beb_sites_raw' : {}
}

def read_paml_file(file, paml_results):
    """
    Input : PAML output file and isoform object
    - Only reads the raw data from the file, no processing
    Output : Returns isoform object
    """
    
    tree_length = {}
    tree_length['M7'] = False
    tree_length['M8'] = False

    with open(file, 'r') as f:
        ready = False
        record_beb = False
        dnds_cap = False

        for line in f.readlines():
            line = line.strip()
            
            ## identify the model
            if "Model" in line:
                model = 'M8a' # this is default because it never shows up in M8a file
                if 'Model 7' in line:
                    model = 'M7'
                if 'Model 8' in line:
                    model = 'M8'
            
            ## find the log-liklihood values
            if "lnL" in line:
                line = line.replace(' ','').split(':')[3].split('+')[0]
                paml_results['logLvals'][model] = float(line)
            

            ## dNdS values
            if dnds_cap:
                if len(line) > 0:
                    # print(line.split('('))
                    for ele in line.split('('):
                        print(ele.split(')'))

                    # for ele in line.split('('):
                    #     if '(' in ele and len(ele) > 0:
                    #         dn = float(ele.split('(')[1])
                    #     if ')' in ele and len(ele) > 0:
                    #         ds = float(ele.split(')')[0])
                    #     if ds != 0:
                    #         dnds = dn/ds
                
            if 'Use runmode = -2 for ML pairwise comparison.)' in line:
                dnds_cap = True

            elif 'Model 7: beta (10 categories)' in line:
                dnds_cap = False

            ## Tree lengths
            if 'tree length =' in line:
                tree_length[model] = float(line.split()[3])

            ## BEB sites
            if record_beb:
                entry = line.split()
                if entry:
                    if entry[0].isdigit():
                        s_num = int(entry[0])
                        s_id = entry[1]
                        s_beb = float(entry[2].replace('*', ''))
                        paml_results['beb_sites_raw'][s_num] = {'ID' : s_id, 'BEB' : s_beb}

            if 'Bayes Empirical Bayes (BEB)' in line:
                ready = True
            if 'Pr(w>1)' in line and ready:
                record_beb = True
            if 'grid ' in line:
                record_beb = False

    print(tree_length)
    return paml_results


print(read_paml_file(file, paml_results))


def average_dn_ds(self):
    """Find the average dN/dS pairwise vales from the table in the PAML file"""
    results_exist = []
    self.dNdS = {}
    self.dN = {}
    self.dS = {}
    
    ## check what outputs there are
    for alg in ['clus', 'mus', 'coff']:
        if os.path.isfile(self.paml_output(alg)):
            results_exist.append(alg)
        self.dNdS['average'] = 'undetermined'
    
    ## for each of those results, grab the dNdS table
    tots = 0
    maxDS = []
    for algnr in results_exist:
        ## find the BEB sites from the paml results file
        with open(self.paml_output(algnr), 'r') as f:
            capture = False
            self.dNdS[algnr] = {'maxDS' : [], 'avDS' : [], 'maxDN' : [], 'avDN' : [], 'maxDNDS' : [], 'avDNDS' : []}
            for line in f.readlines():
                if 'Use runmode = -2 for ML pairwise comparison.)' in line:
                    capture = True
    
                elif 'Model 7: beta (10 categories)' in line:
                    capture = False
                    
                if capture:
                    if len(line) > 0 and "Use runmode" not in line:
                        line = line.split() # makes this into a list of all items
                        for ele in line:
                            if '(' in ele and len(ele) > 0:
                                dn = float(ele.split('(')[1])
                                self.dNdS[algnr]['avDN'].append(dn)
                            if ')' in ele and len(ele) > 0:
                                ds = float(ele.split(')')[0])
                                self.dNdS[algnr]['avDS'].append(ds)
                            if '(' not in ele and ')' not in ele and not any(c.isalpha() for c in ele) and '-2' not in str(ele) and '=' not in str(ele):
                                ## this next bit is to take care of an error in PAML where 
                                ## dN values = 0 caused the pairwise dN/dS to evaluate as 
                                ## -1
                                if '-1' in ele:
                                    ele = '0'
                                self.dNdS[algnr]['avDNDS'].append(float(ele))

        
        if debug:
            print('pairwise dnds vales for ' + algnr + ': ' + str(self.dNdS[algnr]))

        ## get the average dNdS for each aligner that was run
        if len(self.dNdS[algnr]['avDNDS']) > 0:
            ## get the max value for each
            self.dNdS[algnr]['maxDNDS'] = str(max(self.dNdS[algnr]['avDNDS']))
            th = 0
            for x in self.dNdS[algnr]['avDNDS']:
                th += float(x)
            
            self.dNdS[algnr]['avDNDS'] = str(th / len(self.dNdS[algnr]['avDNDS']))          
        else:
            self.dNdS[algnr]['avDNDS'] = 'NA'
        
        ## for dS
        if len(self.dNdS[algnr]['avDS']) > 0:
            self.dNdS[algnr]['maxDS'] = max(self.dNdS[algnr]['avDS'])
        else:
            self.dNdS[algnr]['maxDS'] = 'NA'

        if debug:
            print('Average dNdS for ' + algnr + ': ' + str(self.dNdS[algnr]['avDNDS']))

        ## add to the total
        if self.dNdS[algnr]['avDNDS'] != 'NA':
            tots += float(self.dNdS[algnr]['avDNDS'])
        if self.dNdS[algnr]['maxDS'] != 'NA':
            maxDS.append(self.dNdS[algnr]['maxDS'])
    
    ## remove anything that equals NA
    temp_list = []
    for algnr in results_exist:
        if self.dNdS[algnr]['avDNDS'] != 'NA':
            temp_list.append(algnr)
    results_exist = temp_list

    ## if none exist, set the output to NA
    if len(results_exist) < 1:
        self.dNdS['average'] = 'NA'
        self.dNdS['maxDS'] = 'NA'


    ## after, average all the aligners together
    if self.dNdS['average'] != 'NA':
        self.dNdS['average'] = tots / len(results_exist)
        self.dNdS['maxDS'] = max(maxDS)
    if debug:
        print('made it through dNdS')