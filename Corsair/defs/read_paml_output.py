#!/usr/bin/env python3

import Corsair as cor
from scipy.stats import chi2

def read_paml_output(ctl, iso_name):
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
    results = {
        'logLvals' : {},
        '2delta' : False,
        'pval' : False, # absolutely needed, other code depends on it
        'beb_sites_raw' : False
    }
    with open(file, 'r') as f:
        logflag = False
        for line in f.readlines():
            line = line.strip()
            
            ## log likelihood values
            if 'Model' in line:
                print(line)
                logflag = True
                line = line.split('Model')[1].split(':')[0]
                if line:
                    model = 'M' + str(line.split(' ')[1])
            if logflag and "lnL" in line:
                line = line.replace(' ','').split(':')[3].split('+')[0]
                results['logLvals'][model] = line
    
        

    ## compute pvalues of from log likelihood scores, 2 degrees of freedom
    results['2delta'] = 2 * abs(float(results['logLvals']['M7']) - float(results['logLvals']['M8']))
    results['pval'] = chi2.sf(results['2delta'], 2)

    return results

# def loglike_get(self):
# ## scans the results file for a log-likelihood values
# ## stores properties in a dictionary:
# ## llvals: dictionary of aligner_model : log-likelihood
# self.llvals = {}
# for algr in ('clus','coff','mus'):
#         try:
#         if os.path.exists(self.paml_output(algr)):
#                 with open(self.paml_output(algr), 'r') as f:
#                         flag = False
#                         for line in f.readlines():

#                         if 'Model' in line: ## determines the model
#                                 flag = True
#                                 line = [x.strip() for x in line]
#                                 line = ''.join(line)
#                                 line = line.split('Model')[1]
#                                 line = line.split(':')[0]
#                                 if line: # takes care of some blank lines that get returned
#                                 model = algr + '_M' + line

#                         if flag and 'lnL' in line: ## gets the log-likelihood value
#                                 line = [x.strip() for x in line]
#                                 line = ''.join(line)
#                                 line = line.split('):')[1]
#                                 lnl = line.split('+')[0]
#                                 self.llvals[model] = lnl
#                                 flag = False
#         else:
#                 if debug:
#                 print('PAML output file for ' + algr + ' not detected.')
#         except:
#         print('Can not open paml file')

#         if 'clus_M7_M8_p' not in self.llvals:
#         self.llvals['clus_M7_M8_p'] = 37767



#         if os.path.exists(self.p8ml_output(algr)):
#         with open(self.p8ml_output(algr), 'r') as f:
#                 for line in f.readlines():
#                 if 'lnL' in line: ## gets the log-likelihood value
#                         model = algr + '_M8a'
#                         line = [x.strip() for x in line]
#                         line = ''.join(line)
#                         line = line.split('):')[1]
#                         lnl = line.split('+')[0]
#                         self.llvals[model] = lnl
#                         flag = False

#         ## calculates the 2Delta log-likelihood
#         ## calculates a p-value

#         M7 = algr + '_M7'
#         M8 = algr + '_M8'
#         M8a = algr + '_M8a'
#         M7_M8_2del = algr + '_M7_M8_2del'
#         M7_M8_p = algr + '_M7_M8_p'
#         M8_M8a_2del = algr + '_M8_M8a_2del'
#         M8_M8a_p = algr + '_M8_M8a_p'

#         ## see if there is a M7-M8 comparison
#         try:
#         self.llvals[M7_M8_2del] = 2 * ( abs( float(self.llvals[M7]) - float(self.llvals[M8]) ))
#         self.llvals[M7_M8_p] = chi2.sf ( self.llvals[M7_M8_2del] , 2 )
#         except:
#         pass #@ And a message here I think.

#         ## see if there is a M8-M8a comparison
#         try:
#         self.llvals[M8_M8a_2del] = 2 * ( abs( float(self.llvals[M8]) - float(self.llvals[M8a]) ))
#         self.llvals[M8_M8a_p] = chi2.sf ( self.llvals[M8_M8a_2del] , 1 )
#         except:
#         pass #@ Message about no significant p-values.

# ## determine the current highest p-value
# d = {}
# for algr in ('clus','mus','coff'):
#         try:
#         M7_M8_p = algr + '_M7_M8_p'
#         d[algr] = self.llvals[M7_M8_p] #@ debug step could be to print these.
#         except:
#         pass

# # self.llvals['max_p'] = keywithmaxval(d)
# if d:
#         self.llvals['max_p'] = max(d, key=d.get)
# else:
#         pass