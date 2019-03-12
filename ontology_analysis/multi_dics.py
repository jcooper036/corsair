#! /usr/bin/env python3

import sys
from copy import deepcopy

hit_list = sys.argv[1]
kegg_db = '/Users/Jacob/corsair/ontology_analysis/hsa00001_formated.txt'

A = {}
B = {}
C = {}

with open(kegg_db, 'r') as f:
    
    hold_A = ''
    hold_B = ''
    hold_C = ''
    set_A = set()
    set_B = set()
    set_C = set()

    for line in f:
        line = line.strip()

        ## headers        
        if line.split(',')[0] == 'A':
            if hold_A and set_A:
                A[hold_A] = deepcopy(set_A)
            hold_A = line.split(',')[2]
            set_A = set()
        elif line.split(',')[0] == 'B':
            if hold_B and set_B:
                B[hold_B] = deepcopy(set_B)
            hold_B = line.split(',')[2]
            set_B = set()
        elif line.split(',')[0] == 'C':
            if hold_C and set_C:
                C[hold_C] = deepcopy(set_C)
            hold_C = line.split(',')[2]
            set_C = set()
            
        ## genes
        elif line.split(',')[0] == 'D':
            set_A.add(line.split(',')[2])
            set_B.add(line.split(',')[2])
            set_C.add(line.split(',')[2])

gene_list = []
with open(hit_list, 'r') as f:
    for line in f:
        line = line.strip()
        gene_list.append(line)



def count_gene(gene_list, kegg_dict):
    count_dict = {}
    for key in kegg_dict:
        count_dict[key] = []
        for gene in gene_list:
            if gene in kegg_dict[key]:
                count_dict[key].append(gene)
    
    return count_dict


def kegg_ratio(gene_list, kegg_dict, count_dict):
    ratio_dict = {}
    for key in kegg_dict:
        kegg_r = (len(count_dict[key]) / len(gene_list)) / len(kegg_dict[key])
        ratio_dict[key] = kegg_r
    return ratio_dict


count_A = count_gene(gene_list, A)
count_B = count_gene(gene_list, B)
count_C = count_gene(gene_list, C)

ratio_A = kegg_ratio(gene_list, A, count_A)
ratio_B = kegg_ratio(gene_list, B, count_B)
ratio_C = kegg_ratio(gene_list, C, count_C)


### print
print('######### A ##########')
for key in A:
    print("{} | genes:{}/{} | kegg ratio:{}".format(key, len(count_A[key]), len(A[key]), ratio_A[key]))
print('\n')

print('######### B ##########')
for key in B:
    print("{} | genes:{}/{} | kegg ratio:{}".format(key, len(count_B[key]), len(B[key]), ratio_B[key]))
print('\n')

print('######### C ##########')
for key in C:
    print("{} | genes:{}/{} | kegg ratio:{}".format(key, len(count_C[key]), len(C[key]), ratio_C[key]))
print('\n')
