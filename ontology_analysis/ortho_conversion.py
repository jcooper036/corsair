#! /usr/bin/env python3

import sys

"""Converts hits from one species (ENSEMBL IDS) to gene names from humans"""

hit_file = sys.argv[1]
conversion_table = sys.argv[2]
name_table = sys.argv[3]
outfile = '/Users/Jacob/corsair/ontology_analysis/converted_hit_lists/' + hit_file.split('/')[-1].split('_')[0] + '_converted_hits.txt'

## get the esemble IDS
hits = []
with open(hit_file, 'r') as f:
    for line in f:
        line = line.strip()
        hits.append(line.split('.')[0])

## sort through the hit list
idconvert = {}
with open(conversion_table, 'r') as f:
    for line in f:
        line = line.strip()
        line = line.split(',')
        if len(line) > 3:
            if line[3] == 'ortholog_one2one':
                idconvert[line[2]] = line[0]

## make the conversion to human ENSEMBL
human_ensb = []
for gene in hits:
    if gene in idconvert:
        human_ensb.append(idconvert[gene])

## load in the human gene names and ensb
human_names = {}
with open(name_table, 'r') as f:
    for line in f:
        line = line.strip()
        line = line.split(',')
        human_names[line[1].split('.')[0]] = line[0].replace('-a-', '.')

## get the names
genes = []
for gene in human_ensb:
    if gene in human_names:
        genes.append(human_names[gene])

print("{} hits input, {} hits converted".format(len(hits), len(genes)))

with open(outfile, 'w') as f:
    for gene in genes:
        f.write(gene + '\n')


