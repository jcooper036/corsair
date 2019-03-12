#! /usr/bin/env python3

import sys

infile = sys.argv[1]
outfile = '/Users/Jacob/corsair/ontology_analysis/hit_lists/' + infile.split('_')[2] + '_hits.txt'

hit_genes = []

with open(infile, 'r') as f:
    for line in f:
        line = line.strip('\n')
        line = line.split('\t')
        if len(line) >= 15:
            if line[14] and not ('gene' in [x.lower() for x in line]):
                if float(line[14]) < 0.05:
                    hit_genes.append(line[0])

with open(outfile, 'w') as f:
    for gene in hit_genes:
        f.write(gene + '\n')