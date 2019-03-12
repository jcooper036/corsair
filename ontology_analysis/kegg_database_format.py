#! /usr/bin/env python3

import sys
kegg_file = sys.argv[1]
outfile = kegg_file.split('.')[0] + "_formated.txt"


outlist = []
with open(kegg_file, 'r') as f:
    for line in f:
        line = line.lstrip()
        line = line.rstrip()
        
        ## some line contain only one character, and we don't want those
        ## some also have characters that we don't want
        if not (len(line) > 1):
            continue
        if not any(x in line[0] for x in ['A','B','C','D']):
            continue
        
        ## unfortuneately these need to be handled differently

        if line[0] == 'A':
            outlist.append(
                ['A', line.split(' ')[0][1:], " ".join(line.split(' ')[1:])]
            )
        
        if line[0] == 'B':
            outlist.append(
                ['B', line.split()[1], " ".join(line.split()[2:])]
            )
        
        if line[0] == 'C':
            if '[' in line:
                line = line.split(' [')[0]
            outlist.append(
                ['C', line.split()[1], " ".join(line.split()[2:])]
            )            
        
        if line[0] == 'D':
            ls = line.split()[0:3]
            outlist.append(
                [ls[0], ls[1], ls[2].replace(';', '')]
            )

with open(outfile, 'w') as f:
    for line in outlist:
        line = [str(x) for x in line]
        line = ','.join(line)
        f.write(line + '\n')