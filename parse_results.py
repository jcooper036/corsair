#!/usr/bin/env python3
"""
Usage: parse_results.py <results_file>

Arguments:
    <results_file>    Path to the results file

"""
#imports
import docopt

## Initialize docopt
if __name__ == '__main__':

    try:
        arguments = docopt.docopt(__doc__)
        results_file = str(arguments['<results_file>'])
    except docopt.DocoptExit as e:
        print(e)


## read results file, save lines that have ceratin properties
def read_resutls_table(file):
    """The results table is a tsv. Here is the headers:
    Gene[0] Gene_Ensembl_ID[1] Isoform_Ensembl_ID[2]
    Species[3] Reference_length(aa)[4] ...
    Returns a list of the lines that we want to keep
    """
    retlist = []
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            line = line.split('\t')
            keep = False
            if "Species" in line:
                keep = True
            else:
                if int(line[3]) >= 3:
                    keep = True
            if keep:
                line = '\t'.join(line)
                retlist.append(line)
    return retlist


def print_results_table(results, results_file):
    """Takes the new results and the old name of the file
    to print a new results file based on the name of the
    old one"""
    new_file = results_file.split('.txt')[0] + '.parsed.txt'
    with open(new_file, 'w') as f:
        for line in results:
            f.write(line + '\n')


results = read_resutls_table(results_file)
print_results_table(results, results_file)
