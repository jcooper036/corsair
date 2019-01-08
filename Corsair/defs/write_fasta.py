#!/usr/bin/env python3

def write_fasta(file, fasta_dict):
    """Writes a protein sequence to a file in fasta format"""
    with open(file, 'w+') as f:
        for key in fasta_dict:
            f.write('>' + key + '\n')
            total_length = len(fasta_dict[key])
            total_count = 0
            sixty_count = 0
            while total_count < total_length:
                if sixty_count < 60:
                    f.write(fasta_dict[key][total_count])
                    total_count += 1
                    sixty_count += 1
                else:
                    f.write('\n')
                    sixty_count = 0
            f.write('\n')
        