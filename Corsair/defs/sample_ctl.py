#!/usr/bin/env python3

def sample_ctl():
    """
    Generates a sample control file in the current working directory
    """
    with open ("sample.ctl", 'w+') as f:
        f.write(
            "### CTL Example file ###\n\n" +
            "## this is where all the files get stored. This does NOT have to be where the reference files are located, but any created files will go here\n" +
            "project_directory:corsair/primates/\n\n" +
            "## reference CDS file, must be a cleaned fasta file. See README for more info.\n" +
            "reference_CDS:corsair/primates/reference_CDS.fasta\n\n" +            
            "## this where all the genomes are located\n" +
            "genome_directory:corsair/primates/genomes/\n\n" +
            "## clade tree (newick format). Names don't have to be 4 letters, but they MUST match the prefixes on genome files\n" +
            "tree:(((((Ptro,Ppan),Hsap),Ggor),((Mfas,Mmul),(Caty,Mleu))),Sbol);\n\n" +  
            "## reference species - in the tree listed above, this is the reference\n" +
            "ref_spec:Hsap\n\n" +
            "## minimum species count - this is the fewest species that PAML will be attempted with. min:3, max:species_count, default:0.7*species_count\n" +
            "minimum_species:3\n\n" +
            "## alignment threshold - this is the alignment threshold for identifying genes. min:0.7, max:1, default:0.95\n" +
            "alignment_threshold:0.95\n\n"
        )
