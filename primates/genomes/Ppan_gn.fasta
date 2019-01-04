### CTL Example file ###

## this is where all the files get stored. This does NOT have to be where the reference files are located, but any created files will go here
project_directory:corsair/primates/

## reference CDS file, must be a cleaned fasta file. See README for more info.
reference_CDS:corsair/primates/reference_CDS.fasta

## this where all the genomes are located
genome_directory:corsair/primates/genomes/

## clade tree (newick format). Names don't have to be 4 letters, but they MUST match the prefixes on genome files
tree:(((((Ptro,Ppan),Hsap),Ggor),((Mfas,Mmul),(Caty,Mleu))),Sbol);

## reference species - in the tree listed above, this is the reference
ref_spec:Hsap

## minimum species count - this is the fewest species that PAML will be attempted with. min:3, max:species_count, default:0.7*species_count
minimum_species:3

## alignment threshold - this is the alignment threshold for identifying genes. min:0.7, max:1, default:0.95
alignment_threshold:0.95

