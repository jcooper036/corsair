#! /usr/bin/env python3

"""Usage: corsair.py [options] [--process=N] --ctl=N FILE

Arguments:
    FILE  File that contains gene names

Options:
    --fresh         All data will be new (default: old data used)
    --isoform=N     All isoforms will be used [default: long]
    --process=N     Specifies which process to run [default: all]
    --ctl=N         Indicates the control file name located in the control folder
    --count=N       Changes the number of species necessary for PAML
    --cores=N       Specifices the number of cores to be used
    --debug         Prints more information
"""

## where the genomes are located
genomes_db = '/Volumes/Jacob_2TB_storage/worm_genomes'

#Tools for system
import docopt
import os
import shutil
import subprocess
from joblib import Parallel, delayed
import multiprocessing
import pickle

#Tools for working with sequence
from Bio import SeqIO
from io import StringIO

#Tools for PAML
from Bio import Phylo
from scipy.stats import chi2

## Initialize docopt
if __name__ == '__main__':

    try:
        arguments = docopt.docopt(__doc__)
        fresh = arguments['--fresh']
        alliso = arguments['--isoform']
        debug = arguments['--debug']
        process = arguments['--process']
        ctl = arguments['--ctl']
        speccount = arguments['--count']
        cores = arguments['--cores']
        gene_file = str(arguments['FILE'])

    except docopt.DocoptExit as e:
        print(e)

################################################################
################ definintions
################################################################

def trans(sequence):
    # takes a nucleotide string and returns the translation as an amino acid string

    # codon table
    transtable = {
        'agg' : 'R',
        'aga' : 'R',
        'agc' : 'S',
        'agt' : 'S',
        'aag' : 'K',
        'aaa' : 'K',
        'aac' : 'N',
        'aat' : 'N',
        'aca' : 'T',
        'acc' : 'T',
        'acg' : 'T',
        'act' : 'T',
        'atg' : 'M',
        'ata' : 'I',
        'atc' : 'I',
        'att' : 'I',
        'cgg' : 'R',
        'cga' : 'R',
        'cgc' : 'R',
        'cgt' : 'R',
        'cag' : 'Q',
        'caa' : 'Q',
        'cac' : 'H',
        'cat' : 'H',
        'ccg' : 'P',
        'cca' : 'P',
        'ccc' : 'P',
        'cct' : 'P',
        'ctg' : 'L',
        'cta' : 'L',
        'ctc' : 'L',
        'ctt' : 'L',
        'tgg' : 'W',
        'tga' : 'X',
        'tgc' : 'C',
        'tgt' : 'C',
        'tag' : 'X',
        'taa' : 'X',
        'tac' : 'Y',
        'tat' : 'Y',
        'tcg' : 'S',
        'tca' : 'S',
        'tcc' : 'S',
        'tct' : 'S',
        'ttg' : 'L',
        'tta' : 'L',
        'ttc' : 'F',
        'ttt' : 'F',
        'ggg' : 'G',
        'gga' : 'G',
        'ggc' : 'G',
        'ggt' : 'G',
        'gag' : 'E',
        'gaa' : 'E',
        'gac' : 'D',
        'gat' : 'D',
        'gcg' : 'A',
        'gca' : 'A',
        'gcc' : 'A',
        'gct' : 'A',
        'gtg' : 'V',
        'gta' : 'V',
        'gtc' : 'V',
        'gtt' : 'V'}

    # this is a counting variable
    k = 3

    # this the blank string for the amino acid sequence
    aaseq = ''

    # stops when it gets to the end of the sequence=
    while k <= len(sequence):
        # try to add a codon. will reject because it won't find the key if the variable is blank
        try:
            aaseq += transtable[sequence[(k-3):k].lower()]
        except:
            aaseq += '' ## adds a blank space if no translation is possible

        # count up 3
        k = k + 3

    # return value is a string of amino acids
    return aaseq

def fasta_read(file):
    ## reads a file, outputs a dicitonary
    wholefile = {}
    with open(file,'r') as f:
        for line in f.readlines():
            if ">" in line:
                keyname = [x.strip() for x in line]
                keyname = ''.join(keyname)
                keyname = keyname.split('>')[1]
                wholefile[keyname] = []
            elif '>' not in line:
                line = [x.strip() for x in line]
                line = ''.join(line)
                wholefile[keyname].append(line)
    for key in wholefile:
        wholefile[key] = ''.join(wholefile[key])

    return wholefile

def linkisogene(temp_isodic, genedic):
    # determine all the unique parts of the isoforms
    # requires isoform dictionary, isoform list, gene list
    # returns two dictionaries: updated isoform dictionary, and a gene dictionary

    ## gives each isoform the .parent_gene property
    for x in temp_isodic:
        temp_isodic[x].link_parent(temp_isodic[x].name.split('_')[0])

    # add all the isoforms to the list of isoforms within the gene. This is a little slow if there are only a few genes on the list, since it is looking through every isoform
    for isoforms in temp_isodic:
        if temp_isodic[isoforms].parent_gene in genedic:
            genedic[temp_isodic[isoforms].parent_gene].linkiso(temp_isodic[isoforms].name)
    return temp_isodic, genedic

def keywithmaxval(d):
     """ a) create a list of the dict's keys and values;
         b) return the key with the max value"""
     v=list(d.values())
     k=list(d.keys())
     return k[v.index(max(v))]

def findunique(temp_isodic, genedic, genex):
    # takes the dictionary of isoforms and the gene dictionary
    # returns the dictionary of isoforms and the gene dictionary with unique isoforms

    logg = {}
    for uniso in genedic[genex].all_isoforms:
        logg[uniso] = temp_isodic[uniso].ref_aa


    # get the name of the longest isoform for each gene, or skip the rest of the function if empty
    if logg: #This only runs if the isodict contained the gene name

        maxseq = sorted(sorted(logg.items(), key=lambda x: x[0]),key=lambda x: len(x[1]),reverse=True)[0][0]
        genedic[genex].longiso(maxseq)
        genedic[genex].uniqiso(maxseq)

        ## this loop checks for unique isoforms
        for testiso in genedic[genex].all_isoforms:
            dashes = True
            for knownisos in genedic[genex].uniq_isos:
                if temp_isodic[testiso].ref_aa in temp_isodic[knownisos].ref_aa: ## this is actually a really good line, because it just does a simple check (though not complete) to make sure that one is not just a subset of the other, and without alignment this saves a massive amount of time.
                    dashes = False
            ## save the sequence if not contained in other isoform
            if dashes:
                genedic[genex].uniqiso(testiso)

        if debug:
            print("\n",genedic[genex].all_isoforms)
            print(genedic[genex].uniq_isos, "\n")

    else: #If the gene name was not found, then nothing else is run on it.

        genedic = ''

    return genedic

def blexon(isodic, iso, genomes):
    ## input is the isoform dictionary, an isoform name, gene dictionary, and the genomes list
    ## output is a dictionary with species : sequence values

    blast_dic = {}
    current_gene = isodic[iso].parent_gene

    ##@ load up the pickle file for that gene
    with open('genes/' + current_gene + '/' + current_gene + '_save_file.pkl', 'rb') as input:
        genedic[current_gene] = pickle.load(input)


    ## check if the translation file exists or not, and make it if not
    if debug:
        try:
            print ('Ref aa seq for ' + iso + ' is:\n' + isodic[iso].ref_aa + '\n')
        except:
            print (isodic[iso] + ' does not have a translated aa sequence')

    try:
        with open(isodic[iso].prot_file(), "w") as f:
            f.write('>' + iso + '\n')
            f.write(isodic[iso].ref_aa)
    except:
        print ('Did not write protein file for ' + iso)

    ## check if the sequences is already there, get the scaffold_id if not
    for db in genomes:
        if os.path.exists(isodic[iso].iso_files() + db + '.txt'):
            ## check if sequences are already there
            print(db + ' CDS for ' + iso + ' already obtained.')
        elif db in genedic[current_gene].scaffolds.keys():
            ## check if there is a scaffold recorded for this genes
            print('Previously found the scaffold for ' + current_gene + ':' + genedic[current_gene].scaffolds[db])
        else:
            if debug:
                print('Blasting: ' + db + ' for ' + iso + '.')
            var=str(subprocess.check_output('tblastn -outfmt "6 sseqid" -query "' + isodic[iso].prot_file() + '" -db ' + genomes_db + '/'+ db + '_genome.fasta -max_target_seqs 1| head -n 1',shell=True))
            # at this point, this dictionary contains species : scaffold_id pairs
            genedic[current_gene].add_scaffold(db, (var.replace('b','').replace('\\n','').replace('\'','').replace('|','')))
    if debug:
        print(str(genedic[current_gene].scaffolds))


    for db, scaffold in genedic[current_gene].scaffolds.items():
        #This for loop iterates through the list of scaffolds connected to species and exonerates the gene from the temp scaffold.
        #prints scaffold name
        if debug:
            print(db + ' scaffold for ' + iso + ': ' + scaffold)

        #parses through genome fasta for hit scaffold, writes to temp named by scaffold and gene
        subprocess.check_output('seqret -sequence ' + genomes_db + '/' + db + '_genome.fasta:' + scaffold + ' -outseq temp/' + scaffold + '_' + iso + '.fasta',shell=True)


        #exonerates through scaffold using input protein
        if os.path.exists('temp/' + scaffold + '_' + iso + '.fasta'):
            subprocess.check_output('exonerate --model protein2genome --query ' + isodic[iso].prot_file() + ' --target temp/' + scaffold + '_' + iso + '.fasta --ryo "forcebeams \n>' + db + '\n%tcs\n" --bestn 1 >' + isodic[iso].iso_files() + db + '.txt' ,shell=True)

            #removes temporary scaffold (mitigating measure for parallel)
            os.remove('temp/' + scaffold + '_' + iso + '.fasta')

        else:
            subprocess.check_output('exonerate --model protein2genome --query ' + isodic[iso].prot_file() + ' --target ' + genomes_db + '/' + db + '_genome.fasta --ryo "forcebeams \n>' + db + '\n%tcs\n" --bestn 1 >' + isodic[iso].iso_files() + db + '.txt' ,shell=True)
    #parses exonerate output file, writes alignment to the blast dictionary
    for db in genomes:
        blast_dic[db] = []
        with open(isodic[iso].iso_files() + db + '.txt') as infile:
            copy = False
            for line in infile:
                if line.strip() == 'forcebeams':
                    copy = True
                elif line.strip() == '-- completed exonerate analysis':
                    copy = False
                elif line.strip() == 'C4 Alignment:':
                    copy = False
                elif copy:
                    if '>' not in line:
                        # line = [x.strip() for x in line]
                        # line = ''.join(line) #@ Are these two lines necessary if file is being read line by line? Could just use str(line.strip())
                        blast_dic[db].append(line.strip())
    for key in blast_dic:
        blast_dic[key] = ''.join(blast_dic[key])

    if debug:
        print(blast_dic)

    #percent length required, based on length of gene
    if len(isodic[iso].ref_nt) <= 2400:
        plreq = (len(isodic[iso].ref_nt) / 40000) + .89
    elif len(isodic[iso].ref_nt) > 2400:
        plreq = .949

    if debug:
        print('\nPercent length required: ' + str(plreq*100) + '\n')

    ## add sequences above the cut-off, and then reset the dictionary
    temp_dic = {}
    for db in blast_dic:
        if len(blast_dic[db]) > plreq * len(isodic[iso].ref_nt):
            temp_dic[db] = blast_dic[db]
        else:
            if debug:
                print('Removed ' + db + ' from ' + iso + ' sequence list')

    ## add the reference sequence
    temp_dic[clade] = isodic[iso].ref_nt

    ## write a CDS outfile. needed for tree building
    with open(isodic[iso].CDS_file(), 'w') as f:
        for key in temp_dic:
            f.write('>' + key + '\n')
            f.write(temp_dic[key] + '\n')

    blast_dic = temp_dic

    ##@ save the gene object
    with open('genes/' + current_gene + '/' + current_gene + '_save_file.pkl', 'wb') as output:
        pickle.dump (genedic[current_gene], output, pickle.HIGHEST_PROTOCOL)

    return blast_dic

def align(isodic, iso, aligner): #to be used below after species count check

    infile = 'temp/' + iso + '_' + aligner + '.fasta'
    outfile = isodic[iso].alignment_file(aligner)

    if isodic[iso].blast_dic:
        isodic[iso].blast_trans()

        if debug:
            print(isodic[iso].blast_prot) #@ maybe debug overkill here, especially when running in parallel

        if not os.path.exists(outfile):
            with open(infile, 'w') as f:
                for key in isodic[iso].blast_prot:
                    f.write('>' + str(key) + '\n')
                    f.write(str(isodic[iso].blast_prot[key]) + '\n')

            if aligner=='clus':
                subprocess.check_output('clustalo -i ' + infile +' -o ' + outfile + ' --force',shell=True)
                print('Sequences for ' + iso + ' aligned with ' + aligner + '.')
            elif aligner=='coff':
                coffee = subprocess.check_output('t_coffee ' + infile + ' -outfile ' + outfile + ' -multi_core -quiet -output=fasta',shell=True)
                try:
                    os.remove(iso + '_coff.dnd')
                except:
                    pass
                print('Sequences for ' + iso + ' aligned with ' + aligner + '.')
            elif aligner=='mus':
                subprocess.check_output('muscle -quiet -in ' + infile + ' -out ' + outfile + '',shell=True)
                print('Sequences for ' + iso + ' aligned with ' + aligner + '.')

            ##remove temporary input file
            os.remove(infile)

    ## return just the name of the outfile to be read by the run commands
    return outfile

def back_translate(aa_alignment, nuc_input, iso, aln, isodic):
    # takes a dictionary that contains amino acid strings, that are already
    # aligned, where stop codons are marked with X or *, that are all the same
    # length, and reduces them to
    # the minimum amino continuous amino acid sequence across all species

    # make sure input is divisable by 3
    for key in nuc_input:
        if not (len(nuc_input[key]) % 3) == 0:
            print('\nWARNING: The nucleotide input for back translation for ' + key + ' is not divisable by 3. Please check.\n')

    minaa = {}
    for key in aa_alignment:
        cond_seq = ""
        for i in range(0, len(aa_alignment[key])):
            if not any(aa_alignment[nkey][i] == "*" or aa_alignment[nkey][i] == "X" or aa_alignment[nkey][i] == "-" for nkey in aa_alignment):
                cond_seq += aa_alignment[key][i]
        minaa[key] = cond_seq

    minnuc = {}
    for key in minaa:
        k = 0
        cond_nuc = nuc_input[key]
        while k < len(minaa[key]):
            if not minaa[key][k] == trans(cond_nuc[(k*3):((k*3)+3)]): #@ translating codon by codon, which means that if there is an error message in the translation def then it will pop up every time
                cond_nuc = cond_nuc[:(k*3)] + cond_nuc[(k*3)+3:]
            else:
                k += 1
        else:
            cond_nuc = cond_nuc[:(k*3)]

        minnuc[key] = cond_nuc

    if debug:
        print(str(minnuc))

    ## write a file for the output
    with open(isodic[iso].paml_file(aln), 'w') as f:
        f.write('\t' + str(len(minnuc.keys())) + '\t' + str(len(minnuc[clade])) + '\n') ## changed to 'clade' variable
        for key in minnuc:
            f.write(str(key) + '\n')
            f.write(str(minnuc[key]) + '\n')

    ## and also return the dictionary with the back-translation
    return minnuc

def treebuild(iso, isodic, genomes, newick, speccount):
    #checks if CDS exists, otherwise exits
    if os.path.exists(isodic[iso].CDS_file()):

        #reads the cds file to a new variable (to see what species we have)
        cds = SeqIO.to_dict(SeqIO.parse(isodic[iso].CDS_file(),'fasta'))
        print('Number of sequences for ' + iso + ' : ' + str(len(cds)))

        #works for all groups, dependent on the clade/genomes combination
        handle = StringIO(newick)
        tree = Phylo.read(handle,'newick') #loads tree

        #if species in genomes but not in CDS, then removes from the tree
        for db in genomes:
            if db not in cds:
                tree.prune(db)
        Phylo.write(tree, 'temp/' + iso + '_tree.txt', 'newick') #writes to temp file

        #If there are enough species: The phylo.write function adds branch lengths. New file without them is put into PAML folder.
        if len(cds) >= speccount:
            with open ('temp/' + iso + '_tree.txt','r') as infile, open(isodic[iso].tree_file(),'w') as outfile:
                outfile.write(infile.read().replace(':0.00000',''))
            print('Tree for ' + iso +': built.')

        else:
            print('Not enough species to build tree for ' + iso + '.')
            with open('results.txt','a') as f:
                f.write('\n' + str(iso) + '\tNot enough species for PAML: ' + str(len(cds)))
        os.remove('temp/' + iso + '_tree.txt') #removes temp file

def paml(isodic, iso, aligner):

    ## set up a few checks before trying
    check_input = False
    check_tree = False
    check_output_missing = False
    get_output = False

    ## check for the sequence alignment file
    if os.path.exists(isodic[iso].paml_file(aligner)):
        check_input = True
    else:
        print('There is no sequence input file for ' + iso + ' for ' + aligner )

    ## check for the tree file
    if os.path.exists(isodic[iso].tree_file()):
        check_tree = True
    else:
        print('There is no tree file for ' + iso)

    ## check for the output file, and if it is there, does it contain what it is supposed to?
    if os.path.exists(isodic[iso].paml_output(aligner)):
        isodic[iso].loglike_get()
        key = aligner + '_M8'
        try:
            isodic[iso].llvals[key]
        except:
            os.remove(isodic[iso].paml_output(aligner))
            check_output_missing = True
            print('Previous run for ' + iso + ' was no good, re-running.')
    else:
        check_output_missing = True

    ## use those conditions to actually run PAML or not. will turn on the 'get_output' switch if it runs
    if check_input and check_tree and check_output_missing:
        print ('Running PAML for ' + iso + ' with ' + aligner)

        #change diretories to gene folder
        os.chdir(isodic[iso].iso_files())

        #Write/rewrite control file
        with open('codeml.ctl','w') as f:
            f.write('seqfile = ' + iso + '_' + aligner +'.paml \n')
            f.write('treefile = ' + iso + '_tree.txt \n')
            f.write('outfile = ' + iso + '_' + aligner + '_PAML_out_full.txt \n')
            f.write('noisy = 3 \n')
            f.write('verbose = 1 \n')
            f.write('runmode = 0 \n')
            f.write('seqtype = 1 \n')
            f.write('CodonFreq = 2 \n')
            f.write('ndata = 1 \n')
            f.write('clock = 0  \n')
            f.write('aaDist = 0 \n')
            f.write('model = 0 \n')
            f.write('NSsites = 7 8 \n')
            f.write('icode = 0 \n')
            f.write('Mgene = 0 \n')
            f.write('fix_kappa = 0 \n')
            f.write('kappa = 2 \n')
            f.write('fix_omega = 0 \n')
            f.write('omega = 0.4 \n')
            f.write('fix_alpha = 1 \n')
            f.write('alpha = 0 \n')
            f.write('Malpha = 0 \n')
            f.write('ncatG = 8 \n')
            f.write('getSE = 0 \n')
            f.write('RateAncestor = 1 \n')
            f.write('Small_Diff = .5e-6 \n')
            f.write('cleandata = 1 \n')
            f.write('method = 0 \n')

        #Run PAML
        try:
            test = subprocess.check_output('codeml',shell=True)

        except:
            pass
        #remove all the files codeml leaves behind
        try:
            for f in ['2NG.dN','2NG.dS','2NG.t','4fold.nuc','lnf','rst','rst1','rub']:
                os.remove(f)
        except:
            pass

        #change directories back to Corsair
        os.chdir('../../../..')

        if debug:
            print(os.getcwd())

        print('PAML for ' + iso + '(' + aligner + '): run.')

        get_output = True #@ what is this variable for?

    elif check_input and check_tree and not check_output_missing:
        print('PAML was already run for ' + iso + ' with ' + aligner)
        get_output = True
    else:
        print ('PAML was not run for ' + iso)

def p8ml(isodic, iso, aligner):

    ## set up a few checks before trying
    check_input = False
    check_tree = False
    check_output_missing = False
    get_output = False

    ## check for the sequence alignment file
    if os.path.exists(isodic[iso].paml_file(aligner)):
        check_input = True
    else:
        print('There is no sequence input file for ' + iso + ' for ' + aligner )

    ## check for the tree file
    if os.path.exists(isodic[iso].tree_file()):
        check_tree = True
    else:
        print('There is no tree file for ' + iso)

    ## check for the output file, and if it is there, does it contain what it is supposed to?
    if os.path.exists(isodic[iso].p8ml_output(aligner)):
        isodic[iso].loglike_get()
        key = aligner + '_M8'
        try:
            isodic[iso].llvals[key]
        except:
            os.remove(isodic[iso].p8ml_output(aligner))
            check_output_missing = True
            print('Previous run for ' + iso + ' was no good, re-running.')
    else:
        check_output_missing = True

    ## use those conditions to actually run PAML or not. will turn on the 'get_output' switch if it runs
    if check_input and check_tree and check_output_missing:
        print ('Running PAML for ' + iso + ' with ' + aligner)
        #change diretories to gene folder
        os.chdir(isodic[iso].iso_files())
        #Write/rewrite control file
        with open('codeml.ctl','w') as f:
            f.write('seqfile = ' + iso + '_' + aligner +'.paml \n')
            f.write('treefile = ' + iso + '_tree.txt \n')
            f.write('outfile = ' + iso + '_' + aligner + '_M8a_PAML_out_full.txt \n')
            f.write('noisy = 3 \n')
            f.write('verbose = 1 \n')
            f.write('runmode = 0 \n')
            f.write('seqtype = 1 \n')
            f.write('CodonFreq = 2 \n')
            f.write('ndata = 1 \n')
            f.write('clock = 0  \n')
            f.write('aaDist = 0 \n')
            f.write('model = 0 \n')
            f.write('NSsites = 8 \n')
            f.write('icode = 0 \n')
            f.write('Mgene = 0 \n')
            f.write('fix_kappa = 0 \n')
            f.write('kappa = 2 \n')
            f.write('fix_omega = 1 \n')
            f.write('omega = 1 \n')
            f.write('fix_alpha = 1 \n')
            f.write('alpha = 0 \n')
            f.write('Malpha = 0 \n')
            f.write('ncatG = 8 \n')
            f.write('getSE = 0 \n')
            f.write('RateAncestor = 1 \n')
            f.write('Small_Diff = .5e-6 \n')
            f.write('cleandata = 1 \n')
            f.write('method = 0 \n')

        #Run PAML
        try:
            test = subprocess.check_output('codeml',shell=True)

        except:
            pass

        #remove all the files codeml leaves behind
        try:
            for f in ['2NG.dN','2NG.dS','2NG.t','4fold.nuc','lnf','rst','rst1','rub']:
                os.remove(f)
        except:
            pass
        #print significance to results file with aligner variable

        #change directories back to Corsair
        os.chdir('../../../..')
        print('PAML for ' + iso + '(' + aligner + '): run.')
        get_output = True  #@what is this variable for?
    elif check_input and check_tree and not check_output_missing:
        print('PAML was already run for ' + iso + ' with ' + aligner)
        get_output = True
    else:
        print ('PAML was not run for ' + iso)

def isoform_get(gene_name, genedic, isodic, alliso, fresh):
    ## input: gene name and two boolean operators to tell it what to do
    ## output: dictionary with the gene as an object and dictionary of the isoforms that need to be run

    ## attempt to make gene folder
    try:
        os.makedirs('genes/' + gene_name)
    except:
        pass

    ## adds the gene as a class object to the genedic dictionary
    genedic[gene_name] = gene(gene_name) ## make an object out of the gene name

    ## make a string to search to prevent inclusion
    search_string = str(gene_name)  #@ this might be clade specific
    
    ## make the temp isodic from the isoforms with the gene name
    temp_isodic = {x : isoform(x) for x in cds if search_string in x}

    ## check for actual sequence info
    if temp_isodic:
        ## give nt and aa sequence to those isoforms
        for iso in temp_isodic:
            temp_isodic[iso].ref_nt(cds[iso])
            temp_isodic[iso].translate()

        ## links the isoforms and genes together
        temp_isodic, genedic = linkisogene (temp_isodic, genedic)

        ## identifyies the longest and unique isoforms for each gene
        genedic = findunique (temp_isodic, genedic, gene_name)


        ## make the isodic based on if we need all the unique isoforms or not
        if genedic: #This process only runs if one or more sequences were obtained for the gene
            if alliso == 'unique':
                for x in genedic[gene_name].uniq_isos:
                    isodic[x] = temp_isodic[x]
            elif alliso == 'all':
                for x in temp_isodic:
                    isodic[x] = temp_isodic[x]
            else:
                isodic[genedic[gene_name].long_iso] = temp_isodic[genedic[gene_name].long_iso]
            #print(isodic) #@
    else:
        print("No CDS for " + gene_name + ".")
        with open('results.txt', 'a') as f:
            f.write('\n' + gene_name + '\t' + 'No CDS obtained.')

    ##@ save the genedic entery as a pickle file
    if not os.path.exists('genes/' + gene_name + '/' + gene_name + '_save_file.pkl') or fresh:
        with open('genes/' + gene_name + '/' + gene_name + '_save_file.pkl', 'wb') as output:
            pickle.dump (genedic[gene_name], output, pickle.HIGHEST_PROTOCOL)


    return genedic, isodic

def run_commands(iso, isodic, genedic, genomes, tree, speccount, fresh):

    ## try to make the isoform folder
    try:
        os.makedirs('genes/' + isodic[iso].parent_gene + '/' + iso)
    except:
        pass

    ## try to load the isoform from a save
    if os.path.exists('genes/' + isodic[iso].parent_gene + '/' + iso + '/' + iso + '_save_file.pkl') and not fresh:
        try:
            with open('genes/' + isodic[iso].parent_gene + '/' + iso + '/' + iso + '_save_file.pkl', 'rb') as input:
                isodic[iso] = pickle.load(input)

        except:
            print('ERROR: There is no file loaded for ' + iso + '. Try running "--process=isoform_get"') #@ change this error message

    ## try and make the files folder path
    try:
        os.makedirs(isodic[iso].iso_path + '/' + iso + '_files')
    except:
        pass

    ## make a shortcut for the file path
    iso_path = 'genes/' + isodic[iso].parent_gene + '/' + iso + '/'

    ## run the blast search
    if process =='all' or process =='blast':
        isodic[iso].blast_search(blexon (isodic, iso, genomes))

        if debug:
            print(isodic[iso].blast_dic)

    ## align and back-translate with clustal
    if process =='all' or process =='align':
        alngr = 'clus'

        ## align the sequences
        if len(isodic[iso].blast_dic) > 1:
            isodic[iso].load_alignment(alngr, align(isodic, iso, alngr))
            if debug:
                print(isodic[iso].alignment[alngr])

            ## back translate the sequences, prints a file and makes a property of the isoform (second part to come)
            isodic[iso].load_backtrans(alngr, back_translate(isodic[iso].alignment[alngr], isodic[iso].blast_dic, iso, alngr, isodic))

    ## make the tree and run M7-M8 of PAML for clustal
    if process == 'all' or process == 'paml':

        clus_check = False
        mus_check = False
        coff_check = False
        M8a_check = False

        ## make the tree
        treebuild(iso, isodic, genomes, tree, speccount)

        alngr = 'clus'

        ## run PAML
        #speccheck = isodic[iso].spec_number_check(speccount)

        if isodic[iso].spec_number_check(speccount):
            paml (isodic, iso, alngr)
            isodic[iso].loglike_get()

            if isodic[iso].llvals['clus_M7_M8_p'] < .05:
                clus_check = True

            ## if the p-value for clus is below 0.05, run mus
            if clus_check: #@ Do these all need to be nested? How does it iterate through them?
                alngr = 'mus'

                ## align the sequences
                isodic[iso].load_alignment(alngr, align(isodic, iso, alngr))
                if debug:
                    print(isodic[iso].alignment[alngr])

                ## back translate the sequences, prints a file and makes a property of the isoform (second part to come)
                isodic[iso].load_backtrans(alngr, back_translate(isodic[iso].alignment[alngr], isodic[iso].blast_dic, iso, alngr, isodic))

                ## run PAML
                paml (isodic, iso, alngr)
                isodic[iso].loglike_get()
                if isodic[iso].llvals['mus_M7_M8_p'] < 0.05:
                    mus_check = True

            ## if the p-value for clus and mus is below 0.05, run coff
            if clus_check and mus_check:
                alngr = 'coff'

                ## align the sequences
                isodic[iso].load_alignment(alngr, align(isodic, iso, alngr))

                if debug:
                    print(isodic[iso].alignment[alngr])

                ## back translate the sequences, prints a file and makes a property of the isoform (second part to come)
                isodic[iso].load_backtrans(alngr, back_translate(isodic[iso].alignment[alngr], isodic[iso].blast_dic, iso, alngr, isodic))

                ## run PAML
                paml (isodic, iso, alngr)
                isodic[iso].loglike_get()
                if isodic[iso].llvals['coff_M7_M8_p'] < 0.05:
                    coff_check = True

            ## if p-value passes for all of them, get the highest p-value aligner and run M8a
            if clus_check and mus_check and coff_check:
                print ('Picked ' + isodic[iso].llvals['max_p'] + ' to run M8a for ' + iso)
                p8ml (isodic, iso, isodic[iso].llvals['max_p'])
                isodic[iso].loglike_get()
                M8_M8a_p = isodic[iso].llvals['max_p'] + '_M8_M8a_p'
                try:
                    if isodic[iso].llvals[M8_M8a_p] < 0.05:
                        M8a_check = True
                except:
                    print('Add me to the rerun list stupid.') #@ change this you lazy fuck

            ## if M8-M8a passes, then get all the sites under selection
            if clus_check and mus_check and coff_check and M8a_check:
                for aln in ('clus', 'mus', 'coff'):
                    isodic[iso].site_get(aln)
                isodic[iso].site_analysis()

            isodic[iso].print_info()

            ##populate the species_used results file
            genomes_used = []
            for genome in genomes:
                if genome in isodic[iso].species_used:
                    genomes_used.append('1')
                else:
                    genomes_used.append('0')
            with open('species_collected.txt', 'a') as f:
                f.write('\n' + isodic[iso].name + '\t' + '\t'.join(genomes_used))

    ## save and close the isoform object
    with open('genes/' + isodic[iso].parent_gene + '/' + iso + '/' + iso + '_save_file.pkl', 'wb') as output:
        pickle.dump (isodic[iso], output, pickle.HIGHEST_PROTOCOL)

def run_wrapper(gene_name, alliso, fresh):

    #genedic, isodic = isoform_get(gene_name, alliso, fresh) #@
    if genedic: #This only runs if one or more sequences were obtained for the gene
        for iso in isodic:
            run_commands(iso, isodic, genomes, tree, speccount, fresh)



    else:

        print(str(gene_name) + ' cannot be found in the CDS file.')

        with open('results.txt', 'a') as f:
            f.write('\n' + gene_name + '\t' + 'No CDS obtained.')


################################################################
################ classes
################################################################

class isoform(object):
    """Holds information about an isoform"""

    #################
    ## these are used to Initialize and operate on the isoform class
    #################

    def __init__(self, name):
        self.name = name
        self.alignment = {}
        self.backtrans = {}
        self.BEB_sites = {}
        self.BEB_sites['same_sites'] = []
        self.ref_min_aa = {}

    def rename(self, newname):
        ## renames the isoform
        self.name = newname

    def link_parent(self, parent_gene):
        ## used to link an isoform with a parent gene
        self.parent_gene = parent_gene
        self.iso_path = 'genes/' + self.parent_gene + '/' + self.name + '/'

    def ref_nt(self, ref_nt):
        # nt_seq must be a string. also as is, this can't be over-written
        self.ref_nt = ref_nt

    def translate(self):
        # translates the nt_seq if there is one
        try:
            self.ref_aa = trans(self.ref_nt)
        except:
            print("ERROR: There is no nucleotide sequence for this gene")



    #################
    ## these are some property functions to return common file paths
    ## they help make the handling of everything else much more readable
    #################

    def iso_files(self):
        ## returns the path to the files folder
        return self.iso_path + self.name + '_files/'

    def CDS_file(self):
        ## returns the cds_file path
        return self.iso_path + self.name + '_files/' + self.name + '_CDS.txt'

    def prot_file(self):
        ## returns the protein file path
        return self.iso_path + self.name + '_files/' + self.name + '_prot.txt'

    def alignment_file(self, algnr):
        ## returns the file path with the specified alinger
        try:
            return self.iso_path + self.name + '_files/' + self.name + '_' + algnr + '.fasta'
        except:
            print('Need to specify an aligner')

    def paml_file(self, algnr):
        ## returns the paml alignment file path
        try:
            return self.iso_path + self.name + '_files/' + self.name + '_' + algnr + '.paml'
        except:
            print('Need to specify and aligner')

    def tree_file(self):
        ## returns the tree file path
        return self.iso_path + self.name + '_files/' + self.name + '_tree.txt'

    def paml_output(self, algnr):
        ## returns the file path for the paml output file for a given aligner
        try:
            return self.iso_path + self.name + '_files/' + self.name + '_' + algnr + '_PAML_out_full.txt'
        except:
            print('Need to specify and aligner')

    def p8ml_output(self, algnr):
        ## returns the file path for the paml output file for a given aligner
        try:
            return self.iso_path + self.name + '_files/' + self.name + '_' + algnr + '_M8a_PAML_out_full.txt'
        except:
            print('Need to specify and aligner')

    def results_file(self):
        try:
            return self.iso_path + self.name + '_results.txt'
        except:
            print('Could not properly specifiy results file')

    #################
    ## these are used in the processing of the data
    #################

    def blast_search(self, blast_dic):
        # reads in a dicitonary of species:sequence pairs
        self.blast_dic = blast_dic

    def blast_trans(self):
        # translates the blast dicitonary
        try:
            self.blast_prot = {}
            for key in self.blast_dic:
                self.blast_prot[key] = trans(self.blast_dic[key])
        except:
            print("ERROR: There is no blast dictionary")

    def load_alignment(self, aligner, file1):
        ## gives .aln property, which is a dictionary containing aligner names that are in turn dictionaries with alignments
        self.alignment[aligner] = {}
        self.alignment[aligner] = fasta_read(file1)

    def load_backtrans(self, aligner, bctrns):
        ## gives .backtrans property, which is a dictionary containing aligner names that are in turn dictionaries with alignments
        ## bctrns is a dictionary that contaisn the back-translation in nt form
        self.alignment[aligner] = {}
        self.alignment[aligner] = bctrns

    def spec_number_check(self,speccount): ##This makes sure that analysis isn't run if PAML file isn't there
        if len(self.blast_dic) >= speccount:
            return True
        else:
            return False

    def loglike_get(self):
        ## scans the results file for a log-likelihood values
        ## stores properties in a dictionary:
        ## llvals: dictionary of aligner_model : log-likelihood
        self.llvals = {}
        for algr in ('clus','coff','mus'):
            try:
                if os.path.exists(self.paml_output(algr)):
                        with open(self.paml_output(algr), 'r') as f:
                            flag = False
                            for line in f.readlines():

                                if 'Model' in line: ## determines the model
                                    flag = True
                                    line = [x.strip() for x in line]
                                    line = ''.join(line)
                                    line = line.split('Model')[1]
                                    line = line.split(':')[0]
                                    if line: # takes care of some blank lines that get returned
                                        model = algr + '_M' + line

                                if flag and 'lnL' in line: ## gets the log-likelihood value
                                    line = [x.strip() for x in line]
                                    line = ''.join(line)
                                    line = line.split('):')[1]
                                    lnl = line.split('+')[0]
                                    self.llvals[model] = lnl
                                    flag = False
                else:
                    if debug:
                        print('PAML output file for ' + algr + ' not detected.')
            except:
                print('Can not open paml file')

            if 'clus_M7_M8_p' not in self.llvals:
                self.llvals['clus_M7_M8_p'] = 37767



            if os.path.exists(self.p8ml_output(algr)):
                with open(self.p8ml_output(algr), 'r') as f:
                    for line in f.readlines():
                       if 'lnL' in line: ## gets the log-likelihood value
                            model = algr + '_M8a'
                            line = [x.strip() for x in line]
                            line = ''.join(line)
                            line = line.split('):')[1]
                            lnl = line.split('+')[0]
                            self.llvals[model] = lnl
                            flag = False

            ## calculates the 2Delta log-likelihood
            ## calculates a p-value

            M7 = algr + '_M7'
            M8 = algr + '_M8'
            M8a = algr + '_M8a'
            M7_M8_2del = algr + '_M7_M8_2del'
            M7_M8_p = algr + '_M7_M8_p'
            M8_M8a_2del = algr + '_M8_M8a_2del'
            M8_M8a_p = algr + '_M8_M8a_p'

            ## see if there is a M7-M8 comparison
            try:
                self.llvals[M7_M8_2del] = 2 * ( abs( float(self.llvals[M7]) - float(self.llvals[M8]) ))
                self.llvals[M7_M8_p] = chi2.sf ( self.llvals[M7_M8_2del] , 2 )
            except:
                pass #@ And a message here I think.

            ## see if there is a M8-M8a comparison
            try:
                self.llvals[M8_M8a_2del] = 2 * ( abs( float(self.llvals[M8]) - float(self.llvals[M8a]) ))
                self.llvals[M8_M8a_p] = chi2.sf ( self.llvals[M8_M8a_2del] , 1 )
            except:
                pass #@ Message about no significant p-values.

        ## determine the current highest p-value
        d = {}
        for algr in ('clus','mus','coff'):
            try:
                M7_M8_p = algr + '_M7_M8_p'
                d[algr] = self.llvals[M7_M8_p] #@ debug step could be to print these.
            except:
                pass

        # self.llvals['max_p'] = keywithmaxval(d)
        if d:
            self.llvals['max_p'] = max(d, key=d.get)
        else:
            with open('results.txt','a') as f:
                f.write('\n' + self.name + '\tPAML did not run correctly.')

    def site_get(self, algnr):
        # gets the positions of the amino acids under selection in the origonal reference sequence
        # requires: self.llvals['max_p']
        # creates property: self.BEB_sites

        uncor_sites = []  ## this is the list of uncorrected BEB sites
        flag1=False
        flag2=False

        with open(self.paml_output(algnr),'r') as f:
            for line in f.readlines():
                if '+-' in line and not 'Pr(w>1)' in line:
                    line = line.split()
                    line[2] = line[2].replace('*', '')
                    if float(line[2]) > 0.90:
                        uncor_sites.append(line[0])
        if debug:
            print('Sites under selection, before correction (pp > 0.90) for ' + algnr + ': ' + str(uncor_sites))

        ## extracts the D.mel reduced sequence form the PAML input file
        flag1 = False
        count = 0
        red_seq = []
        with open(self.paml_file(algnr),'r') as f:
            for line in f.readlines():
                if clade in line:
                    flag1 = True
                if flag1 and count < 2:
                    line = [x.strip() for x in line]
                    line = ''.join(line)
                    red_seq.append(line)
                    count += 1

        ## two new properties that contain the reduced reference sequences
        self.ref_min_nt = {algnr : red_seq[1]}
        self.ref_min_aa[algnr] = trans(self.ref_min_nt[algnr])

        ## use UPPER CASE to denote the BEB amino acids
        indicies = [] ## since the index position will be sites - 1
        for i in uncor_sites:
            indicies.append(int(i) - 1)
        self.ref_min_aa[algnr] = ("".join(c.upper() if i in indicies else c for i, c in enumerate(self.ref_min_aa[algnr].lower())))

        if debug:
            print(self.ref_min_aa[algnr])
            print(self.ref_aa)
            print('\n')

        ## iterate through the length of the main sequence. Check if the two are
        ## equal. If not, add a '-' to the reduced sequence, and check the next
        ## position. The result is that self.ref_min_aa has dashes and BEB sites
        ## in UPPER CASE, and is the same length as teh ref_aa sequence.
        for i in range(0,len(self.ref_aa)):
            try:
                if debug:
                    print(str(i) + ' - ' + self.ref_min_aa[i] + ' , ' + self.ref_aa[i])
                if self.ref_min_aa[algnr][i].lower() != self.ref_aa[i].lower():
                    if i == 0:
                        self.ref_min_aa[algnr] = '-' + self.ref_min_aa[algnr]
                    else:
                        self.ref_min_aa[algnr] = self.ref_min_aa[algnr][:i] + '-' + self.ref_min_aa[algnr][i:]
            except:
                pass

        # get the list of sites that are in UPPER CASE
        self.BEB_sites[algnr] = []
        for i, c in enumerate(self.ref_min_aa[algnr]):
            if c.isupper():
                self.BEB_sites[algnr].append(i+1)


        if debug:
            print ('Sites under selection (pp > 0.90) for ' + algnr + ': ' + str(self.BEB_sites[algnr]))

    def site_analysis(self):
        ## reduces to the common sites for all aligners, and generates some other data
        ## makes a list called self.BEB_sites['same_sites']
        if len(self.BEB_sites['same_sites']) <= 0: #@ this condition prevents an error if the pickle file exists already
            for i in self.BEB_sites[keywithmaxval(self.BEB_sites)]:
                copy = True
                for aln in ('clus', 'mus', 'coff'):
                    if i not in self.BEB_sites[aln]:
                        copy = False
                if copy:
                    self.BEB_sites['same_sites'].append(i)
            if debug:
                print(self.BEB_sites)

    def print_info(self):

        # print(self.alignment['clus'][clade])
        # print(len(self.alignment['clus'][clade]))
        # print(len(trans(self.alignment['clus'][clade])))

        pamlerror = False
        ## make sure that these variable exists, return NA's if they don't
        if self.llvals['max_p']:

            try:
                M7_M8 = self.llvals[self.llvals['max_p'] + '_M7_M8_p']
            except:
                M7_M8 = 'NA'

            try:
                M8a = self.llvals[self.llvals['max_p'] + '_M8_M8a_p']
            except:
                M8a = 'NA'


            try:
                pa = len(trans(self.alignment['clus'][clade]))/len(self.ref_aa) * 100 #length of the ref_aa of the max pvalue divided by the reference
            except:
                pa = 'NA'

            if self.llvals['clus_M7_M8_p'] == 37767:
                pamlerror = True

        else:
            M7_M8 = 'NA'
            M8a = 'NA'

            try:
                pa = len(trans(self.alignment['clus'][clade]))/len(self.ref_aa) * 100
            except:
                pa = 'NA'

        ## make a variable of the species that were used for analysis
        self.species_used = []
        for key in self.alignment[self.llvals['max_p']]:
            self.species_used.append(key)

        numspec = len(self.species_used)

        ## make sure that these variable exists, return NA's if they don't
        if len(self.BEB_sites['same_sites']) == 0 or self.BEB_sites['same_sites'] == 'NA':
            self.BEB_sites['same_sites_number'] = 0
            self.BEB_sites['same_sites'] = 'NA'
            self.BEB_sites['average_sites'] = 'NA'
        else:
            self.BEB_sites['same_sites_number'] = len(self.BEB_sites['same_sites'])
            self.BEB_sites['average_sites'] = int((len(self.BEB_sites['clus'])
                                                  +len(self.BEB_sites['mus'])
                                                  +len(self.BEB_sites['coff']))/3)


        ## print command for all the output variables
        print('\n' + self.name
              + '\t' + str(M7_M8)
              + '\t' + str(M8a)
              + '\t' + str(self.BEB_sites['same_sites_number'])
              + '\t' + str(self.BEB_sites['same_sites']) + '\n')

        if pamlerror:
            with open('results.txt','a') as f:
                f.write('\n' + self.name + '\tPAML did not run correctly.')
        else:
            with open('results.txt', 'a') as f: #@
                f.write('\n' + self.name
                        + '\t' + str(numspec)
                        + '\t' + str(len(self.ref_aa))
                        + '\t' + str(pa)
                        + '\t' + str(M7_M8)
                        + '\t' + str(M8a)
                        + '\t' + str(self.BEB_sites['average_sites'])
                        + '\t' + str(self.BEB_sites['same_sites_number'])
                        + '\t' + str(self.BEB_sites['same_sites'])
                        + '\t' + str(self.llvals['max_p']))


        with open(self.results_file(), 'w') as f:
            algnr = self.llvals['max_p']

            ##write results to file
            f.write('Gene name:\t' + self.name + '\n')
            f.write('Aligner with the least significant p-value:\t' + algnr + '\n')
            f.write('Number of species analyzed:\t' + str(len(self.species_used)) + '\t(' + ', '.join(self.species_used) + ')\n')
            f.write('M7-M8 pval:\t' + str(M7_M8) + '\n')
            f.write('M8-M8a pval:\t' + str(M8a) + '\n')
            f.write('Number of sites common with all aligners:\t' + str(self.BEB_sites['same_sites_number']) + '\n')
            f.write('Site positions:\t' + str(self.BEB_sites['same_sites']).strip('[]') + '\n\n')


            ## show the two sequences on top of each other
            site_write = True
            try:
                f.write('Clustal aligned sequence:  \t' + self.ref_min_aa['clus'] + '\n')
            except: ## for this one, if it was never made, make it here. This next bit is slightly modified from the site_get function above
                site_write = False

                ## extracts the D.mel reduced sequence form the PAML input file
                flag1 = False
                count = 0
                red_seq = []
                with open(self.paml_file('clus'),'r') as k:
                    for line in k.readlines():
                        if clade in line:
                            flag1 = True
                        if flag1 and count < 2:
                            line = [x.strip() for x in line]
                            line = ''.join(line)
                            red_seq.append(line)
                            count += 1

                ## two new properties that contain the reduced reference sequences
                self.ref_min_nt = {'clus' : red_seq[1]}
                self.ref_min_aa['clus'] = trans(self.ref_min_nt['clus'])

                ## iterate through the length of the main sequence. Check if the two are
                ## equal. If not, add a '-' to the reduced sequence, and check the next
                ## position. The result is that self.ref_min_aa has dashes and BEB sites
                ## in UPPER CASE, and is the same length as teh ref_aa sequence.
                for i in range(0,len(self.ref_aa)):
                    try:
                        if debug:
                            print(str(i) + ' - ' + self.ref_min_aa[i] + ' , ' + self.ref_aa[i])
                        if self.ref_min_aa[algnr][i].lower() != self.ref_aa[i].lower():
                            if i == 0:
                                self.ref_min_aa['clus'] = '-' + self.ref_min_aa['clus']
                            else:
                                self.ref_min_aa['clus'] = self.ref_min_aa['clus'][:i] + '-' + self.ref_min_aa['clus'][i:]
                    except:
                        pass
                self.ref_min_aa['clus'] = self.ref_min_aa['clus'].lower()

                ## finally, write it to the output file
                f.write('Clustal aligned sequence:  \t' + self.ref_min_aa['clus'] + '\n')

            try: ## try for Muscle
                f.write('Muscle aligned sequence:   \t' + self.ref_min_aa['mus'] + '\n')
            except: ## for these two, don't try and make them, just set to false
                site_write = False

            try: ## try for T-Coffee
                f.write('T-Coffee aligned sequence: \t' + self.ref_min_aa['coff'] + '\n')
            except:
                site_write = False

            f.write('Reference sequence:        \t' + self.ref_aa.lower() + '\n')

            ## generate the alignment for sites with a few aa on either side
            if site_write and self.BEB_sites['same_sites_number'] > 0 :
                f.write('\nAlignment of sites under selection:\n- Alignments are centered around the site, with 4 amino acids on either side. \n- There may be more than one site per alignment if they are close together. \n\n')
                for site in self.BEB_sites['same_sites']:
                    lower_bound = int(site - 5)
                    upper_bound = int(site + 4)
                    f.write('Site position:\t' + str(site) + '\n')
                    ref_min_aa_algnr = 'ref_min_aa_' + algnr
                    f.write('PAML sequence:\t' + self.ref_min_aa[algnr][lower_bound : upper_bound] + '\n')
                    f.write('Reference seq:\t' + self.ref_aa.lower()[lower_bound : upper_bound] + '\n\n')


class gene(object):
    """Holds information about a gene"""

    def __init__(self, gene_name):
        self.gene_name = gene_name
        self.all_isoforms = []
        self.uniq_isos = []
        self.scaffolds = {}

    def linkiso(self, string1):
        self.all_isoforms.append(string1)

    def longiso(self, string1):
        self.long_iso = string1

    def uniqiso(self, string1):
        self.uniq_isos.append(string1)

    def add_scaffold(self, species, scaffold_number):
        ## stores the scaffold ID for each species for that gene, to be references for other isoforms
        self.scaffolds[species] = scaffold_number

################################################################
################ run commands
################################################################

## Load the control file
if os.path.exists(ctl):

    ### Make a list of genes in gene file
    with open(gene_file) as f:
        genes = f.read().splitlines()

    with open(ctl) as f:
        settings={}
        for line in f:
            ls = line.split('\t')
            settings[ls[0]]=ls[1].strip('\n')

    if debug:
        print('Settings: ')
        print(settings)

    try: #checks for a clade name
        clade = settings['cladename']
    except:
        pass

    try: #checks for a number of species (this needs to be integrated well into the program)
        speccount = int(settings['species'])
    except:
        pass

    try: #tree must be specified
        tree = settings['tree']
        genomes = settings['tree'].replace('(','').replace(')','').replace(';','').split(',')
        genomes.remove(clade) #@ this line doesn't work if the clade is not the species name. For example, if the clade name is 'Dmel_sub', it can't find it and it isnt removed
    except:
        print('No tree for ' + clade + ' specified.')

    # try: #CDS file must be specified
    if not os.path.exists(settings['CDS']):
        print('ERROR : Cannot find CDS master file')
    cds = fasta_read(settings['CDS'])
    # except:
    #     print('ERROR : No CDS file for ' + clade + ' specified.')
    # if debug:
    #     print(str(len(genes)), "genes")
else:
    print('No control file for ' + clade + '.')

## make sure it knows this is a number
if cores:
    n_cores = int(cores)

## make the output files
with open('results.txt', 'w') as f:
    f.write('Gene name\tSpecies\tReference Length\tPercent aligned\tM7-M8 pval\tM8-M8a pval\tSite average\tNumber of sites\tSite positions\tMax aligner')

with open('species_collected.txt', 'w') as f:
    f.write('Gene Name\t' + '\t'.join(genomes))


## Checks if fresh, and deletes gene folder if so
if fresh and process == 'all':
    for gene_name in genes:
        shutil.rmtree('genes/' + gene_name, ignore_errors=True)
elif fresh and not process == 'all':
    print('The fresh setting only works for the all process setting. Files must be deleted by hand for the other processes.')


## Determines the number of cores to use for parallel
if cores:
    n_cores = int(cores)
else:
    n_cores = multiprocessing.cpu_count()

## Runs Corsair.
genedic = {}
isodic = {}

## uses the gene name to build the genedic and isodic dictionaries
for gene_name in genes:
    genedic, isodic = isoform_get(gene_name, genedic, isodic, alliso, fresh)

## runs the main commands for all isoforms in parallel
Parallel(n_jobs=n_cores)(delayed(run_commands)(iso,isodic,genedic,genomes,tree,speccount,fresh) for iso in isodic.keys())
