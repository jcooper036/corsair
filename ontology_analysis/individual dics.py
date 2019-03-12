#! /usr/bin/env python2

import pprint
import pickle 

kegg_db = open("hsaKEGG.txt", "r")
CLASS_A = {}
CLASS_B = {}
CLASS_C = {}

def class_A(db):
    with open("hsaKEGG.txt",'r') as db:
            for line in db:
                x = line.split(',')
                if x[0] == "A":
                    A = str(x[1:])
                    genelist = []
                    CLASS_A.update( { A : genelist } )
                if x[0] == "D":
                    B = str(x[2]).upper()
                    CLASS_A.update( A = genelist.append(B) )
    CLASS_A.pop("A",None)                
    pp = pprint.PrettyPrinter(stream=open('testdicA.txt','w'))
    pp.pprint(CLASS_A)


class_A(kegg_db)


def class_B(db):
    with open("hsaKEGG.txt",'r') as db:
            for line in db:
                x = line.split(',')
                if x[0] == "B":
                    AB = str(x[1:])
                    genelist = []
                    CLASS_B.update( { AB : genelist } )
                elif x[0] == "D":
                    B = str(x[2]).upper()
                    CLASS_B.update( AB = genelist.append(B) )
    CLASS_B.pop("AB",None)
    pp = pprint.PrettyPrinter(stream=open('testdicB.txt','w'))
    pp.pprint(CLASS_B)

class_B(kegg_db)


def class_C(db):
    with open("hsaKEGG.txt",'r') as db:
            for line in db:
                x = line.split(',')
                if x[0] == "C":
                    A = str(x[1:])
                    genelist = []
                    CLASS_C.update( { A : genelist } )
                if x[0] == "D":
                    B = str(x[2]).upper()
                    CLASS_C.update( A = genelist.append(B) )
    CLASS_C.pop("A",None)                 
    pp = pprint.PrettyPrinter(stream=open('testdicC.txt','w'))
    pp.pprint(CLASS_C)

class_C(kegg_db)

CLASS_B_COUNTS = {}

bats = open("BatHits.txt" , "r")

def counting_B(input):
    inputlist = []
    countdic = {}
    notfound = []
    hits = []
    with open ("BatHits.txt", "r") as input:
        for GENE_ID in input:
            inputlist.append(str(GENE_ID).rstrip())
        for category in CLASS_B:
            for gene in CLASS_B[category]:
                if gene in inputlist:
                   hits.append(gene)
            countdic.update( { category : hits } )
    for gene in inputlist:
        if gene not in hits:
            notfound.append(gene)
    print notfound
        
    #for category in countdic:
        #print category + str(len(countdic[category]))
                                  
    pp = pprint.PrettyPrinter(stream=open('countdicB.txt','w'))
    pp.pprint(countdic)

counting_B(bats)


                
                
            
    


    
