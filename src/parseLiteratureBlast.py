#!/usr/bin/env python2.6

"""
nknguyen at soe ucsc edu
Sep 1
"""

import os, sys, re
from Bio.Blast import NCBIXML

def reformatTitle(title):
    title = title.replace("_", " ")
    title = title.replace(";", " ")
    return title

def checkKeywords(title, keywords):
    t = title.lower()
    for keyword in keywords:
        if re.search(keyword, t):
            return True
    return False

def checkNumSamples(title):
    items = title.split(';')
    samples = items[0].lstrip('>').split(',')
    sample2patient={'SBC8': 'B', 'asBD':'B', 'asBR':'B', 
                    'SBC7': '1', 'as1D':'1', 'as1R':'1',
                    'adapt11D': '11', 'as11R': '11'}
    patients = []
    for s in samples:
        if s not in sample2patient:
            patient = s
        else:
            patient = sample2patient[s]
        if patient not in patients:
            patients.append(patient)
    return len(patients)

#Read in blast-output file:
#rh = open("all-sorted.xml")
rh = open(sys.argv[1])
records = NCBIXML.parse( rh)

clone2hits = {} #key = cloneName, val = [ (list of hit papers, identity) ]
clones = []

#minPositives=int(sys.argv[3])
minPositives=float(sys.argv[3])
for record in records:
    clone = record.query
    #HACK
    #numsams = checkNumSamples(clone)
    #if numsams < 4:
    #    continue

    clones.append(clone)
    for aln in record.alignments:
        for hit in aln.hsps:
            #if hit.positives < minPositives:
            if float(hit.positives)/len(hit.query) < minPositives:
                continue
            if clone in clone2hits:
                clone2hits[ clone ].append( (reformatTitle(aln.title), hit.identities, hit.query, hit.match) )
            else:
                clone2hits[ clone ] = [ (reformatTitle(aln.title), hit.identities, hit.query, hit.match) ]

outfile = sys.argv[2]
allfh = open(outfile, 'w')
allfh.write("#MinPos %f\n" %minPositives)
cutoff = 1

autoimmuneKeywords=['arthritis', 'ankylosing', 'spondy', 'autoreactive', 'autoantigen', 'reactive arthritis', 'rheumatoid arthritis', 'multiple sclerosis', 'self', 'cross-reactive', 'mimicry', 'synovial', 'crohn', 'psoriasis', 'inflammatory bowel disease', 'ibd', 'ulcerative colitis', 'uveitis']
b27Keywords=['b27']
pathogenKeywords=['chlamydia', 'salmonella', 'yersinia', 'shigella', 'campylobacter', 'vipr1', 'ebv', 'epstein-barr', 'lmp2']

numAuto = 0 #number of clones with at least one hit passed cutoff and matches one of the autoimmuneKeywords
numB27 = 0
numPathogen = 0

#Print out clones that have match(es) with at least minPositives
for clone in clones:
    items = clone.split(';')
    id = '.'.join( [items[0], items[1]] )
    
    if clone in clone2hits:
        hits = clone2hits[ clone ]
        hitstr = '\t'.join([ hit[0] for hit in hits ])
        #allfh.write( '\n%s\t%s\n' %(clone, hitstr) )
        allfh.write("\n>%s\n" %clone)
        matchAuto = False
        matchB27 = False
        matchPathogen = False
        for i, hit in enumerate(hits):
            allfh.write("\t%d/ %s\n" %(i, hit[0]))
            allfh.write("\t\t%s\n" %hit[2])
            allfh.write("\t\t%s\n" %hit[3])
        
            #Check to see if any keywords matched:
            if not matchAuto and checkKeywords(hit[0], autoimmuneKeywords):
                matchAuto = True
            if not matchB27 and checkKeywords(hit[0], b27Keywords):
                matchB27 = True
            if not matchPathogen and checkKeywords(hit[0], pathogenKeywords):
                matchPathogen = True
        
        if matchAuto:
            numAuto += 1
        elif matchB27:
            numB27 += 1
        elif matchPathogen:
            numPathogen += 1

def getPc(count, total):
    if total == 0:
        return 0
    return 100.0*count/total

total = len(clones)
numhits = len(clone2hits)
allfh.write("\n### Summary ###\n")
allfh.write("Total\tNumHits\t% hits/total\tnumAuto\t% auto/total\t% auto/hits\tnumB27\t% b27/total\t% b27/hits\tnumPathogen\t% pathogen/total\t% pathogen/hits\n")
allfh.write("%d\t%d\t%f\t%d\t%f\t%f\t%d\t%f\t%f\t%d\t%f\t%f\n" %(total, numhits, getPc(numhits, total), numAuto, getPc(numAuto, total), getPc(numAuto, numhits), numB27, getPc(numB27, total), getPc(numB27, numhits), numPathogen, getPc(numPathogen, total), getPc(numPathogen, numhits)) )
allfh.close()




