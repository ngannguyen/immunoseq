#!/usr/bin/env python2.6

"""
nknguyen at soe ucsc edu
Nov 26 2012
"""

import os, sys, re
from Bio.Blast import NCBIXML

#Read in blast-output file:
rh = open(sys.argv[1])
records = NCBIXML.parse( rh)

clone2hits = {} #key = cloneName, val = [ (list of hit papers, identity) ]
clones = []

minPositives=float(sys.argv[3])
for record in records:
    clone = record.query
    querylen = record.query_length

    clones.append(clone)
    for aln in record.alignments:
        for hit in aln.hsps:
            #if hit.positives < minPositives:
            if float(hit.positives)/len(hit.query) < minPositives:
                continue
            if len(hit.query) < querylen:
                continue
            if clone in clone2hits:
                clone2hits[ clone ].append( ( aln.title, hit.identities, hit.query, hit.match, hit.sbjct ) )
            else:
                clone2hits[ clone ] = [ (aln.title, hit.identities, hit.query, hit.match, hit.sbjct) ]

outfile = sys.argv[2]
allfh = open(outfile, 'w')
allfh.write("#MinPos %f\n" %minPositives)
cutoff = 1

#Print out clones that have match(es) with at least minPositives
for clone in clones:
    items = clone.split(';')
    id = '.'.join( [items[0], items[1]] )
    
    if clone in clone2hits:
        hits = clone2hits[ clone ]
        hitstr = '\t'.join([ hit[0] for hit in hits ])
        #allfh.write( '\n%s\t%s\n' %(clone, hitstr) )
        allfh.write("\n>%s\n" %clone)
        for i, hit in enumerate(hits):
            allfh.write("\t%d/ %s\n" %(i, hit[0]))
            allfh.write("\t\t%s\n" %hit[2])
            allfh.write("\t\t%s\n" %hit[3])
            allfh.write("\t\t%s\n" %hit[4])
        
def getPc(count, total):
    if total == 0:
        return 0
    return 100.0*count/total

total = len(clones)
numhits = len(clone2hits)
allfh.write("\n### Summary ###\n")
allfh.write("Total\tNumHits\t% hits/total\n")
allfh.write("%d\t%d\t%f\n" %(total, numhits, getPc(numhits, total)))
allfh.close()




