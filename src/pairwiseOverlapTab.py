#!/usr/bin/env python

'''
Jun 04 2012
Latex Table of pairwise overlap statistics of two samples.
Input: output tab-separated outfile of overalp.py stats
Output: latex table
'''

import os, re, sys
import immunoseq.lib.immunoseqLib as iseqlib

def tabHeader(f):
    f.write('\\begin{sidewaystable}\n')
    f.write('\\centering\n')
    f.write("\\scalebox{1}{%\n")
    f.write("\\begin{tabular}{r|r|r|r|r|r|r|r|r}\n")
    f.write("\\hline\n")
    f.write("\%Cutoff & Clones1 & Overlap1 & \%OClones1 & \%OReads1 & Clones2 & Overlap2 & \%OClones2 & \%OReads2\\\\\n")
    f.write("\\hline\n")

def tab(f, cutoff2line, field2index):
    cutoffs = sorted(cutoff2line.keys())
    fields = ['%Cutoff', 'Clones1', 'Overlap1', '%1overlap2', '%reads1overlap2', 'Clones2', 'Overlap2', '%2overlap1', '%reads2overlap1']

    for cutoff in cutoffs:
        line = cutoff2line[cutoff]
        items = line.split("\t")
        for i, field in enumerate(fields):
            index = field2index[field]
            item = items[index]
            f.write("%s" %item)
            if i < len(fields) -1:
                f.write(" & ")
        f.write("\\\\\n")
        f.write("\\hline\n")


#---- main -----

field2index = {'%Cutoff':0, 'Clones1':1, 'Clones2':2, 'Overlap1':3, 'Overlap2':4, '%1overlap2':5, '%2overlap1':6, '%reads1overlap2':7, '%reads2overlap1':8}
cutoff2line = {}

for line in sys.stdin:
    line = line.strip('\n')
    if len(line) == 0:
        continue
    if line[0] == '#':
        line = line.lstrip('#')
        items = line.split('\t')
        for i, item in enumerate(items):
            field2index[item] = i
    else:
        items = line.split('\t')
        cutoffIndex = field2index['%Cutoff']
        cutoff = float(items[cutoffIndex])
        cutoff2line[cutoff] = line


iseqlib.writeDocumentStart(sys.stdout)
tabHeader(sys.stdout)
tab(sys.stdout, cutoff2line, field2index)
label=''
captionStr = ''
iseqlib.sidewaystableCloser(sys.stdout, captionStr, label)
iseqlib.writeDocumentEnd(sys.stdout)


