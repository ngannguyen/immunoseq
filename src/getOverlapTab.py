#!/usr/bin/env python

'''
nknguyen soe ucsc edu
Feb 14 2012
Aggregate overlap statistics from different experiments into a table
Input: input directory containing overlap-statistic files, one for each experiment
'''

import os, sys, re

class Exp():
    def __init__(self, name):
        self.name = name
        self.cutoffs = []
        self.clones1 = [] #total number of clones in sample 1 passed the cutoffs
        self.clones2 = [] #total number of clones in sample 2 passed the cutoffs
        self.oclones = [] #number of clones (that are passed the cutoffs and are) in both sample 1 and sample 2
        self.reads1o2 = [] #Percentage of sample 1 reads in the overlapped clones
        self.reads2o1 = [] #Percentage of sample 2 reads in the overlapped clones
        self.avrOreads = [] #average of reads1o2 and reads2o1

    def addCutoffStats(self, items):
        self.cutoffs.append( float(items[0]) )
        self.clones1.append( int(items[1]) )
        self.clones2.append( int(items[2]) )
        self.oclones.append( int(items[3]) )
        self.reads1o2.append( float(items[6]) )
        self.reads2o1.append( float(items[7]) )
        self.avrOreads.append( (float(items[6]) + float(items[7]))/2.0 )

class FileFormatError(Exception):
    pass

class InconsistentCutoffsError(Exception):
    pass

def getOverlapReadsTab(exps, file):
    if len(exps) == 0:
        return
    
    cutoffs = exps[0].cutoffs
    for exp in exps:
        for i, c in enumerate(exp.cutoffs):
            if c != cutoffs[i]:
                raise InconsistentCutoffsError("Input files don't have the same cutoffs.")
    
    f = open(file, 'w')
    f.write("Exp\tSamplingSize\t%s\n" %("\t".join(["%.3f" %c for c in cutoffs])) )
    for exp in exps:
        nameitems = exp.name.split('-')
        if len(nameitems) != 2:
            raise ValueError("Wrong filename format. Required: experimentName-samplingsize\n")
        f.write("%s\t%s\t%s\n" %(nameitems[0], nameitems[1], "\t".join([ "%.2f" %p for p in exp.avrOreads ]) ))

    f.close()

def readFile(file):
    name = os.path.basename(file).rstrip(".txt")
    exp = Exp(name)
    f = open(file, 'r')
    for line in f:
        line = line.rstrip('\n')
        if len(line) == 0 or line[0] == '#':
            continue
        items = line.split('\t')
        if len(items) < 8:
            raise FileFormatError( "Wrong file format, file %s. 8 fields were expected\n" %file )
        exp.addCutoffStats(items)
    f.close()
    return exp

def main():
    indir = sys.argv[1]
    if not os.path.isdir(indir):
        raise ValueError("Input directory %s is not a directory\n" %indir)

    exps = []
    for file in os.listdir(indir):
        exp = readFile( os.path.join(indir, file) )
        exps.append(exp)
    orfile = "overlapReads.txt"
    getOverlapReadsTab(exps, orfile)

if __name__ == '__main__':
    main()
