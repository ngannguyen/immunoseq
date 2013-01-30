#!/usr/bin/env python2.6

'''
05/01/2012 nknguyen soe ucsc edu
Input: input directory of fasta files
Output: Plot showing sequence length distribution
'''

import os, sys, re
import matplotlib.pyplot as pyplot
from immunoseq.lib.immunoseqLib import *
from immunoseq.lib.lendistLib import *

def getLenDist(seqs):
    len2freq = {} #key = len, values = (readFreq, uniqFreq)
    for s in seqs.values():
        l = len(s.seq)
        if l not in len2freq:
            len2freq[l] = (s.count, 1)
        else:
            len2freq[l][0] += s.count
            len2freq[l][1] += 1
    
    total = len(seqs)
    totalreads = sum([s.count for s in seqs.values()])
    if total > 0:
        for l, c in len2freq.iteritems():
            len2freq[l][0] = c[0]*100.0/totalreads
            len2freq[l][1] = c[1]*100.0/total
    return len2freq

def getSampleLenDist(samples):
    sample2len2freq = {}
    lens = []
    for sample in samples:
        len2freq = getLenDist( sample.seqs )
        sample2len2freq[sample.name] = len2freq
        for l in len2freq.keys():
            if l not in lens:
                lens.append(l)

    for s, l2f in sample2len2freq.iteritems():
        for l in lens:
            if l not in l2f:
                l2f[l] = (0.0, 0.0)
    return sample2len2freq, sorted(lens)

def main():
    parser = initOptions3()
    initPlotOptions( parser )
    
    options, args = parser.parse_args()
    checkPlotOptions( options, parser )
    checkOptions3(parser, args, options)

    mincount = 1
    samples = readfiles(options.indir, mincount, 1)
    
    uniq = True
    outfile = os.path.join(options.outdir, "lenDist")
    drawLenDist(samples, outfile, uniq)
    outfileRead = os.path.join(options.outdir, "lenDist-read")
    drawLenDist(samples, outfileRead, not uniq)


if __name__ == '__main__':
    main()

