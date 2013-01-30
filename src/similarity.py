#!/usr/bin/env python

'''
nknguyen
Jun 13 2012
Calculates similarity index
Input: input directory containing fasta files
Output: Morisita - Horn similarity index
'''

import os, re, sys
from math import log
import immunoseq.lib.immunoseqLib as iseqlib

def morisitaHorn(data1, data2):
    X = sum(data1)
    Y = sum(data2)

    if len(data1) != len(data2):
        raise ValueError("Error while calculating MorisitaHorn similarity index. The two input data lists must have the same length!\n")
    
    if X == 0 or Y == 0:
        return -1

    sumXiYi = 0
    sumXiSqr = 0
    sumYiSqr = 0

    for i, x in enumerate(data1):
        y = data2[i]
        sumXiYi += x*y
        sumXiSqr += x*x
        sumYiSqr += y*y

    numerator = 2*sumXiYi
    denominator = ( float(sumXiSqr)/(X*X) + float(sumYiSqr)/(Y*Y) )*X*Y

    return float(numerator)/denominator

def log10Transform(data):
    logdata = []
    for d in data:
        if d > 0:
            d = log(d, 10)
        logdata.append(d)
    return logdata        

def pairwiseSimilarity(sample1, sample2, aa2v2j1, aa2v2j2):
    seqs1 = sample1.seqs
    seqs2 = sample2.seqs
    
    data1 = []
    data2 = []

    for header1, seq1 in seqs1.iteritems():
        count1 = seq1.count
        data1.append( count1 )

        #check to see if this sequence is in sample 2 
        #seq2 = iseqlib.findSeq( seq1, aa2v2j2 )
        oseqs2 = iseqlib.findSeqs( seq1, aa2v2j2 )

        count2 = 0
        if len(oseqs2) > 0:
            count2 = sum([seq2.count for seq2 in oseqs2])
        data2.append( count2 )

    #Fill in data2 with sequences that are uniquely of sample 2
    for header2, seq2 in seqs2.iteritems():
        count2 = seq2.count
        #seq1 = iseqlib.hasSeq( seq2, aa2v2j1 )
        if not iseqlib.hasSeq(seq2, aa2v2j1):
            data2.append( count2 )
            data1.append( 0 )

    if sample1.name == sample2.name:
        print sample1.name
        print len(sample1.seqs)
        print len(data1)
        print len(data2)
        print data1 == data2
        for i,d in enumerate(data1):
            d2 = data2[i]
            if d != d2:
                print "Different values: %d, %d" %(d, d2)


    #data1 = log10Transform(data1)
    #data2 = log10Transform(data2)
    similarity = morisitaHorn(data1, data2)
    return similarity

def addOptions( parser ):
    parser.add_option('-i', '--indir', dest='indir', help='Input directory')
    parser.add_option('-o', '--outdir', dest='outdir', help='Output directory')
    parser.add_option('-c', '--mincount', dest='mincount', type='int', default=1, help='Minimum read count. Default=%default')
    #parser.add_option('-f', '--freq', dest='freq', type='float', default=0.0, help='Minimum frequency. Default=%default %')
    #parser.add_option('-s', '--minsam', dest='minsam', type='int', default=1, help='Minimum number of samples. Default=%default')

def main():
    parser = iseqlib.initOptions()
    addOptions(parser)

    options, args = parser.parse_args()
    
    mode = 1
    samples = iseqlib.readfiles(options.indir, options.mincount, mode)
    sample2aa2v2j = iseqlib.getsample2aa2v2j(samples)
    
    #### MORISITA-HORN SIMILARITY INDEX ####
    outfile = os.path.join(options.outdir, "similarity.txt")
    f = open(outfile, 'w')
    #f.write("Sample1\tSample2\tMorisitaHorn Similarity Index\n")
    names = sorted([s.name for s in samples])
    name2sample = {}
    for s in samples:
        name2sample[ s.name ] = s

    f.write("\t%s\n" % '\t'.join(names) )
    #for i in xrange( len(names) -1 ):
    for i in xrange( len(names) ):
        s1 = name2sample[ names[i] ]
        aa2v2j1 = sample2aa2v2j[s1.name]
        f.write("%s" %names[i])
        #f.write("%s" % '\t'.join( ['' for t in xrange(i + 2)] ))
        #for j in xrange( i+1, len(names) ):
        for j in xrange( len(names) ):
            s2 = name2sample[ names[j] ]
            aa2v2j2 = sample2aa2v2j[s2.name]
            similarity = pairwiseSimilarity(s1, s2, aa2v2j1, aa2v2j2)
            #f.write("%s\t%s\t%f\n" %(s1.name, s2.name, similarity))
            f.write("\t%f" %similarity)
        f.write("\n")
    f.close()

if __name__ == '__main__':
    main()
