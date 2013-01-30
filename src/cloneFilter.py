#!/usr/bin/env python2.6

#nknguyen soe ucsc edu
#Tue Jun 19 16:53:10 PDT 2012
#Filter low-frequency clones
#Input: input directory of fasta files
#Output: fasta files of clones that pass filtering requirements

import os, sys, re, random, copy
import immunoseq.lib.immunoseqLib as iseqlib
from sonLib.bioio import system
from optparse import OptionParser
from math import log

import matplotlib.pyplot as pyplot
from matplotlib.ticker import *
from matplotlib.font_manager import FontProperties
import numpy as np

class Sample():
    def __init__(self, name):
        self.name = name
        #self.seqs = []
        #self.total = 0

class Seq():
    def __init__(self, seq, count, vs, js, name):
        self.name = name
        self.seq = seq
        self.count = count
        self.freq = -1
        self.vs = vs
        self.js = js
        vstr = ','.join(vs)
        jstr = ','.join(js)
        #self.header =  '|'.join([seq, vstr, jstr]) #header example: CASSLRRGGKPGELFF|TRBV5-4|TRBJ2-2
        self.header =  '|'.join([seq, vstr, jstr, name]) #header example: CASSLRRGGKPGELFF|TRBV5-4|TRBJ2-2

    def setFreq(self, total):
        if total <= 0:
            raise ValueError("Total count <= 0: %d\n" %total)
        self.freq = self.count*100.0/total
    
    def setCount(self, count):
        self.count = count

    def updateCount(self, count):
        self.count += count

    def __cmp__(self, other):
        return cmp(self.seq, other.seq)

def getSeq(header, seq):
    items = header.split(";size=")
    
    if len(items) < 2:
        raise ValueError("Wrong header format: %s\n" %header)
    name = items[0].split(';')[0]

    count = int(items[-1])
    #Find v and j:
    genestr = items[0].split(';')[-1]
    clusters = genestr.split(',,')
    vs = []
    js = []
    for cluster in clusters:
        citems = cluster.split('|')
        currVs = citems[1].split(',')
        for v in currVs:
            if v not in vs:
                vs.append(v)
        
        currJs = citems[2].split(',')
        for j in currJs:
            if j not in js:
                js.append(j)
    return Seq(seq, count, sorted(vs), sorted(js), name)

def addSeq( header, seq, seqs, minCount ):
    seq = getSeq(header, seq)
    if seq.count >= minCount:
        seqs[seq.header] = seq
        #seqs.add( seq )
        #seqs.insert( iseqlib.BinaryTreeNode(seq) )

def readFile(file, minCount):
    #seqs = iseqlib.BinaryTree()
    #seqs = iseqlib.Seqs()
    seqs = {} #key = sequence, val = Seq
    f = open(file, 'r')
    header = ''
    seq = ''
    for line in f:
        line = line.strip()
        if line[0] == '>':
            if header != '' and seq != '' and not re.search("(undefined)", header):
                addSeq( header, seq, seqs, minCount )
            header = line.lstrip('>')
            seq = ''
        else:
            seq = line
    if seq != '' and header != '':
        addSeq( header, seq, seqs, minCount )
    f.close()
    return seqs

def readfiles(indir, minCount, logTransform):
    files = os.listdir( indir )
    samples = []
    for file in files:
        path = os.path.join(indir, file)
        if os.path.isdir(path):
            continue
        print file
        sampleName = file.split('.')[0]
        sample = Sample(sampleName)
        seqs = readFile( path, minCount )
        if logTransform:
            for seq in seqs.values():
                seq.count = log(seq.count, 2)
        #Set frequencies
        total = sum( [s.count for s in seqs.values()] )
        #sample.total = total
        for s in seqs.values():
            s.setFreq(total)
        sample.seqs = seqs
        samples.append(sample)
    return samples

def freqFilter(samples, minPercent, maxPercent):
    newsamples = []
    for sample in samples:
        newsample = Sample(sample.name)
        seqs = {}
        for header, seq in sample.seqs.iteritems():
            if minPercent <= seq.freq and seq.freq <= maxPercent:
                seqs[header] = seq
        newsample.seqs = seqs
        newsamples.append(newsample)

    return newsamples

################## SAMPLING ###################
#Different samples have different level of sequencing, and hence will create bias. For example, large sample (lots of sequences got amplified and sequenced) will inherently have more overlapping with other samples than a smaller sample. THe purpose of sampling is to make sure that every sample has the same number of starting sequences to avoid the bias
def samplingSample_uniq(sample, size):
    newseqs = {}
    sys.stderr.write("Sampling sample: %s with %d uniq sequences. Sampling size: %d\n" %(sample.name, len(sample.seqs.keys()), size))
    chosenHeaders = random.sample( sample.seqs.keys(), size)
    for header in chosenHeaders:
        newseqs[header] = copy.copy( sample.seqs[header] )
    total = sum( [s.count for s in newseqs.values()] )
    for s in newseqs.values():
        newcount = int( s.freq*total/100.0 )
        s.setCount(newcount)
        #s.setFreq(total)
    sample.seqs = newseqs
    return

def samplingSample_weightedUniq(sample, size):
    newseqs = {}
    freq2seqs = {}
    seqs = sorted( sample.seqs.values(), reverse=True, key=lambda seq:seq.freq )
    print sample.name
    print size
    print len(seqs)

    for seq in seqs:
        if seq.freq not in freq2seqs:
            freq2seqs[ seq.freq ] = [seq]
        else:
            freq2seqs[ seq.freq ].append(seq)
    
    freqs = [ seq.freq for seq in seqs ]
    i = 0
    numseqToAdd = size
    while numseqToAdd > 0 and i < len(freqs):
        freq = freqs[i]

        freqseqs = freq2seqs[freq]
        numseq = len( freqseqs )

        print "Freq: %f, numseqToAdd: %d, numseq: %d\n" %(freq, numseqToAdd, numseq)

        selectedSeqs = freqseqs
        if numseqToAdd < numseq: #number of sequences to add is less than number of sequences with current frequency: randomly select
            selectedSeqs = random.sample( freqseqs, numseqToAdd )
            
        for s in selectedSeqs:
            newseqs[ s.header ] = copy.copy( s )
        numseqToAdd -= len(selectedSeqs)
        i += 1

    #Normalize the count:
    total = sum( [ seq.count for seq in newseqs.values() ] )
    for s in newseqs.values():
        newcount = int( s.freq*total/100.0 )
        s.setCount(newcount)
    sample.seqs = newseqs
    print len(newseqs)
    return

#def samplingSample_weightedUniq(sample, size, numDraws):
#    #NOT WORKING!
#    sys.stderr.write("Sampling sample: %s with %d uniq sequences. Sampling size: %d\n" %(sample.name, len(sample.seqs.keys()), size))
#    newseqs = {}
#    seqs = sorted( sample.seqs.values(), reverse=True, key=lambda seq: seq.freq )
#    if len(seqs) < size:
#        return
#    freqs = [ seq.freq/100.0 for seq in seqs ]
#    if sum(freqs) > 1:
#        sys.stderr.write("Sample %s, sum of freqs is %f > 1\n" %(sample.name, sum(freqs) ))
#
#   #numDraws = max( [sum([ seq.count for seq in seqs ]), len(seqs)*2] )
#    counts = np.random.multinomial( numDraws, freqs )
#    
#    numseq = 0
#    #Look for non-zero cells:
#    indices = []
#    for i, c in enumerate(counts):
#        if c > 0:
#            indices.append(i)
#    #Select 'size' cells from indices:
#    selectedIndices = random.samle( indices, size )
#
#    for i in indices:
#        c = counts[i]
#        seq = Seq( seqs[i].seq, c, seqs[i].vs, seqs[i].js )
#        newseqs[seq.header] = seq
#        numseq += 1
#
#    total = sum( [s.count for s in newseqs.values()] )
#    for s in newseqs.values():
#        s.setFreq(total)
#    sample.seqs = newseqs
#    return

def samplingSample(sample, size):
    sys.stderr.write("Begin to sample %d sequences from sample %s\n" %(size, sample.name))

    seqs = sample.seqs
    #newseqs = iseqlib.Seqs()
    newseqs = {}
    indexList = []#index list represents all the sequences of the sample
    for i, s in enumerate(seqs.values()):
        indexList.extend( [i for j in xrange(s.count)] )
    sys.stderr.write("\tDone with the indexList...\n")
    sampleTotalReads = sum( [s.count for s in seqs.values()] )
    if sampleTotalReads < size:
        return

    headers = seqs.keys()
    chosenIndices = random.sample(indexList, size)
    sys.stderr.write("\tDone sampling from the indexList, size %d...\n" %len(chosenIndices))

    for i in chosenIndices:
        header = headers[i]
        if header in newseqs:
            newseqs[header].updateCount(1)
        else:
            newseqs[header] = copy.copy( seqs[header] )
            newseqs[header].setCount(1)
    sys.stderr.write("\tDone filling in newseqs...\n")

    #Update sample with new sequences:
    #Set frequencies
    total = sum( [s.count for s in newseqs.values()] )
    for s in newseqs.values():
        s.setFreq(total)
    sample.seqs = newseqs
    return

def sampling(samples, size, uniq, weightedUniq):
    '''Randomly pick 'size' sequences from each sample'''
    #maxTotal = max( [ sum([seq.count for seq in sample.seqs.values()])  for sample in samples ] )

    for sample in samples:
        if weightedUniq:
            samplingSample_weightedUniq(sample, size)
        elif uniq:
            samplingSample_uniq(sample, size)
        else:
            samplingSample(sample, size)
        sys.stderr.write("Done sampling sample %s\n" %sample.name)
    return
################# END SAMPLING ################

################ PRINT OUT SAMPLED FASTAS ###############
def printFasta(samples, outdir):
    for sample in samples:
        outfile = os.path.join(outdir, "%s.fa" %sample.name)
        f = open(outfile, 'w')
        for header, seq in sample.seqs.iteritems():
            #f.write(">%s;freq=%f|%s|%s;size=%d\n" %(sample.name, seq.freq, ",".join(seq.vs), ",".join(seq.js), seq.count))
            f.write(">%s;freq=%f|%s|%s;size=%d\n" %(seq.name, seq.freq, ",".join(seq.vs), ",".join(seq.js), seq.count))
            f.write("%s\n" %seq.seq)
        f.close()

def checkOptions(parser, args, options):
    if options.indir == '':
        parser.error("Input directory was not specified\n")
    if options.outdir == '':
        parser.error("Output directory was not specified\n")

def addOptions(parser):
    parser.add_option('-i', '--indir', dest='indir', help="Required argument. Input directory")
    parser.add_option('-o', '--outdir', dest='outdir', help="Required argument. Output directory")
    parser.add_option('-c', '--minCount', dest='minCount', default=1, type='int', help='Minimun clone count to be included. Default=%default')
    parser.add_option('-p', '--minPercent', dest='minPercent', default=0.0, type='float', help='Minimum clone frequency to be included. Default=%default')
    parser.add_option('-m', '--maxPercent', dest='maxPercent', default=100.0, type='float', help='Maximum clone frequency to be included. Default=%default')
    parser.add_option('-l', '--logTransform', dest='logTransform', action='store_true', default=False, help='If specified, will perform a log 2 transform on the read counts of each clones')

def main():
    usage=('Usage: %prog [options]')
    parser = OptionParser(usage=usage)
    addOptions(parser)
    options, args = parser.parse_args()
    checkOptions(parser, args, options)

    samples = readfiles(options.indir, options.minCount, options.logTransform)
    if options.minPercent > 0.0 or options.maxPercent < 100.0:
        samples = freqFilter(samples, options.minPercent, options.maxPercent)
    
    printFasta(samples, options.outdir)

if __name__ == '__main__':
    main()

