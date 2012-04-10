#!/usr/bin/env python2.6

#nknguyen soe ucsc edu
#March 20 2012
#Input: input directory of fasta files
#Output: fasta files of a sampled subset of the original files with a specific size

import os, sys, re, random, copy
import immunoseq.lib.immunoseqLib as iseqlib
from sonLib.bioio import system
from optparse import OptionParser

import matplotlib.pyplot as pyplot
from matplotlib.ticker import *
from matplotlib.font_manager import FontProperties

class Sample():
    def __init__(self, name):
        self.name = name
        #self.seqs = []
        #self.total = 0

class Seq():
    def __init__(self, seq, count, vs, js):
        self.seq = seq
        self.count = count
        self.freq = -1
        self.vs = vs
        self.js = js
        vstr = ','.join(vs)
        jstr = ','.join(js)
        self.header =  '|'.join([seq, vstr, jstr]) #header example: CASSLRRGGKPGELFF|TRBV5-4|TRBJ2-2

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
    return Seq(seq, count, sorted(vs), sorted(js))

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
            if header != '' and seq != '':
                addSeq( header, seq, seqs, minCount )
            header = line.lstrip('>')
            seq = ''
        else:
            seq = line
    if seq != '' and header != '':
        addSeq( header, seq, seqs, minCount )
    f.close()
    return seqs

def readfiles(indir, minCount):
    files = os.listdir( indir )
    samples = []
    for file in files:
        path = os.path.join(indir, file)
        sampleName = file.split('.')[0]
        sample = Sample(sampleName)
        seqs = readFile( path, minCount )

        #Set frequencies
        total = sum( [s.count for s in seqs.values()] )
        #sample.total = total
        for s in seqs.values():
            s.setFreq(total)
        sample.seqs = seqs
        samples.append(sample)
    return samples

################## SAMPLING ###################
#Different samples have different level of sequencing, and hence will create bias. For example, large sample (lots of sequences got amplified and sequenced) will inherently has more overlapping with other samples than a smaller sample. THe purpose of sampling is to make sure that every sample has the same number of starting sequences to avoid the bias
def samplingSample_uniq(sample, size):
    newseqs = {}
    chosenHeaders = random.sample( sample.seqs.keys(), size)
    for header in chosenHeaders:
        newseqs[header] = copy.copy( sample.seqs[header] )
    total = sum( [s.count for s in newseqs.values()] )
    for s in newseqs.values():
        s.setFreq(total)
    sample.seqs = newseqs
    return

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
    sys.stderr.write("\tDone sampling from the indexList...\n")

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

def sampling(samples, size, uniq):
    '''Randomly pick 'size' sequences from each sample'''
    for sample in samples:
        if uniq:
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
            f.write(">%s;freq=%f|%s|%s;size=%d\n" %(sample.name, seq.freq, ",".join(seq.vs), ",".join(seq.js), seq.count))
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
    parser.add_option('-s', '--sampling', dest='sampling', default=-1, type='int', help='Normalized size wished to look at. Default=-1, meaning no normalization')
    parser.add_option('-u', '--samplingUniq', dest='samplingUniq', action='store_true', default=False, help='If specified, sampling on the set of unique sequences instead of all sequences')

def main():
    usage=('Usage: %prog [options]')
    parser = OptionParser(usage=usage)
    addOptions(parser)
    options, args = parser.parse_args()
    checkOptions(parser, args, options)

    samples = readfiles(options.indir, options.minCount) 
    if options.sampling > 0:
        sampling(samples, options.sampling, options.samplingUniq)
    printFasta(samples, options.outdir)

if __name__ == '__main__':
    main()

