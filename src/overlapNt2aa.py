#!/usr/bin/env python2.6

'''
05/09/2012 nknguyen soe ucsc edu
Comparing number of uniq nucleotide sequences each 
Input: input directory of fasta files
Output: 
'''

import os, sys, re, copy, gzip
import matplotlib.pyplot as pyplot
import immunoseq.lib.immunoseqLib as iseqlib

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
import cPickle as pickle

from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import getTempDirectory
from sonLib.bioio import setLogLevel
import numpy as np


class Sample():
    def __init__(self, name):
        self.name = name
        #self.seqs = []
        #self.total = 0

class Seq():
    def __init__(self, name, seq, count, vs, js, freq):
        self.name = name
        self.seq = seq
        self.count = count
        self.freq = freq
        self.vs = sorted(vs)
        self.js = sorted(js)
        vstr = ','.join(vs)
        jstr = ','.join(js)
        self.aa = iseqlib.nt2aa(self.seq)
        #self.header =  '|'.join([seq, vstr, jstr]) #header example: CASSLRRGGKPGELFF|TRBV5-4|TRBJ2-2
        self.header =  '|'.join([self.aa, vstr, jstr]) #header example: CASSLRRGGKPGELFF|TRBV5-4|TRBJ2-2
    
    def getFastaHeader(self):
        return ">%s;freq=%f|%s|%s;size=%d" %(self.name, self.freq, ','.join(self.vs) , ','.join(self.js), self.count )

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
    name = items[0].split(';')[0]
    genestr = items[0].split(';')[-1]

    freq = -1.0
    if re.search("freq=", genestr):
        freq = float( genestr.split('|')[0].lstrip('freq=') )

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
    
    seq =  Seq(name, seq, count, sorted(vs), sorted(js), freq)
    return seq

def addSeq( header, seq, seqs, minCount ):
    seq = getSeq(header, seq)
    if seq.count >= minCount:
        if seq.header not in seqs:
            seqs[seq.header] = {seq.seq: seq}
        else:
            seqs[seq.header][seq.seq] = seq

def readFile(file, minCount):
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
        if os.path.isdir(path) or file.split('.')[-1] != 'fa':
            continue
        sampleName = file.split('.')[0]
        sample = Sample(sampleName)
        seqs = readFile( path, minCount )

        #Set frequencies
        #total = sum( [s.count for s in seqs.values()] )
        #sample.total = total
        #for s in seqs.values():
        #    s.setFreq(total)
        sample.seqs = seqs
        samples.append(sample)
    return samples

############### SET UP JOBTREE #####################
class Setup(Target):
    """Setting up the samplings
    """
    def __init__(self, options):
        Target.__init__(self)
        self.options = options

    def run(self):
        setLogLevel("DEBUG")
        globalTempDir = self.getGlobalTempDir() 
        for i in xrange( self.options.samnum ) :
            samplingdir = os.path.join(globalTempDir, "%d" %i)
            system("mkdir -p %s" %samplingdir)
            self.addChildTarget( RunSampling(samplingdir, self.options) )
        
        self.setFollowOnTarget( RunDrawNt2aaDist(globalTempDir, self.options) )
            
class RunSampling(Target):
    """Sampling 'samsize' sequences from each sample, and compute the nt2aa dist from the sampled seqs
    """
    def __init__(self, outdir, options):
        Target.__init__(self)
        self.outdir = outdir
        self.options = options

    def run(self):
        #Sampling from the samples:
        globalTempDir = self.getGlobalTempDir()
        if not self.options.samsize:
            fadir = self.options.indir
        else:
            fadir = os.path.join(globalTempDir, 'fa')
            system("mkdir -p %s" %fadir)
            #system("sampling.py -i %s -o %s -u -s %d -c %d" %(self.options.indir, fadir, self.options.samsize, self.options.mincount) )
            system("sampling.py -i %s -o %s -s %d -c %d" %(self.options.indir, fadir, self.options.samsize, self.options.mincount) )
        
        #DEBUG:
        #samplingNum = self.outdir.split('/')[-1]
        #system("cp -R %s %s" %(fadir, os.path.join(self.options.outdir, samplingNum) ))
        #logger.info("Done sampling")
        self.setFollowOnTarget( RunGetSamplingData(fadir, self.outdir) ) 

class RunGetSamplingData(Target):
    def __init__(self, indir, outdir):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir

    def run(self):
        mincount = 1
        samples = readfiles( self.indir, 1 )
        sam2nt2aa = getNumNtPerAa(samples)
        pickleFile = os.path.join(self.outdir, "sam2nt2aa.pickle") 
        pickle.dump(sam2nt2aa, gzip.open(pickleFile, "wb") )
        
        #print samples
        samplePickleFile = os.path.join(self.outdir, 'samples.pickle')
        pickle.dump(samples, gzip.open(samplePickleFile, "wb") )

class RunDrawNt2aaDist(Target):
    """
    """
    def __init__(self, indir, options):
        Target.__init__(self)
        self.indir = indir
        self.options = options

    def run(self):
        samnum = self.options.samnum
        sam2nt2aa = {}
        sam2nt2vals = {}

        for i in xrange(samnum):
            #### HACK ####: print aa -> nt:
            if i == 0:
                samplePickleFile = os.path.join( self.indir, "%d" %i, "samples.pickle")
                samples = pickle.load( gzip.open(samplePickleFile) )
                for sample in samples:
                    fh = open( os.path.join(self.options.outdir, "aa2nts-%s.txt" %sample.name), "w")
                    for header, nt2seq in sample.seqs.iteritems():
                        if len(nt2seq) == 1:
                            continue
                        aa = header.split('|')[0]
                        fh.write(">%s\n" %header)
                        for nt in nt2seq.keys():
                            fh.write("%s\n" %nt)
                fh.close()


            #### END HACK ####

            pickleFile = os.path.join( self.indir, "%d" %i, "sam2nt2aa.pickle" )
            currSam2nt2aa = pickle.load( gzip.open(pickleFile, "rb") )
            for sam, nt2aa in currSam2nt2aa.iteritems():
                if sam not in sam2nt2aa:
                    sam2nt2aa[sam] = nt2aa
                    sam2nt2vals[sam] = {}
                    for nt, aa in nt2aa.iteritems():
                        sam2nt2vals[sam][nt] = [aa]
                else:
                    for nt, aa in nt2aa.iteritems():
                        if nt not in sam2nt2aa[sam]:
                            sam2nt2aa[sam][nt] = aa
                            sam2nt2vals[sam][nt] = [aa]
                        else:
                            sam2nt2aa[sam][nt] += aa
                            sam2nt2vals[sam][nt].append(aa)
        #Normalize
        if samnum > 1:
            for sam, nt2aa in sam2nt2aa.iteritems():
                for nt, aa in nt2aa.iteritems():
                    sam2nt2aa[sam][nt] = aa/samnum
        
        drawDist(sam2nt2aa, self.options)
        
        #print stats to text:
        xdata = []
        for sam, nt2vals in sam2nt2vals.iteritems():
            for nt, aas in nt2vals.iteritems():
                if nt not in xdata:
                    xdata.append(nt)
                if len(aas) < samnum:
                    for i in xrange( samnum - len(aas) ):
                        sam2nt2vals[sam][nt].append(0.0)
        for sam, nt2vals in sam2nt2vals.iteritems():
            for x in xdata:
                if x not in nt2vals:
                    sam2nt2vals[sam][x] = [ 0.0 for i in xrange(samnum) ]
        textfile = os.path.join(self.options.outdir, "nt2aa.txt")
        printNt2aa(sam2nt2vals, sorted(xdata), textfile)

############### DONE JOBTREE ######################
def printNt2aa(sam2nt2vals, xvals, textfile):
    f = open(textfile, 'w')
    f.write( "Sample\t%s\n" %( '\t'.join([ str(x) for x in xvals ]) ) )
    for sam, nt2vals in sam2nt2vals.iteritems():
        f.write("%s" %sam)
        for x in xvals:
            vals = nt2vals[x]
            meanCount = np.mean(vals)
            stdCount = np.std(vals)
            f.write("\t%f (%f)" %(meanCount, stdCount))
        f.write("\n")
    f.close()

#========== DRAWING ===
def getNumNtPerAaSample(sample):
    numNt2aa = {} #key = number of nucleotide sequences, val = number of amino acid from that many nuc. sequences
    for header, nt2seq in sample.seqs.iteritems():
        numNt = len(nt2seq)
        if numNt not in numNt2aa:
            numNt2aa[numNt] = 1
        else:
            numNt2aa[numNt] += 1
    
    #Convert aa count into frequency:
    total = sum( numNt2aa.values() )
    if total > 0:
        for n, a in numNt2aa.iteritems():
            numNt2aa[n] = a*100.0/total

    return numNt2aa

def getNumNtPerAa(samples):
    sample2numNt2aa = {}
    for sample in samples:
        sample2numNt2aa[sample.name] = getNumNtPerAaSample(sample)
    return sample2numNt2aa

def drawDistData(axes, sam2nt2aa):
    #sam2nt2aa = getNumNtPerAa(samples)
    
    lines = []
    labels = []
    colors = iseqlib.getColors6()
    if len(colors) < len(sam2nt2aa.keys()):
        colors.extend( iseqlib.getColors6light() )
    if len(colors) < len(sam2nt2aa.keys()):
        colors.extend( iseqlib.getColors6dark() )
    #colors = iseqlib.getColors0()
    #lightColors = getColors6light()
    markersize = 10.0
    c = -1
    axes.set_yscale('log')
    for s in sorted( sam2nt2aa.keys() ):
        nt2aa = sam2nt2aa[s]
        
        c += 1
        xdata = sorted( nt2aa.keys() )
        ydata = [ nt2aa[x] for x in xdata ]
        
        line = axes.plot(xdata, ydata, color=colors[c], marker='o', markeredgecolor=colors[c], markersize = markersize, linestyle='-', linewidth=2)
        #axes.plot(xdata, ydata, color=lightColors[c], linestyle='-', linewidth=0.5)
        lines.append(line)
        labels.append(s)
        print s
        print xdata
        print ydata
    
    xticks = xrange(0,8,1)
    xticklabels = [ str(x) for x in xticks]
    axes.xaxis.set_ticks(xticks)
    axes.xaxis.set_ticklabels( xticklabels )
    #axes.set_xlim(0.5, 5.5)
    #axes.set_ylim(-0.1, 40)

    axes.set_title('Nucleotide sequences to amino acid sequence', size="xx-large")
    iseqlib.editSpine( axes )
    axes.set_xlabel("Number of nucleotide sequences", size='x-large')
    #axes.set_ylabel("Number of amino acid sequences", size='x-large')
    axes.set_ylabel("Percentage of amino acid sequences", size='x-large')
    legend = axes.legend( lines, labels, numpoints=1, loc='best', ncol=1)
    legend.__drawFrame = False
    axes.yaxis.grid(b=True, color="#BDBDBD", linestyle='-', linewidth=0.005)
    axes.xaxis.grid(b=True, color="#BDBDBD", linestyle='-', linewidth=0.005)

def drawDist(sam2nt2aa, options):
    options.out = os.path.join(options.outdir, 'nt2aa')
    fig, pdf = iseqlib.initImage(10.0, 10.0, options)
    axes = iseqlib.setAxes(fig)
    drawDistData(axes, sam2nt2aa)
    iseqlib.writeImage(fig, pdf, options)

def addOptions( parser ):
    parser.add_option('-i', '--indir', dest='indir', help='Input directory')
    parser.add_option('-o', '--outdir', dest='outdir', help='Output directory')
    parser.add_option('-c', '--mincount', dest='mincount', default=1, type='int', help='Minimum read counts. Default = %default')
    parser.add_option('-s', '--sampling', dest='sampling', default=False, action='store_true', help='If specified, will compute the stats by sampling from the orginal samples. Options --samplingsize and --samplingnum are required if this option is chosen.')
    parser.add_option('--samplingsize', dest='samsize', type='int', help='Size wish to sampling from the sample.')
    parser.add_option('--samplingnum', dest='samnum', type='int', default=1, help='Number of sampling to perform. Default=%default')
    #parser.add_option('-g', '--group2samples', dest='group2samples', help='File specifying which samples belong to which group. Format: <group> <comma,separated,list,of,samples,belong,to,this,group>. If specified, will draw average of each group insteads of each sample separately')
    #parser.add_option('-t', '--ttest', dest='ttest', action='store_true', default=False, help='If specified, will perform ttest')

def main():
    parser = iseqlib.initOptions()
    Stack.addJobTreeOptions( parser )
    iseqlib.initPlotOptions( parser )
    addOptions(parser)
    
    options, args = parser.parse_args()
    iseqlib.checkPlotOptions( options, parser )

    mincount = options.mincount
    #samples = readfiles(options.indir, mincount)

    i = Stack( Setup(options) ).startJobTree(options)
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" %i)
    
    #drawDist(samples, options)
    
    #minsam = 3
    #printSharedSeqFreqAll( samples, minsam, "counts-atleast3sams" )
    #filterByNumSampleAll(samples, minsam, "fasta-atleast3sams")


if __name__ == '__main__':
    from immunoseq.src.overlapNt2aa import *
    main()

