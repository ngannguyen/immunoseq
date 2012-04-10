#!/usr/bin/env python2.6

#nknguyen soe ucsc edu
#Feb 21 2012
#Input: input directory of fasta files
#Output: pair-wise overlaping statistics

import os, sys, re, random, copy
import immunoseq.lib.immunoseqLib as iseqlib
from sonLib.bioio import system

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

def getPc(count, total):
    if total == 0:
        return 0
    return 100.0*count/total

def printPairwiseOverlap(reads1, reads2, clones1, clones2, stats1, stats2, cutoffs, outfile):
    f = open(outfile, "w")
    f.write("#%Cutoff\tClones1\tClones2\tOverlap1\tOverlap2\t%1overlap2\t%2overlap1\t%reads1overlap2\t%reads2overlap1\n")
    for i,c in enumerate(cutoffs):
        oc1 = stats1["oclones"][i]
        or1 = stats1["oreads"][i]
        t1 = clones1[i]
        r1 = reads1[i]
        
        oc2 = stats2["oclones"][i]
        or2 = stats2["oreads"][i]
        t2 = clones2[i]
        r2 = reads2[i]
        
        f.write("%.3f\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n" %( c, t1, t2, oc1, oc2, getPc(oc1, t1), getPc(oc2, t2), getPc(or1, r1), getPc(or2, r2) ))
    f.close()

def getPairwiseOverlap( seqs1, seqs2, outfile, cutoffs, mode, discrete ):
    #Initialize stats:
    stats1 = {"oclones":[], "oreads":[]}
    stats2 = {"oclones":[], "oreads":[]}
    for i in xrange( len(cutoffs) ):
        for k in stats1:
            stats1[k].append(0)
            stats2[k].append(0)

    #get number of clones that pass the cutoff:
    reads1, clones1, total1 = getNumClones(seqs1, cutoffs, discrete)
    reads2, clones2, total2 = getNumClones(seqs2, cutoffs, discrete)

    if total1 == 0 or total2 == 0:
        return
    #get overlap:
    for k in seqs1:
        if k not in seqs2: #not found in repertoire 2
            continue
        s1 = seqs1[k]
        s2 = seqs2[k]
        for i, cutoff in enumerate(cutoffs):
            if mode == 1:
                if s1.freq >= cutoff and s2.freq >= cutoff:
                    if not discrete or i == len(cutoffs) -1 or (discrete and s1.freq < cutoffs[i+1] and s2.freq < cutoffs[i+1]):
                        stats1["oclones"][i] += 1
                        stats1["oreads"][i] += s1.count
                        stats2["oclones"][i] += 1
                        stats2["oreads"][i] += s2.count
            else:
                if s1.freq >= cutoff:
                    if not discrete or i == len(cutoffs) -1 or (discrete and s1.freq < cutoffs[i+1]):
                        stats1["oclones"][i] += 1
                        stats1["oreads"][i] += s1.count
                if s2.freq >= cutoff:
                    if not discrete or i == len(cutoffs) -1 or (discrete and s2.freq < cutoffs[i+1]):
                        stats2["oclones"][i] += 1
                        stats2["oreads"][i] += s2.count

    printPairwiseOverlap(reads1, reads2, clones1, clones2, stats1, stats2, cutoffs, outfile)
    #return reads1, reads2, clones1, clones2, stats1, stats2

def getNumClones(seqs, cutoffs, discrete):
    reads = [ 0 for c in cutoffs ] 
    clones = [ 0 for c in cutoffs ]
    total = sum([ s.count for s in seqs.values() ])
    for s in seqs.values():
        for i,c in enumerate(cutoffs):
            if s.freq >= c:
                if not discrete or i == len(cutoffs) -1 or (discrete and s.freq < cutoffs[i + 1]):
                    clones[i] += 1
                    reads[i] += s.count
    return reads, clones, total

################### PRINT SHARED SEQUENCES ##############################
def printSharedSeqSummary(outdir, name1, name2, seqs1, seqs2):
    outfile = os.path.join(outdir, "%s-%s" %(name1, name2))
    f = open(outfile, 'w')
    f.write("#%s\n"%name1)
    for seq in sorted(seqs1, key=lambda s:s.freq, reverse=True):
        f.write("%s\t%f\n" %(seq.seq, seq.freq))
    f.write("#%s\n"%name2)
    for seq in sorted(seqs2, key=lambda s:s.freq, reverse=True):
        f.write("%s\t%f\n" %(seq.seq, seq.freq))
    f.close()

def printSharedFasta(outdir, name1, name2, seqs):
    outdir = os.path.join(outdir, name1)
    system("mkdir -p %s" %outdir)
    outfile = os.path.join(outdir, "%s_%s.fa" %(name1, name2))
    f = open(outfile, 'w')
    seqs = sorted(seqs, key=lambda s: s.count, reverse=True)
    for i, s in enumerate(seqs):
        f.write(">%s_%s;%d|%s|%s;size=%d\n" %(name1, name2, i, ','.join(s.vs), ','.join(s.js), s.count) )#name;id|vs|js;size=###
        f.write("%s\n" %s.seq)
    f.close()

def printPairwiseOverlapSequences(name1, name2, seqs1, seqs2, outdir, cutoffs, mode, discrete):
    cutoff2seqs1 = {}
    cutoff2seqs2 = {}
    #Initialize cutoff2seqs:
    for c in cutoffs:
        cutoff2seqs1[c] = []
        cutoff2seqs2[c] = []

    for k in seqs1:
        if k not in seqs2:
            continue
        s1 = seqs1[k]
        s2 = seqs2[k]
        for i, cutoff in enumerate(cutoffs):
            if mode == 1:
                if s1.freq >= cutoff and s2.freq >= cutoff:
                    if not discrete or i == len(cutoffs) -1 or (discrete and s1.freq < cutoffs[i+1] and s2.freq < cutoffs[i+1]):
                        cutoff2seqs1[cutoff].append(s1)
                        cutoff2seqs2[cutoff].append(s2)
            else:
                if s1.freq >= cutoff:
                    if not discrete or i == len(cutoffs) -1 or (discrete and s1.freq < cutoffs[i+1]):
                        cutoff2seqs1[cutoff].append(s1)
                if s2.freq >= cutoff:
                    if not discrete or i == len(cutoffs) -1 or (discrete and s2.freq < cutoffs[i+1]):
                        cutoff2seqs2[cutoff].append(s2)
    #Print the sequences
    for c in cutoffs:
        cOutdir = os.path.join(outdir, "%.3f" %c)
        system("mkdir -p %s" %cOutdir)
        printSharedFasta(cOutdir, name1, name2, cutoff2seqs1[c])
        printSharedFasta(cOutdir, name2, name1, cutoff2seqs2[c])
        #Print frequencies of shared sequences:
        freqdir = os.path.join(outdir, "freq", "%.3f" %c)
        system("mkdir -p %s" %freqdir)
        printSharedSeqSummary(freqdir, name1, name2, cutoff2seqs1[c], cutoff2seqs2[c])
  

################### END PRINT SHARED SEQUENCES ##############################

################## PAIRWISE OVERLAP PLOT (ROBINS et al.) #####################

def drawOverlapPlotData(data, axes, labels):
    colors = iseqlib.getColors0()
    lines = []
    for i, d in enumerate(data):
        xdata = d[0]
        ydata = d[1]
        l = axes.plot(xdata, ydata, color=colors[i], linestyle='-')
        lines.append(l)
    legend = pyplot.legend(lines, labels, numpoints=1, loc='best')

def drawOverlapPlot(data, file, labels, options):
    options.out = file
    fig, pdf = iseqlib.initImage(8.0, 10.0, options)
    axes = iseqlib.setAxes(fig)
    drawOverlapPlotData(data, axes, labels)
    iseqlib.writeImage(fig, pdf, options)


#def pairwiseOverlapPlot(seqs1, seqs2, outfile): 
def pairwiseOverlapPlot(seqs1, seqs2): 
    #sorting seqs1
    seqs1list = [(header, seq) for header, seq in seqs1.iteritems()]
    seqs1list = sorted( seqs1list, key=lambda item:item[1].count, reverse=True )
    #sorting seqs2
    seqs2list = [(header, seq) for header, seq in seqs2.iteritems()]
    seqs2list = sorted( seqs2list, key=lambda item:item[1].count, reverse=True )
    
    top1 = {}
    top2 = {}
    xdata = []
    ydata = []
    currx = 0
    curry = 0
    
    for i in xrange( min([len(seqs1list), len(seqs2list)]) ):
        (h1, s1) = seqs1list[i]
        (h2, s2) = seqs2list[i]
        #n1 = s1.count
        #n2 = s2.count
        top1[h1] = s1
        top2[h2] = s2
        if h1 in top2: #current sequence1 is also in sequences 2
            #curry += min( [s1.count, top2[h1].count] )
            curry += 1
        if h1 != h2 and h2 in top1: #current sequence2 is also in sequences 1
            #curry += min( [s2.count, top1[h2].count] )
            curry += 1
        #currx += n1*n2
        currx = i*i
        xdata.append( currx )
        ydata.append(curry)
    
    #write to output file
    #f = open(outfile, 'w')
    #for i, x in enumerate(xdata):
    #    f.write("%d\t%d\n" %(x, ydata[i]))
    #f.close()
    #return
    return xdata, ydata

################## END PAIRWISE OVERLAP PLOT (ROBINS et al.) #####################

################# GET FASTA FILES OF UNIQUE SEQUENCES, SHARED WITH SAME GROUP SEQUENCES, SHARED WITH OTHER GROUPS SEQUENCES, SHARED SEQUENCES ###################
def getSeq2sample2count(samples): #key = seq|vs|js, val = {sample: (count, freq)}
    seq2sample2count = {}
    for sample in samples:
        for header, seq in sample.seqs.iteritems():
            if header not in seq2sample2count:
                seq2sample2count[header] = {sample.name: (seq.count, seq.freq)}
            else:
                seq2sample2count[header][sample.name] = (seq.count, seq.freq)
    return seq2sample2count

def uniqAndSharedSequences(outdir, sample, seq2sample2count, group2samples, group):
    uf = open( os.path.join(outdir, 'uniq.fa'), 'w' )
    sgf = open( os.path.join(outdir, 'sameGroupShared.fa'), 'w' )
    dgf = open( os.path.join(outdir, 'diffGroupShared.fa'), 'w' )
    sf = open( os.path.join(outdir, 'shared.fa'), 'w' )
    
    total = sum([ seq.count for seq in sample.seqs.values()])
    ucount = 0
    sgcount = 0
    dgcount = 0
    scount = 0

    for header, seq in sample.seqs.iteritems():
        if header not in seq2sample2count:
            raise ValueError("Sequence %s of sample %s is not in seq2sample2count\n" %(header, sample.name))
        sample2count = seq2sample2count[header]
        fastaheader = "%s;freq=%.4f|%s|%s;size=%d" %(sample.name, seq.freq, ','.join(seq.vs), ','.join(seq.js), seq.count)
        samegroup = False
        diffgroup = False
        
        for g, samples in group2samples.iteritems():
            if g == group: #same group
                for s in samples:
                    if s != sample.name and s in sample2count:
                        samegroup = True
                        break
            else: #diff group
                for s in samples:
                    if s != sample.name and s in sample2count:
                        diffgroup = True
                        break
        if samegroup or diffgroup:
            sf.write(">%s\n" % fastaheader)
            sf.write("%s\n" % seq.seq)
            scount += seq.count
            if samegroup and not diffgroup:
                sgf.write(">%s\n" % fastaheader)
                sgf.write("%s\n" % seq.seq)
                sgcount += seq.count
            if diffgroup:
                dgf.write(">%s\n" % fastaheader)
                dgf.write("%s\n" % seq.seq)
                dgcount += seq.count
        else: #unique sequence
            uf.write(">%s\n" % fastaheader)
            uf.write("%s\n" % seq.seq)
            ucount += seq.count

    uf.write(">%s;others||;size=%d\n%s\n" %(sample.name, total - ucount, "OTHERS"))
    sgf.write(">%s;others||;size=%d\n%s\n" %(sample.name, total - sgcount, "OTHERS"))
    dgf.write(">%s;others||;size=%d\n%s\n" %(sample.name, total - dgcount, "OTHERS"))
    sf.write(">%s;others||;size=%d\n%s\n" %(sample.name, total - scount, "OTHERS"))

    uf.close()
    sgf.close()
    dgf.close()
    sf.close()

################# END OF GET FASTA FILES OF UNIQUE SEQUENCES, SHARED WITH SAME GROUP SEQUENCES, SHARED WITH OTHER GROUPS SEQUENCES, SHARED SEQUENCES ###################


################# CONVERTING SEQS (SORTED LIST) TO A (BALANCED) BINARY TREE ################
#Converting the sorted list into a balanced binary tree
#def insertMiddleNode(bt, seqs, leftindex, rightindex):
#    middleindex = (leftindex + rightindex)/2
#    btnode = BinaryTreeNode( seqs[middleindex] )
#    bt.insert(btnode)
#
#def convertSeqsToBinaryTree(seqs):
#    bt = iseqlib.BinaryTree()
#    leftindex = 0
#    rightindex = len(seqs) -1
#     
#    while True:
#        middleindex = (leftindex + rightindex)/2 
#        bt.insert( seqs[middleindex] )
#
#
#    return bt


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

def readGroup2samples(file):
    group2samples = {}
    sample2group = {}
    f = open(file, 'r')
    for line in f:
        items = line.strip().split()
        samples = items[1].split(',')
        group2samples[items[0]] = samples
        for s in samples:
            sample2group[s] = items[0]
    f.close()
    return group2samples, sample2group

def checkOptions(parser, args, options):
    if options.indir == '':
        parser.error("Input directory was not specified\n")
    if options.outdir == '':
        parser.error("Output directory was not specified\n")
    if options.mode != 1 and options.mode !=2:
        parser.error("Invalid mode %d. Please select either 1 or 2\n" %options.mode)
    options.cutoffs = [ float(p) for p in options.minPercentage.strip(',').split(',') ]

def addOptions(parser):
    parser.add_option('-i', '--indir', dest='indir', help="Required argument. Input directory")
    parser.add_option('-o', '--outdir', dest='outdir', help="Required argument. Output directory")
    parser.add_option('-p', '--minPercentage', dest='minPercentage', default="0,0.001,0.01,0.05,0.1,0.5,1,5", help='Comma separated list of percentage cutoffs. Clone must have at least this percentage of total reads to be included. Default=%default')
    parser.add_option('-m', '--mode', dest='mode', type='int', default=1, help="Mode 1: for each sample, only look at clones that pass the minCount and minPercentage. Mode 2: clones must pass minCount and minPercentage in at least one sample. Default=1")
    parser.add_option('-c', '--minCount', dest='minCount', default=1, type='int', help='Minimun clone count to be included. Default=%default')
    parser.add_option('-d', '--discrete', dest='discrete', default=False, action="store_true", help='If specified, calculate stats discretely, not cumulatively; i.e looks at clones lie within a range from A to B frequencies. Default = %default')
    parser.add_option('-s', '--sampling', dest='sampling', default=-1, type='int', help='Normalized size wished to look at. Default=-1, meaning no normalization')
    parser.add_option('-u', '--samplingUniq', dest='samplingUniq', action='store_true', default=False, help='If specified, sampling on the set of unique sequences instead of all sequences')
    parser.add_option('-f', '--fasta', dest='fasta', action='store_true', default=False, help='If specified, will produce fasta files of shared sequences for each pair of samples')
    parser.add_option('--noPairwiseOverlap', dest='noPairwiseOverlap', action='store_true', default=False, help='Calculate pairwise shared sequences. Default = %default')
    parser.add_option('--overlapPlot', dest='overlapPlot', action='store_true', default=False, help='Default=%default')
    parser.add_option('--group2samples', dest='group2samples', help='Only specified if wish to print out, for each sample, the uniq sequences, sequences that are shared with at least another sample in the same group, and sequences that are shared with at least another sample of another group, and sequences that are shared with at least another sample, regardless of the group')

def main():
    parser = iseqlib.initOptions()
    addOptions(parser)
    iseqlib.initPlotOptions(parser)
    options, args = parser.parse_args()
    checkOptions(parser, args, options)
    iseqlib.checkPlotOptions(options, parser)

    samples = readfiles(options.indir, options.minCount) 
    if options.sampling > 0:
        sampling(samples, options.sampling, options.samplingUniq)
    
    if options.overlapPlot:
        overlapPlotDir = os.path.join(options.outdir, "overlapPlots")
        system("mkdir -p %s" %(overlapPlotDir))

    if not options.noPairwiseOverlap:
        statdir = os.path.join(options.outdir, "overlapStats")
        system("mkdir -p %s" %statdir)
    
    #Pairwise overlap statistics and sequences
    for i in xrange( len(samples) - 1 ):
        sample1 = samples[i]
        for j in xrange(i+1, len(samples)):
            sample2 = samples[j]
            pair = "%s_%s" %(sample1.name, sample2.name)
            if not options.noPairwiseOverlap:
                outfile = os.path.join(statdir, "%s.txt" %pair)
                getPairwiseOverlap(sample1.seqs, sample2.seqs, outfile, options.cutoffs, options.mode, options.discrete)
            if options.fasta:
                fastaDir = os.path.join(options.outdir, "pairwiseSharedSeqs")
                system("mkdir -p %s" %fastaDir)
                printPairwiseOverlapSequences(sample1.name, sample2.name, sample1.seqs, sample2.seqs, fastaDir, options.cutoffs, options.mode, options.discrete)
                
    #HARLAN ROBINS OVERLAP PLOT (Y axis = number of shared sequences, X axis = n1j x n2j where nij (i=1,2) is the number of sequences of sample i in the top jth clones. 
    if options.overlapPlot:
        for i in xrange( len(samples) -1 ):
            sample1 = samples[i]
            data = []
            labels = []
            for j in xrange(len(samples) -1):
                if i == j:
                    continue
                sample2 = samples[j]
                #Robins overlap plot:
                xdata, ydata = pairwiseOverlapPlot(sample1.seqs, sample2.seqs)
                data.append( [xdata, ydata] )
                labels.append(sample2.name)
                ofile = os.path.join(overlapPlotDir, '%s-%s' %(sample1.name, sample2.name) )
                drawOverlapPlot([ [xdata, ydata] ], ofile, [ "%s-%s" %(sample1.name, sample2.name)], options)
            #ofile = os.path.join(overlapPlotDir, '%s' %(sample1.name) )
            #drawOverlapPlot(data, ofile, labels, options)

    #====== IF GROUP2SAMPLEs ========
    if options.group2samples:
        group2samples, sample2group = readGroup2samples(options.group2samples)
        seq2sample2count = getSeq2sample2count(samples) #key = seq, val = {sample: (count, freq)}
        for sample in samples:
            outdir = os.path.join(options.outdir, "uniqAndSharedSeqs", sample.name)
            system("mkdir -p %s" %outdir)
            uniqAndSharedSequences(outdir, sample, seq2sample2count, group2samples, sample2group[sample.name])


if __name__ == '__main__':
    main()

