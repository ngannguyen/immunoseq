#!/usr/bin/env python2.6

'''
05/01/2012 nknguyen soe ucsc edu
Input: input directory of fasta files
    Each fasta file is a concatenation of all samples fasta files

Output: Plot showing number of shared sequences versus number of samples containing the sequences for Each input fasta file
'''

import os, sys, re, copy
import matplotlib.pyplot as pyplot
import immunoseq.lib.immunoseqLib as iseqlib
from sonLib.bioio import system

class Sample():
    def __init__(self, name):
        self.name = name
        self.subsamples = []
        #self.seqs = []
        #self.total = 0

class Seq():
    def __init__(self, name, seq, count, vs, js, freq, translate):
        self.name = name
        self.samples = name.split(',')
        self.seq = seq
        self.count = count
        self.freq = freq
        self.vs = sorted(vs)
        self.js = sorted(js)
        #vstr = ','.join(vs)
        #jstr = ','.join(js)
        vfams = getGeneFamilies(self.vs)
        #jfams = getGeneFamilies(self.js)
        vstr = ','.join(vfams)
        jstr = ','.join(js)
        
        if translate:
            seq = iseqlib.nt2aa(seq) 
        self.header =  '|'.join([seq, vstr, jstr]) #header example: CASSLRRGGKPGELFF|TRBV5|TRBJ2-2  (aminoacid|vfamilies|jgenes)
    
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

def getGeneFamilies(genes):
    fams = []
    for gene in genes:
        fam = gene.split('-')[0]
        if fam not in fams:
            fams.append(fam)
    return fams

def getSeq(header, seq, translate):
    items = header.split(";size=")
    if len(items) < 2:
        raise ValueError("Wrong header format: %s\n" %header)
    count = int(items[-1])
    #Find v and j:
    name = items[0].split(';')[0]
    genestr = items[0].split(';')[-1]

    freq = 0.0
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
    
    seq =  Seq(name, seq, count, sorted(vs), sorted(js), freq, translate)
    return seq

def addSeq( header, seq, seqs, minCount, translate ):
    seq = getSeq(header, seq, translate)
    if seq.count >= minCount:
        if seq.header not in seqs:
            seqs[seq.header] = {seq.name: seq}
        else:
            seqs[seq.header][seq.name] = seq

def getSubSamples(seqs):
    samples = []
    for header, name2seq in seqs.iteritems():
        for name in name2seq.keys():
            if name not in samples:
                samples.append(name)
    return samples

def readFile(file, minCount, translate):
    seqs = {} #key = sequence, val = Seq
    f = open(file, 'r')
    header = ''
    seq = ''
    for line in f:
        line = line.strip()
        if line[0] == '>':
            if header != '' and seq != '':
                addSeq( header, seq, seqs, minCount, translate )
            header = line.lstrip('>')
            seq = ''
        else:
            seq = line
    if seq != '' and header != '':
        addSeq( header, seq, seqs, minCount, translate )
    f.close()
    subsamples = getSubSamples(seqs)
    return seqs, subsamples

def setFreq(sample):
    '''Set the clone frequencies, frequencies are relative to each subsamples
    '''
    for subsample in sample.subsamples:
        total = 0
        for header, name2seq in sample.seqs.iteritems():
            if subsample in name2seq:
                total += name2seq[subsample].count
        for header, name2seq in sample.seqs.iteritems():
            if subsample in name2seq:
                name2seq[subsample].setFreq(total)
    return

def readfiles(indir, minCount, translate):
    files = os.listdir( indir )
    samples = []
    for file in files:
        path = os.path.join(indir, file)
        if os.path.isdir(path) or file.split('.')[-1] != 'fa':
            continue
        sampleName = file.split('.')[0]
        sample = Sample(sampleName)
        seqs, subsamples = readFile( path, minCount, translate )

        #Set frequencies
        #total = sum( [s.count for s in seqs.values()] )
        #sample.total = total
        #for s in seqs.values():
        #    s.setFreq(total)
        #setFreq(sample)
        sample.seqs = seqs
        sample.subsamples = sorted(subsamples)
        samples.append(sample)
    return samples

def sharedSeqFreq(sample):
    numsam2header2freqs = {}
    for header, name2seq in sample.seqs.iteritems():
        numsam = len( name2seq.keys() )
        #freqs = sorted( [ seq.freq for seq in name2seq.values() ] )
        n2freq = {}
        for n, s in name2seq.iteritems():
            #n2freq[n] = s.freq
            n2freq[n] = s.count #HACK

        if numsam not in numsam2header2freqs:
            numsam2header2freqs[numsam] = {header:n2freq}
        else:
            numsam2header2freqs[numsam][header] = n2freq
    return numsam2header2freqs

def printSharedSeqFreq(sample, minsam, outfile):
    f = open(outfile, 'w')
    numsam2header2freqs = sharedSeqFreq(sample)
    for numsam, header2freqs in numsam2header2freqs.iteritems():
        if numsam >= minsam:
            for header, n2f in header2freqs.iteritems():
                #freqs = n2f.values()
                #f.write("%s\t%s\n" %(header, ','.join( [str(freq) for freq in freqs] )) )
                names = []
                freqs = []
                for n, freq in n2f.iteritems():
                    names.append(n)
                    freqs.append( str(freq) )
                f.write("%s\t%s\t%s\n" % (header, ','.join(names) , ','.join(freqs) ))

    f.close()
    return

def printSharedSeqFreqAll( samples, minsam, outdir ):
    for s in samples:
        file = os.path.join( outdir, "%s-freqs.txt" %s.name )
        printSharedSeqFreq(s, minsam, file)

def filterByNumSample(sample, minsam, outfile, minfreq):
    #sample = calcFreq(sample) #HACK
    f = open(outfile, 'w')
    for header, name2seq in sample.seqs.iteritems():
        if len(name2seq.keys()) >= minsam:
            for name, seq in name2seq.iteritems():
                if seq.freq < minfreq:
                    continue
                h = seq.getFastaHeader()
                f.write("%s\n" %h)
                f.write("%s\n" %seq.seq)
    f.close()
    return

def filterByNumSampleAll(samples, minsam, outdir, minfreq):
    for s in samples:
        outfile = os.path.join(outdir, "%s.fa" % s.name)
        filterByNumSample(s, minsam, outfile, minfreq)

#========== DRAWING ===
def getSharedSeqDist(samples):
    sample2dist = {}
    for sample in samples:
        #sample = calcFreq(sample) #HACK
        numsam2count = {}
        for header, name2seq in sample.seqs.iteritems():
            #numsam = len( name2seq.values()[0].samples )
            numsam = len( name2seq.keys() ) #number of samples having this sequence
            #if numsam == 1:
            #    continue
            #    print sample.name
            #    print header
            #    print name2seq[ name2seq.keys()[0] ].name
            #    exit(1)
            if numsam not in numsam2count:
                numsam2count[numsam] = 1
            else:
                numsam2count[numsam] += 1
        sample2dist[sample.name] = numsam2count
    return sample2dist

def drawDistData(axes, samples):
    sample2dist = getSharedSeqDist(samples)
    
    lines = []
    labels = []
    colors = iseqlib.getColors6()
    #lightColors = getColors6light()
    markersize = 10.0
    c = -1
    for s in sorted( sample2dist.keys() ):
        numsam2count = sample2dist[s]
        
        c += 1
        xdata = sorted( numsam2count.keys() )
        ydata = [ numsam2count[x] for x in xdata ]
        totaly = sum(ydata)
        pcydata = [(100.0*y)/totaly for y in ydata]
        
        #line = axes.plot(xdata, ydata, color=colors[c], marker='o', markeredgecolor=colors[c], markersize = markersize, linestyle='-', linewidth=2)
        line = axes.plot(xdata, pcydata, color=colors[c], marker='o', markeredgecolor=colors[c], markersize = markersize, linestyle='-', linewidth=2)
        #axes.plot(xdata, ydata, color=lightColors[c], linestyle='-', linewidth=0.5)
        lines.append(line)
        labels.append(s)
        print s
        print xdata
        print ydata
        print pcydata
    
    axes.set_yscale('log')
    #axes.set_xlim(8.5, 23.5)
    #xticks = xrange(8, 24, 2)
    #xticklabels = [ str(x) for x in xticks]
    #axes.xaxis.set_ticks(xticks)
    #axes.xaxis.set_ticklabels( xticklabels )

    axes.set_title('Shared sequences', size="xx-large")
    iseqlib.editSpine( axes )
    axes.set_xlabel("Number of samples", size='x-large')
    axes.set_ylabel("Number of clones", size='x-large')
    legend = axes.legend( lines, labels, numpoints=1, loc='best', ncol=1)
    legend.__drawFrame = False
    axes.yaxis.grid(b=True, color="#BDBDBD", linestyle='-', linewidth=0.005)
    axes.xaxis.grid(b=True, color="#BDBDBD", linestyle='-', linewidth=0.005)

def drawDist(samples, options):
    options.out = os.path.join(options.outdir, "sharedSeqsDist")
    fig, pdf = iseqlib.initImage(10.0, 10.0, options)
    axes = iseqlib.setAxes(fig)
    drawDistData(axes, samples)
    iseqlib.writeImage(fig, pdf, options)

######### CLONE MATRIX ###############
def calcFreq(sample):
    newsample = copy.copy(sample)
    name2total = {}
    for header, name2seq in sample.seqs.iteritems():
        for name, seq in name2seq.iteritems():
            if name in name2total:
                name2total[name] += seq.count
            else:
                name2total[name] = seq.count
    for header, name2seq in sample.seqs.iteritems():
        for name, seq in name2seq.iteritems():
            seq.setFreq( name2total[name] )
    return newsample

def getTopClones(sample, num):
    sample2topheaders = {}
    for subsample in sample.subsamples:
        seqs = []
        for header, name2seq in sample.seqs.iteritems():
            if subsample in name2seq:
                seqs.append(name2seq[subsample])
        seqs = sorted(seqs, key=lambda seq: seq.count, reverse=True)
        if len(seqs) < num:
            sample2topheaders[subsample] = [s.header for s in seqs]
        else:
            sample2topheaders[subsample] = [s.header for s in seqs[:num]]
    return sample2topheaders

def printCloneMatrix(outfile, sample, minsam, freqTransform):
    sample = calcFreq(sample)
    f = open(outfile, 'w')
    f.write( "Clone\t%s\n" % '\t'.join(sample.subsamples) )
    
    #HACK
    #top clones:
    num=20
    s2topheaders = getTopClones(sample, num) 
    topheaders = []
    for s, hs in s2topheaders.iteritems():
        for h in hs:
            if h not in topheaders:
                topheaders.append(h)
    #END HACK

    #for header, name2seq in sample.seqs.iteritems():
    for header in topheaders:
        name2seq = sample.seqs[header]
        rowvals = []
        if len(name2seq) < minsam:
            continue
        for name in sample.subsamples:
            count = 0
            if name in name2seq:
                count = name2seq[name].count
                if freqTransform:
                    count = int ( name2seq[name].freq*freqTransform + 0.5 )
            rowvals.append( str(count) )
        f.write("%s\t%s\n" %(header, "\t".join(rowvals)))
    f.close()

def printCloneMatrixAll(outdir, samples, minsam, freqTransform):
    for sample in samples:
        outfile = os.path.join(outdir, "cloneM-%s.txt" %sample.name)
        printCloneMatrix(outfile, sample, minsam, freqTransform)

def addOptions( parser ):
    parser.add_option('-i', '--indir', dest='indir', help='Input directory')
    parser.add_option('-o', '--outdir', dest='outdir', help='Output directory')
    parser.add_option('-f', '--freq', dest='freq', type='float', default=0.0, help='Minimum frequency. Default=%default')
    parser.add_option('-c', '--count', dest='count', type='int', default=1, help='Minimum read count. Default=%default')
    parser.add_option('-s', '--minsam', dest='minsam', type='int', default=1, help='Minimum number of samples. Default=%default')
    parser.add_option('-m', '--clonematrix', dest='clonematrix', action='store_true', default=False, help='If specified, print out the clonesize matrix where rows are clones, columns are samples, and each cell is the read-count of corresponding clone of corresponding sample. Default=%default')
    parser.add_option('--freqTransform', dest='freqTransform', type='float', help='If specified, the cell values of clonematrix is the transform of frequency*freqTransform. Default= no transform, just print out read counts.')
    parser.add_option('-d', '--drawdist', dest='drawdist', action='store_true', default=False, help='If specified, draw distribution of number of clones versus number of samples shared. Default=%default')
    parser.add_option('--fasta', dest='fasta', action='store_true', default=False, help='If specified, print fasta file of sequences that are shared by at least "minsam" number of samples. Default=%default')
    parser.add_option('-t', '--translate', dest='translate', action='store_true', default=False)

def main():
    parser = iseqlib.initOptions()
    iseqlib.initPlotOptions( parser )
    addOptions(parser)
    
    options, args = parser.parse_args()
    iseqlib.checkPlotOptions( options, parser )

    samples = readfiles(options.indir, options.count, options.translate)

    if options.drawdist:
        drawDist(samples, options) 
    
    #printSharedSeqFreqAll( samples, minsam, "counts-atleast3sams" )
    if options.fasta:
        faoutdir = os.path.join(options.outdir, "fasta-atleast%dsams" %options.minsam)
        system("mkdir -p %s" %faoutdir)
        filterByNumSampleAll(samples, options.minsam, faoutdir, options.freq)
    
    if options.clonematrix:
        printCloneMatrixAll(options.outdir, samples, options.minsam, options.freqTransform)



if __name__ == '__main__':
    main()

