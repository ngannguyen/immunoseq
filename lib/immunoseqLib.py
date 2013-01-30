#nknguyen soe ucsc edu
#Library of common functions used in the immunoseq pipeline

import os, re, sys
import copy as cp
import random as rnd
from numpy import *
from matplotlib.ticker import LogLocator
from optparse import OptionParser
#from sonLib.bioio import system

###########################################################################
################## READ FASTA FILES WITH HEADER HAS V,J, SIZE= FORMAT #####
###########################################################################
class Sample():
    def __init__(self, name):
        self.name = name
        #self.seqs = []
        #self.total = 0

    def setTotal(self, total):
        self.total = total

class Seq():
    def __init__(self, seq, count, vs, js, ds):
        self.seq = seq
        self.count = count
        self.freq = -1
        self.vs = sorted(vs)
        self.js = sorted(js)
        self.ds = sorted(ds)
        self.vstr = ','.join(self.vs)
        self.jstr = ','.join(self.js)
        self.dstr = ','.join(self.ds)
        self.header =  '|'.join([seq, self.vstr, self.jstr, self.dstr]) #header example: CASSLRRGGKPGELFF|TRBV5-4|TRBJ2-2|TRBD1-1

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
    #name = items[0].split(';')[0]
    #genestr = items[0].split(';')[-1]
    genestr = items[0].split(';')[1]
    clusters = genestr.split(',,')
    vs = []
    js = []
    ds = [] 
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

        if len(citems) == 4:
            d = citems[3]
            if d != '(undefined)' and d not in ds:
                ds.append(d)
    
    return Seq(seq, count, vs, js, ds)

def addSeq( header, seq, seqs, minCount, mode ):
    #mode == 1: each sample: seqs is { header: Seq }
    #mode == 2: seqs is {aaSeq : { v: {j : Seq} } }
    seq = getSeq(header, seq)
    if seq.count < minCount:
        return 0

    if mode == 1:
        if seq.header not in seqs:
            seqs[seq.header] = seq
        else:
            seqs[seq.header].count += seq.count
        return seqs[seq.header].count

    #Mode 2:
    j2seq = {}
    for j in seq.js:
        j2seq[j] = seq
    v2j = {}
    for v in seq.vs:
        v2j[v] = j2seq

    if seq.seq not in seqs:
        seqs[seq.seq] = v2j
    else:
        currv2j = seqs[seq.seq]
        for v in seq.vs:
            currv2j[v] = j2seq
    return seq.count

def readFile(file, minCount, mode):
    #mode == 1: each sample: seqs is { header: Seq }
    #mode == 2: seqs is {aaSeq : { v: {j : Seq} }}
    seqs = {} 
    f = open(file, 'r')
    header = ''
    seq = ''
    total = 0
    for line in f:
        line = line.strip()
        if len(line) == 0:
            continue
        if line[0] == '>':
            if header != '' and seq != '':
                count = addSeq( header, seq, seqs, minCount, mode )
                total += count
            header = line.lstrip('>')
            seq = ''
        else:
            seq = line
    if seq != '' and header != '':
        count = addSeq( header, seq, seqs, minCount, mode )
        total += count
    f.close()
    return seqs, total

def readfiles(indir, minCount, mode):
    #mode == 1: each sample: seqs is { header: Seq }
    #mode == 2: seqs is {aaSeq : { v: {j : Seq} }}
    ext = "fa"
    files = getfiles(indir, ext)
    samples = []
    for file in files:
        path = os.path.join(indir, file)
        sampleName = file.split('.')[0]
        sample = Sample(sampleName)
        seqs, total = readFile( path, minCount, mode )

        #Set frequencies
        if mode == 1:
            for s in seqs.values():
                s.setFreq(total)
        else:
            for aa, v2j in seqs.iteritems():
                for v, j2seq in v2j.iteritems():
                    for j, seq in j2seq.iteritems():
                        seq.setFreq(total)

        sample.seqs = seqs
        sample.setTotal(total)
        samples.append(sample)
    return samples

#######
def hasSeq(seq, aa2v2j):
    if seq.seq in aa2v2j:
        v2j = aa2v2j[seq.seq]
        for v in seq.vs:
            if v in v2j:
                j2seq = v2j[v]
                for j in seq.js:
                    if j in j2seq:
                        return True
    return False

def findSeq(seq, aa2v2j):
    if seq.seq in aa2v2j:
        v2j = aa2v2j[seq.seq]
        for v in seq.vs:
            if v in v2j:
                j2seq = v2j[v]
                for j in seq.js:
                    if j in j2seq:
                        return j2seq[j]
    return None 

def findSeqs(seq, aa2v2j):
    seqs = []
    if seq.seq in aa2v2j:
        v2j = aa2v2j[seq.seq]
        for v in seq.vs:
            if v in v2j:
                j2seq = v2j[v]
                for j in seq.js:
                    if j in j2seq:
                        currseq = j2seq[j]
                        if currseq not in seqs:
                            seqs.append(currseq)
    return seqs 

def getsample2aa2v2j(samples):
    sample2aa2v2j = {}
    for sample in samples:
        aa2v2j = {}
        seqs = sample.seqs
        for header, seq in seqs.iteritems():
            aa = seq.seq
            v2j = {}
            if aa in aa2v2j:
                v2j = aa2v2j[aa]

            for v in seq.vs:
                j2seq = {}
                if v in v2j:
                    j2seq = v2j[v]
                for j in seq.js:
                    j2seq[j] = seq
                v2j[v] = j2seq
            aa2v2j[aa] = v2j

        sample2aa2v2j[sample.name] = aa2v2j
    return sample2aa2v2j
#######

######### FILTERING BY GENES #####################
def filterSampleByGenes(sample, vs, js):
    if not vs and not js:
        return sample
    newsample = Sample(sample.name)
    seqs = {}
    for header, seq in sample.seqs.iteritems():
        if not vs:#if not selected v list, include all Vs
            checkV = True
        else:#check to see if has at least 1 gene in the selected V list
            checkV = False
            for v in seq.vs:
                if v in vs:
                    checkV = True
                    break

        if not js:
            checkJ = True
        else:
            checkJ = False
            for j in seq.js:
                if j in js:
                    checkJ = True
                    break
        if checkV and checkJ:
            seqs[header] = seq
    
    #set frequencies
    total = sum([s.count for s in seqs.values()])
    if total > 0:
        for s in seqs.values():
            s.setFreq(total)
    newsample.seqs = seqs
    return newsample

######### END FILTERING BY GENES #####################


#======= SAMPLING =============
#Different samples have different level of sequencing, and hence will create bias. For example, large sample (lots of sequences got amplified and sequenced) will inherently have more overlapping with other samples than a smaller sample. THe purpose of sampling is to make sure that every sample has the same number of starting sequences to avoid the bias
def samplingSample_uniq(sample, size):
    newsample = Sample(sample.name)
    newseqs = {}
    sys.stderr.write("Sampling sample: %s with %d uniq sequences. Sampling size: %d\n" %(sample.name, len(sample.seqs.keys()), size))
    chosenHeaders = rnd.sample( sample.seqs.keys(), size)
    for header in chosenHeaders:
        newseqs[header] = cp.copy( sample.seqs[header] )
    total = sum( [s.count for s in newseqs.values()] )
    for s in newseqs.values():
        s.setFreq(total)
    newsample.seqs = newseqs
    return newsample

#def samplingSample_weightedUniq(sample, size, numDraws):
#    sys.stderr.write("Sampling sample: %s with %d uniq sequences. Sampling size: %d\n" %(sample.name, len(sample.seqs.keys()), size))
#    newsample = Sample(sample.name)
#    newseqs = {}
#    seqs = sorted( sample.seqs.values(), reverse=True, key=lambda seq: seq.freq )
#    if len(seqs) < size:
#        return
#    freqs = [ seq.freq/100.0 for seq in seqs ]
#    if sum(freqs) > 1:
#        sys.stderr.write("Sample %s, sum of freqs is %f > 1\n" %(sample.name, sum(freqs) ))
#
#    #numDraws = max( [sum([ seq.count for seq in seqs ]), len(seqs)*2] )
#    counts = np.random.multinomial( numDraws, freqs )
#    
#    numseq = 0
#    for i, c in enumerate(counts):
#        if c > 0:
#            seq = Seq( seqs[i].seq, c, seqs[i].vs, seqs[i].js )
#            newseqs[seq.header] = seq
#            numseq += 1
#        if numseq >= size:
#            break
#    total = sum( [s.count for s in newseqs.values()] )
#    for s in newseqs.values():
#        s.setFreq(total)
#    newsample.seqs = newseqs
#    return newsample

def samplingSample_weightedUniq(sample, size):
    sys.stderr.write("Begin to sample %d uniq sequences from sample %s\n" %(size, sample.name))
    seqs = sample.seqs
    if len(seqs) < size:
        return
    
    newsample = Sample(sample.name)
    headers = seqs.keys()
    indexList = []#index list represents all the sequences of the sample
    for i, h in enumerate(headers):
        s = seqs[h]
        indexList.extend( [i for j in xrange(s.count)] )
    sys.stderr.write("\tDone with the indexList...\n")

    newseqs = {}
    numadded = 0
    while numadded < size:
        #i = rnd.sample(indexList, 1)
        #header = headers[i[0]]
        j = rnd.randint(0, len(indexList) -1)
        header = headers[ indexList[j] ]
        if header in newseqs:
            newseqs[header].updateCount(1)
        else:
            numadded += 1
            newseqs[header] = cp.copy( seqs[header] )
            newseqs[header].setCount(1)
    sys.stderr.write("\tDone filling in newseqs...\n")

    #Update sample with new sequences:
    #Set frequencies
    total = sum( [s.count for s in newseqs.values()] )
    for s in newseqs.values():
        s.setFreq(total)
    newsample.seqs = newseqs
    return newsample

def samplingSample(sample, size):
    sys.stderr.write("Begin to sample %d sequences from sample %s\n" %(size, sample.name))
    newsample = Sample(sample.name)
    seqs = sample.seqs
    sampleTotalReads = sum( [s.count for s in seqs.values()] )
    if sampleTotalReads < size:
        return

    headers = seqs.keys()
    indexList = []#index list represents all the sequences of the sample
    for i, h in enumerate(headers):
        s = seqs[h]
        indexList.extend( [i for j in xrange(s.count)] )
    sys.stderr.write("\tDone with the indexList...\n")
    chosenIndices = rnd.sample(indexList, size)
    sys.stderr.write("\tDone sampling from the indexList...\n")

    newseqs = {}
    for i in chosenIndices:
        header = headers[i]
        if header in newseqs:
            newseqs[header].updateCount(1)
        else:
            newseqs[header] = cp.copy( seqs[header] )
            newseqs[header].setCount(1)
    sys.stderr.write("\tDone filling in newseqs...\n")

    #Update sample with new sequences:
    #Set frequencies
    total = sum( [s.count for s in newseqs.values()] )
    for s in newseqs.values():
        s.setFreq(total)
    newsample.seqs = newseqs
    return newsample

def sampling(samples, size, uniq):
    '''Randomly pick 'size' sequences from each sample'''
    maxTotal = max( [ sum([seq.count for seq in sample.seqs.values()])  for sample in samples ] )

    newsamples = []
    for sample in samples:
        if uniq:
            newsample = samplingSample_uniq(sample, size)
            #newsample = samplingSample_weightedUniq(sample, size, maxTotal)
        else:
            newsample = samplingSample(sample, size)
        sys.stderr.write("Done sampling sample %s\n" %sample.name)
        newsamles.append(newsample)
    return newsamples
################# END SAMPLING ################

#============== UTILITIES =========
def readList(file, sep):
    l = []
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) > 0:
            l.extend( line.split(sep) )
    f.close()
    return l

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

def sortDictByValue( dictionary ):
    items = [ (v, k) for k,v in dictionary.items() ]
    items.sort()
    items.reverse()
    return [ (k, v) for v, k in items ]

def nt2aa(nt):
    #Translate nucleotide sequence to amino acid sequence
    codon2aa = {
    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
    'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
    'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
    'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
    'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
    'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
    'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
    'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S', 
    'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
    'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 
    'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', 
    'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
    'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
    'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 
    'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G' 
    }
    aaseq = ''
    for i in xrange(0, len(nt)/3):
        codon = nt[i*3 : i*3+3].upper()
        aaseq = aaseq + codon2aa[codon]
    return aaseq

def aaletters():
    aas = 'ACDEFGHIKLMNPQRSTVWY'
    #aas = ['F', 'S', 'Y', 'C', 'L', 'W', 'H', 'P', 'R', 'Q', 'I', 'T', 'N', 'M', 'V', 'A', 'K', 'D', 'E', 'G']
    return aas

def rc(nt):
    #Return a reverse complement of the input sequence.
    complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'a':'t', 't':'a', 'c':'g', 'g':'c'}
    rcnt = ''
    for i in xrange( len(nt) -1, -1, -1):
        if nt[i] in complement:
            rcnt += complement[nt[i]] 
        else:
            rcnt += nt[i]
    return rcnt

def getfiles(indir, ext):
    files = []
    for file in os.listdir(indir):
        items = file.split('.')
        if items[-1] == ext:
            files.append(file)
    return files

def getPc(count, total):
    if total == 0:
        return 0
    return 100.0*count/total

#============== LATEX =========
def prettyInt( number ):
    numstr = str( number )
    prettyStr = ''
    for i in range( len(numstr) ):
        if i > 0 and (len(numstr) - i) %3 == 0:
            prettyStr = prettyStr + ','
        prettyStr = prettyStr + numstr[i]

    return prettyStr

def prettyFloat( num ):
    if num == 0:
        return "0"
    else:
        return "%.2e" %num

def tabHeader(f, colnames):
    f.write("\\begin{table}\n")
    f.write("\\centering\n")
    f.write("\\scalebox{1.0}{%\n")
    f.write("\\begin{tabular}{c%s}\n" %("|c"*(len(colnames)-1)))
    f.write("\\hline\n")
    f.write("%s \\\\\n" %(" & ".join(colnames)))
    f.write("\\hline\n")

def writeDocumentStart( f ):
    f.write( "\\documentclass[11pt]{article}\n" )

    f.write( "\\usepackage{epsfig}\n" )
    f.write( "\\usepackage{multirow}\n" )
    f.write( "\\usepackage{graphicx}\n" )
    f.write( "\\usepackage{array}\n" )
    f.write( "\\usepackage{color, colortbl}\n" )
    f.write( "\\usepackage{rotating}\n" )
    f.write( "\\usepackage[table]{xcolor}\n" )
    f.write( "\\definecolor{LightCyan}{rgb}{0.88,1,1}")
    f.write( "\n" )

    f.write( "\\newcommand{\\figref}[1]{Figure~\\ref{fig:#1}}\n" )
    f.write( "\\newcommand{\\tabref}[1]{Table~\\ref{tab:#1}}\n" )
    f.write( "\n" )

    f.write("\\textwidth=6.5in\n")
    f.write("\\textheight=9in\n")
    f.write("\\oddsidemargin=0in\n")
    f.write("\\evensidemargin=0in\n")
    f.write("\\topmargin=0in\n")
    f.write("\\topskip=0in\n")
    f.write("\\headheight=0in\n")
    f.write("\\headsep=0in\n")
    f.write("\n")

    f.write("\\begin{document}\n")
    return

def writeDocumentEnd( f ):
    f.write( "\\end{document}\n" )

def tableCloser(f, captionStr, label):
    f.write("\\end{tabular}\n")
    f.write("}\n")
    f.write("\\caption{%s}\n" %captionStr)
    #f.write("\\end{center}\n")
    f.write("\\label{%s}" %label)
    f.write("\\end{table}\n\n")

def sidewaystableCloser(f, captionStr, label):
    f.write("\\end{tabular}\n")
    f.write("}\n")
    f.write("\\caption{%s}\n" %captionStr)
    #f.write("\\end{center}\n")
    f.write("\\label{%s}" %label)
    f.write("\\end{sidewaystable}\n\n")

#======= ERROR classes ===========
class FileFormatError( Exception ):
    pass

class InputOptionError( Exception ):
    pass

#======= initialize options =======
def initOptions():
    usage = "usage: %prog [options]\n"
    parser = OptionParser( usage=usage )
    return parser

def initOptions2():
    parser = initOptions()
    parser.add_option('-i', '--indir', dest='indir', help="Required argument: Input directory")
    parser.add_option('-o', '--outdir', dest='outdir', default='.', help="Output directory. Default = current directory")
    return parser

def initOptions3():
    parser = initOptions2()
    parser.add_option('--group2samples', dest='group2samples', help='Format: <groupname> <Comma separated list of samples of this group>')
    return parser

def checkOptions2(parser, args, options):
    if not options.indir:
        parser.error("Input directory was not specified\n")
    if not os.path.exists(options.indir):
        parser.error("Input directory %s does not exists\n" %options.indir)
    if not os.path.exists(options.outdir):
        parser.error("Output directory %s does not exists\n" %options.outdir)
        #system("mkdir -p %s" %options.outdir)

def checkOptions3(parser, args, options):
    checkOptions2(parser, args, options)
    if options.group2samples:
        options.group2samples, options.sample2group = readGroup2samples(options.group2samples)

#======= matplotlib functions ============
def initPlotOptions( parser ):
    parser.add_option( '--dpi', dest='dpi', default=300, type='int', help='Dots per inch. Default=%default')
    parser.add_option( '--outFormat', dest='outFormat', default='pdf', type='string',
                             help='Output plot format [pdf|png|all|eps]. default=%default' )

def checkPlotOptions( options, parser ):
    if options.dpi < 72:
        parser.error('dpi must >= than screen res, 72. Got (%d)' % options.dpi )
    if options.outFormat not in ('pdf', 'png', 'eps', 'all'):
        parser.error('Unrecognized output format: %s. Choose one from: pdf png eps all.' % options.outFormat )
    options.out = 'plot'

def initImage( width, height, options ):
    """
    initImage takes a width and height and returns
    both a fig and pdf object. options must contain outFormat,
    and dpi
    """
    import matplotlib.backends.backend_pdf as pltBack
    import matplotlib.pyplot as plt
    pdf = None
    if options.outFormat == 'pdf' or options.outFormat == 'all':
        pdf = pltBack.PdfPages( options.out + '.pdf' )
    fig = plt.figure( figsize=( width, height ), dpi=options.dpi, facecolor='w' )
    return ( fig, pdf )

def writeImage( fig, pdf, options ):
    """
    writeImage assumes options contains outFormat and dpi.
    """
    if options.outFormat == 'pdf':
        fig.savefig( pdf, format='pdf' )
        pdf.close()
    elif options.outFormat == 'png':
        fig.savefig( options.out + '.png', format='png', dpi=options.dpi )
    elif options.outFormat == 'all':
        fig.savefig( pdf, format='pdf' )
        pdf.close()
        fig.savefig( options.out + '.png', format='png', dpi=options.dpi )
        fig.savefig( options.out + '.eps', format='eps' )
    elif options.outFormat == 'eps':
        fig.savefig( options.out + '.eps', format='eps' )

def initImage2( width, height, outformat, outname, dpi ):
    """
    initImage takes a width and height and returns
    both a fig and pdf object. options must contain outFormat,
    and dpi
    """
    import matplotlib.backends.backend_pdf as pltBack
    import matplotlib.pyplot as plt
    pdf = None
    if outformat == 'pdf' or outformat == 'all':
        pdf = pltBack.PdfPages( outname + '.pdf' )
    fig = plt.figure( figsize=( width, height ), dpi=dpi, facecolor='w' )
    return ( fig, pdf )

def writeImage2( fig, pdf, outformat, out, dpi ):
    if outformat == 'pdf':
        fig.savefig( pdf, format='pdf' )
        pdf.close()
    elif outformat == 'png':
        fig.savefig( out + '.png', format='png', dpi=dpi )
    elif outformat == 'all':
        fig.savefig( pdf, format='pdf' )
        pdf.close()
        fig.savefig( out + '.png', format='png', dpi=dpi )
        fig.savefig( out + '.eps', format='eps' )
    elif outformat == 'eps':
        fig.savefig( out + '.eps', format='eps' )

def getColors0():
    colors = [ "#1f77b4", "#aec7e8", # blues 
               "#9467bd", "#c5b0d5", # lavenders
               "#ff7f0e", "#ffbb78", # oranges
               "#2ca02c", "#98df8a", # greens
               "#d62728", "#ff9896", # reds
               "#4A766E", "#48D1CC", #
               "#DAA520", "#DBDB70", #gold
               "#302B54", "#5959AB", #purple - blueish
               "#F6A8B6", "#F6CCDA"] #pink

    return colors

#AS color:
def getColors3():
    colors = [ "#FF3333", 
               "#00008B", 
               "#D02090", 
               "#FF6EB4", 
               "#EEAD0E", 
               "#CD2626", 
               "#FF7D40", 
               "#0000FF" ]
    return colors

def getColors6():
    #colors = ["#E31A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#1B9E77", "#FFFF33", "#A65628", "#CE1256"] #red, blue, green, purple, orange, green-ish, yellow, brown, pink 
    colors = ["#377EB8", "#E31A1C", "#4DAF4A", "#984EA3", "#FF7F00", "#1B9E77", "#FFFF33", "#A65628", "#CE1256"] #blue, red, green, purple, orange, green-ish, yellow, brown, pink 
    return colors

def getColors6light():
    colors = ["#A6D7FE", "#FE8E8F", "#B8FEB5", "#F6BDFE", "#FEBF80", "#95FEDF", "#FFFFB3", "#D8885A", "#D7B5D8"] #blue, red, green, purple, orange, green-ish
    return colors

def getColors6dark():
    colors = ["#275880", "#B30000", "#367A33", "#6A3772", "#B25900", "#136E53", "#B2B324", "#743C1C", "#900D3D"]
    return colors

########## HACK COLOR ############
def sample2color( colors ):
    s2ci = {"SBC1": 0, "SBC2": 0, "SBC3": 0, "SBC4": 0, "SBC5": 0, "SBC6":0,
           "SBC7": 1, "SBC8": 1,
           "Patient-01-R":2, "Patient-B-R":2, "Patient-10":2, "Patient-11":2, "Patient-12":2, "Patient-13":2,
           "Patient-01-D":3, "Patient-B-D":3, "Patient-8-D":3, "Patient-15-D":3,
           "SEQ_A1a":4, "SEQ_A1b":4, "SEQ_A1":4, "SEQ_A2":4, "SEQ_B":4,
           "PtX-1-Month-Post":5, "PtX-2-Months-Post":5, "PtX-Before-transplant":5,
           "female":4,
           "male2":5,
           "male1_CD45RA+RO-_naive_day1":7, "male1_CD45RA+RO-_naive_day8":7, "male1_CD45RO+RA-_memory_day1":7, "male1_CD45RO+RA-_memory_day8":7, "male1_blooddraw1":7, "male1_blooddraw2":7,
           "iRepertoire-RNA-AS": 0,
           "iRepertoire-DNA-AS": 1,
           "iRepertoire-RNA-Control": 0,
           "iRepertoire-DNA-Control": 1,
           "adaptiveTCR-RNA-AS": 2,
           "adaptiveTCR-DNA-AS": 3,
           "adaptiveTCR-RNA-Control": 2,
           "adaptiveTCR-DNA-Control": 3,
           "adaptiveTCR-healthy": 4,
           "adaptiveTCR-PtX": 5,
           "warren": 6,
           "adaptF28cd8":1, "adaptF57cd8":1, "adaptM35cd8":1, "asBDR":1,
           "adaptF28cd4":2, "adaptF57cd4":2, "adaptM35cd4":2,
           #"as8D":0, "as10R":0, "as11R":0, "as12R":0, "as13R":0, "as15D":0, "as1DR":0,
           "as10R":0, "as11R":0, "as12R":0, "as13R":0, "as1R":0,
           "adaptCD8R":1, "adaptMA8R":1,
           #"as8D":1, "as15D":1, "as1D":1, "as16D":1, "as11D":1, "as1Ddraw2":1,
           "as8D":0, "as15D":2, "as1D":0, "as16D":0, "as11D":0, "as1Ddraw2":0,
           "as1Ddraw2notcd8": 2,
           #"asBD":4, "asBDdraw2":4, "as20D":4,
           "asBD":1, "asBDdraw2":1, "as20D":1,
           "asBR":1,
           "diffGroupShared":0, "sameGroupShared":1, "shared":2, "uniq":3,
           'cd4':0, 'cd8':1, 'f57cd4':0, 'f57cd8':1, 'm35cd4':0, 'm35cd8':1,
           'adapt1D':1, 'adaptBD':1, 'adapt11D':1,
           'adapt1R':2, 'adaptBR':2, 'adapt11R':2,
           'rna':0, 'dna':1,
           'adaptControls':2, 'asControls':4, 'asPatients':1,
           "adaptF28":2, "adaptF57":2, "adaptM35":2,
           'adaptAS':2, 'adaptCD':2, 'adaptMA':2,'adaptBSD1n':2,
           #"female":8,
           #"male2":6,
           #"male1_CD45RA+RO-_naive_day1":7, "male1_CD45RA+RO-_naive_day8":7, "male1_CD45RO+RA-_memory_day1":7, "male1_CD45RO+RA-_memory_day8":7, "male1_blooddraw1":7, "male1_blooddraw2":7
          }
     
    #s2ci = {"SBC1": 20, "SBC2": 20, "SBC3": 20, "SBC4": 20, "SBC5": 20, "SBC6":20, #RNA
    #       "SBC7": 22, "SBC8": 22, #DNA
    #       "Patient-01-R":20, "Patient-B-R":20, "Patient-10":20, "Patient-11":20, "Patient-12":20, "Patient-13":20, #RNA
    #       "Patient-01-D":22, "Patient-B-D":22, "Patient-8-D":22, "Patient-15-D":22, #DNA
    #       "SEQ_A1a":4, "SEQ_A1b":4, "SEQ_A1":4, "SEQ_A2":4, "SEQ_B":4, #
    #       "PtX-1-Month-Post":6, "PtX-2-Months-Post":6, "PtX-Before-transplant":6 #
    #      }
    
    markers = ['o', 'd', '^', 'p', 'v', '*', 's']
    s2mi = {"SBC1":1 , "SBC2": 2, "SBC3": 3, "SBC4": 4, "SBC5": 5, "SBC6":6, #RNA
           "SBC7": 1, "SBC8": 2, #DNA
           "Patient-01-R":1, "Patient-B-R":2, "Patient-10":3, "Patient-11":4, "Patient-12":5, "Patient-13":6, #RNA
           "Patient-01-D":1, "Patient-B-D":2, "Patient-8-D":3, "Patient-15-D":4, #DNA
           "SEQ_A1a":0, "SEQ_A1b":1, "SEQ_A1":2, "SEQ_A2":3, "SEQ_B":4, #
           "PtX-1-Month-Post":0, "PtX-2-Months-Post":1, "PtX-Before-transplant":2, #
           "female":2,
           "male2":2,
           "male1_CD45RA+RO-_naive_day1":0, "male1_CD45RA+RO-_naive_day8":1, "male1_CD45RO+RA-_memory_day1":2, "male1_CD45RO+RA-_memory_day8":3, "male1_blooddraw1":4, "male1_blooddraw2":5,
           "iRepertoire-RNA-AS": 1,
           "iRepertoire-DNA-AS": 1,
           "iRepertoire-RNA-Control": 2,
           "iRepertoire-DNA-Control": 2,
           "adaptiveTCR-RNA-AS": 1,
           "adaptiveTCR-DNA-AS": 1,
           "adaptiveTCR-RNA-Control": 2,
           "adaptiveTCR-DNA-Control": 2,
           "adaptiveTCR-healthy": 3,
           "adaptiveTCR-PtX": 4,
           "warren": 5,
           "adaptF28cd8":1, "adaptF57cd8":2, "adaptM35cd8":3, "asBDR":0,
           "adaptF28cd4":1, "adaptF57cd4":2, "adaptM35cd4":3,
           #"as8D":1, "as10R":2, "as11R":3, "as12R":4, "as13R":5, "as15D":6, "as1DR":0,
           "as10R":1, "as11R":2, "as12R":3, "as13R":4, "as1R":5,
           "adaptCD8R":1, "adaptMA8R":2,
           "as8D":1, "as15D":3, "as1D":5, "as16D":4, "as11D":2, "as1Ddraw2":6,
           "as1Ddraw2notcd8": 6,
           "asBD":5, "asBDdraw2":6, "as20D":3,
           "asBR":5,
           "diffGroupShared":0, "sameGroupShared":1, "shared":2, "uniq":3,
           'cd4':0, 'cd8':1, 'f57cd4':0, 'f57cd8':0, 'm35cd4':1, 'm35cd8':1,
           'adapt1D':0, 'adaptBD':1, 'adapt11D':2,
           'adapt1R':0, 'adaptBR':1, 'adapt11R':2,
           'rna':0, 'dna':4,
           'adaptControls':0, 'asControls':1, 'asPatients':2,
           "adaptF28":0, "adaptF57":1, "adaptM35":2,
           'adaptAS':3, 'adaptCD':4, 'adaptMA':5,'adaptBSD1n':6,
          }
    
    lightcolors = getColors6light()
    s2c = {}
    s2m = {}
    s2cLight = {}
    for s, i in s2ci.iteritems():
        s2c[ s ] = colors[i]
        s2cLight[ s ] = lightcolors[i]
    for s, i in s2mi.iteritems():
        s2m[ s ] = markers[i]
    #return s2c, s2m, s2cLight
    return s2c, s2m, s2c

####Proper name for paper ###
def properName(name):
    mydict = {'as1D':'AS1', 'as1Ddraw2': 'AS1',
              'as11D': 'AS2', 'as16D':'AS3', 'as8D':'AS4', 'as15D':'AS5',
              #'asBD': 'H1_draw1', 'asBDdraw2':'H1', 'as20D':'H2', 
              'asBD': 'H1', 'asBDdraw2':'H1', 'as20D':'H2', 
              'adaptF57':'EH1', 'adaptM35':'EH2',
              'adapt1D':'AS1', 'adapt1Ddraw2': 'AS1',
              'adapt11D': 'AS2', 'adapt16D':'AS3', 'adapt8D':'AS4', 'adapt15D':'AS5',
              'adaptBD': 'H1', 'adaptBDdraw2':'H1', 'adapt20D':'H2', 
              'controls':'B27+,Healthy', 'patients':'B27+,AS',
              }
    if name not in mydict:
        return name
    else:
        return mydict[name]

def properName2name(propername):
    #mydict = {'as1D':'AS1', 'as1Ddraw2': 'AS1', 'as11D': 'AS2', 'as16D':'AS3', 'as8D':'AS4', 'as15D':'AS5', 'asBD': 'H1', 'asBDdraw2':'H1', 'as20D':'H2', 'adaptF57':'EH1', 'adaptM35':'EH2'}
    #rvdict = {'AS3': 'as16D', 'AS2': 'as11D', 'AS1': 'as1Ddraw2', 'EH2': 'adaptM35', 'AS5': 'as15D', 'AS4': 'as8D', 'H2': 'as20D', 'H1': 'asBDdraw2', 'EH1': 'adaptF57'}
    rvdict = {'AS3': 'adapt16D', 'AS2': 'adapt11D', 'AS1': 'adapt1Ddraw2', 'EH2': 'adaptM35', 'AS5': 'adapt15D', 'AS4': 'adapt8D', 'H2': 'adapt20D', 'H1': 'adaptBDdraw2', 'EH1': 'adaptF57'}
    #rvdict = {'AS3': 'as16D', 'AS2': 'as11D', 'AS1': 'as1D', 'EH2': 'adaptM35', 'AS5': 'as15D', 'AS4': 'as8D', 'H2': 'as20D', 'H1': 'asBD', 'EH1': 'adaptF57'}
    if propername not in rvdict:
        return propername
    else:
        return rvdict[propername]

########### END HACK COLOR #########



def setAxes( fig ):                                                                     
    #Set single axes
    return fig.add_axes( [0.12, 0.1, 0.83, 0.85] )                                      

def setAxes2( fig, numSamples, samplesPerPlot ):
    #Set multiple axes, depending on the number of samplesPerPlot and total Samples
    numaxes = numSamples / samplesPerPlot 
    if numSamples % samplesPerPlot != 0:
        numaxes += 1
    
    axesList = []
    axleft = 0.12
    axright = 0.97
    axwidth = axright - axleft
    axbottom = 0.15
    axtop = 0.93
    axheight = axtop - axbottom
    margin = 0.1

    h = ( axheight - ( margin * (numaxes - 1) ) )/float( numaxes )

    bottom = axtop - h
    for i in range( numaxes ):#add axes from top down
        axesList.append( fig.add_axes( [axleft, bottom, axwidth, h] ) )
        bottom = bottom - (margin + h)
    
    return axesList

def setAxes3(fig, range1, range2):
    axleft = 0.12
    axright = 0.95
    axwidth = axright - axleft
    axbottom = 0.15
    axtop = 0.95
    axheight = axtop - axbottom
    margin = 0.05

    h1 = (axheight - margin)*(range1/(range1 + range2))
    h2 = axheight - margin - h1

    ax2 = fig.add_axes( [axleft, axbottom, axwidth, h2] )
    ax = fig.add_axes( [axleft, axbottom + h2 + margin, axwidth, h1] )
    return ax, ax2

def drawDiscontinueSign(topAxes, bottomAxes, top, bottom):
    d = 2
    topAxes.plot( (-0.6, -0.4), (top + d, top - d), color="k", clip_on=False )
    bottomAxes.plot( (-0.6,-0.4), (bottom +d, bottom -d), color="k", clip_on=False)

def editSpine2(topAxes, bottomAxes):
    topAxes.spines['bottom'].set_visible(False)
    topAxes.spines['top'].set_visible(False)
    topAxes.spines['right'].set_visible(False)
    topAxes.yaxis.set_ticks_position('left')
    topAxes.xaxis.set_ticks_position('none')

    bottomAxes.spines['top'].set_visible(False)
    bottomAxes.spines['right'].set_visible(False)
    bottomAxes.xaxis.tick_bottom()
    bottomAxes.yaxis.set_ticks_position( 'left' )

def editSpine( axes ):                                                                 
    for loc, spine in axes.spines.iteritems():                                          
        if loc in [ 'left', 'bottom' ]:                                                 
            spine.set_position( ('outward', 10) )                                       
        elif loc in ['right', 'top']:                                                   
            spine.set_color('none')                                                     
        else:                                                                           
            raise ValueError( 'Unknown spine location %s\n' % loc )                     
                                                                                        
def setTicks( axes ):                                                                   
    axes.xaxis.set_ticks_position( 'bottom' )                                           
    axes.yaxis.set_ticks_position( 'left' )                                             
    minorLocator = LogLocator( base=10, subs = range(1, 10) )                           
    axes.xaxis.set_minor_locator( minorLocator )                                        
                                                
def bihist( y1, y2, axes, bins, orientation, color=None ):
    n1, bins1, patch1 = axes.hist( y1, bins=bins, orientation=orientation, color=color ) #Top hist
    n2, bins2, patch2 = axes.hist( y2, bins=bins, orientation=orientation, color=color ) #bottom hist
    #set ymax:
    if orientation == 'vertical':
        ymax = max( [ i.get_height() for i in patch1 ] )
        for i in patch2:
            i.set_height( -i.get_height() )
        ymin = min( [ i.get_height() for i in patch2 ] )
    elif orientation == 'horizontal':
        ymax = max( [ i.get_width() for i in patch1 ] )
        for i in patch2:
            i.set_width( -i.get_width() )
        ymin = min( [ i.get_width() for i in patch2 ] )
    
    #axes.set_ylim( ymin*1.1, ymax*1.1 )
    return ymin, ymax

