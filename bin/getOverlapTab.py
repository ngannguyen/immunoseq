#!/usr/bin/env python2.6

'''
nknguyen soe ucsc edu
Feb 14 2012
Aggregate overlap statistics from different experiments into a table
Input: input directory containing overlap-statistic files, one for each pair of samples
Output: 1/Mode1: Table of overlapping statistics across different cutoffs where: Rows = pairs of samples. Cols = "Sample1\tSample2\tSamplingSize\tCutoff 1\t Cutoff 2\t ...
                 Plot: xaxis: cutoffs (log scale). yaxis: Percentage of reads 
        2/Mode2: For each cutoff: table of different overlapping statistics ( number of clones, number of overlapped clones, percentage of overlapped reads )
                 Cols: Sample1\tSample2\tClones1\tClones2\tOclones1\tOclones2\t%clones1o2\t%clones2o1\t%reads1o2\t%reads2o1
                 Rows: Pairs of samples
'''

import os, sys, re
import immunoseq.lib.immunoseqLib as iseqlib

class Exp():
    def __init__(self, name):
        self.name = name
        nameitems = name.split('_')
        if len(nameitems) < 3:
            #raise ValueError("Wrong filename format. Required: experimentName-samplingsize\n")
            self.size = 0
        else:
            self.size = int(nameitems[1].rstrip("k"))
        
        self.exp = '_'.join( nameitems[:2] )
        if len(nameitems) < 2: 
            print nameitems
            raise ValueError("Wrong filename format. Required: sample1_sample2[-samplingsize] where '-samplingsize' is optional. Filename: %s\n" %name)
        self.sample1 = nameitems[0]
        self.sample2 = nameitems[1]

        self.cutoffs = []
        self.clones1 = [] #total number of clones in sample 1 passed the cutoffs
        self.clones2 = [] #total number of clones in sample 2 passed the cutoffs
        
        self.oclones1 = [] #percentage of clones in sample 1 that passed the cutoffs and are in sample 2
        self.oclones2 = [] #percentage of clones in sample 2 that passed cutoffs and present in sample 1
        self.avrOclones = []

        self.stdOclones1 = [] 
        self.stdOclones2 = [] 
        self.stdAvrOclones = [] 
        
        self.reads1o2 = [] #Percentage of sample 1 reads in the overlapped clones
        self.reads2o1 = [] #Percentage of sample 2 reads in the overlapped clones
        self.avrOreads = [] #average of reads1o2 and reads2o1
        
        self.stdReads1o2 = [] #Standard deviation of Percentage of sample 1 reads in the overlapped clones
        self.stdReads2o1 = [] #Standard deviation of Percentage of sample 2 reads in the overlapped clones
        self.stdAvrOreads = [] #Standard deviation average of reads1o2 and reads2o1

    def addCutoffStats(self, items):
        self.cutoffs.append( float(items[0]) )
        self.clones1.append( float(items[1]) )
        self.clones2.append( float(items[2]) )

        self.oclones1.append( float(items[5]) )
        self.oclones2.append( float(items[6]) )
        
        self.reads1o2.append( float(items[7]) )
        self.reads2o1.append( float(items[8]) )
        self.avrOreads.append( (float(items[7]) + float(items[8]))/2.0 )
        
        if len(items) == 17:
            self.stdOclones1.append( float(items[13]) )
            self.stdOclones2.append( float(items[14]) )
            self.stdAvrOclones.append( (float(items[13]) + float(items[14]))/2.0 )
            
            self.stdReads1o2.append( float(items[15]) )
            self.stdReads2o1.append( float(items[16]) )
            self.stdAvrOreads.append( (float(items[15]) + float(items[16]))/2.0 )

class FileFormatError(Exception):
    pass

class InconsistentCutoffsError(Exception):
    pass

def sortByOrder(exps, sampleOrder):
    sortedexps = []
    for s in sampleOrder:
        for e in exps:
            if e.name == s:
                sortedexps.append(e)
                break
    return sortedexps

#def groupAvr(group2samples, sample, exps):
#    group2avr = {} #key = group, val = Experiment object with average stats
#    for group, samples in group2samples.iteritems():
#        expname = "%s_%s" %(sample, group)
#        exp = Experiment(group)
#        for s in samples:
#            for e in exps:
#                if s == e.sample1:
#                   for c in e.cutoffs:
#                        if c not in exp.cutoffs:
#                            exp.cutoffs.append(c)
#   

#====== FIG ======
def drawOverlapReadsData(axes, exps, sample, sampleOrder):
    if len( exps ) <= 0:
        return
    if sampleOrder:
        exps = sortByOrder(exps, sampleOrder)
    else:
        exps = sorted( exps, key = lambda exp: (exp.exp, exp.size) )

    colors = iseqlib.getColors6()
    lightcolors = iseqlib.getColors6light()
    markers = ['o', '*', 's', 'd', '^', 'p', 'v']
    textsize = 'large'
    lines = []
    expnames = []
   
    cutoffs = exps[0].cutoffs
    xticklabels = [ str(c) for c in cutoffs ]
    xdata = cutoffs
    for i, x in enumerate(xdata):
        if x == 0 and i < len(xdata) - 1:
            xdata[i] = xdata[i + 1]/10.0

    prevexp = ''
    colorIndex = 0
    markerIndex = -1
    axes.set_xscale('log')
    #axes.set_yscale('log')

    maxy = 0
    miny = 'inf'
    for exp in exps:
        expnames.append( exp.name )
        ydata, stddata = getStats(exp, sample)

        maxy = max( [max(ydata), maxy] )
        miny = min( [min(ydata), miny] )
        if prevexp != '' and prevexp != exp.exp:
            colorIndex += 1
            markerIndex = 0
        else:
            markerIndex += 1
        prevexp = exp.exp
        markersize = 12.0
        if markers[markerIndex] == '*':
            markersize += 2.0
        elif markers[markerIndex] == 's':
            markersize -= 2.0

        l = axes.plot(xdata, ydata, color=colors[colorIndex], marker=markers[markerIndex], markeredgecolor=colors[colorIndex], markersize=markersize, linestyle='none')
        axes.plot( xdata, ydata, color=lightcolors[colorIndex], linestyle='-', linewidth=0.2 )
        lines.append(l)
        
    axes.set_title('Percentage of Reads in Overlapped Clones', size='xx-large') 
    axes.set_xlabel('Minimum clone size (percentage of total reads)', size=textsize)
    axes.set_ylabel('Percentage of reads overlapped', size=textsize)
    axes.xaxis.set_ticklabels( xticklabels )
    axes.xaxis.set_ticks( xdata )
    minx = min(xdata)
    maxx = max(xdata)
    axes.set_xlim( 7*minx/10, maxx*1.5 )

    #axes.set_ylim(-2, 102 )
    axes.set_ylim(-0.5, 8 )
    axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    legend = axes.legend( lines, expnames , numpoints=1, loc='best', ncol=1)
    legend._drawFrame = False

def drawOverlapReads(exps, options):
    sample = options.sample
    options.out = os.path.join(options.outdir, "overlapPlot-%s" % sample)
    fig, pdf = iseqlib.initImage( 10.0, 12.0, options )
    axes = iseqlib.setAxes(fig)
    drawOverlapReadsData(axes, exps, sample, options.sampleOrder)
    iseqlib.writeImage(fig, pdf, options)

#====== LATEX =====
def tabHeader(f, cutoffs):
    f.write("\\begin{table}\n") 
    f.write("\\centering\n")
    f.write("\\scalebox{0.7}{%\n")
    f.write("\\begin{tabular}{r|r|r|%s}\n" %( '|'.join(['r' for c in cutoffs]) ) )
    f.write("\\hline\n")
    f.write( "S1 & S2 & SamplingSize & %s\\\\\n" %( '&'.join( ["%.3f" %c for c in cutoffs] ) ) )
    f.write("\\hline\n")

def getStats(exp, type):
    #clone1, clone2, cloneAvr, read1, read2, readAvr
    type2data = {'clone1':(exp.oclones1, exp.stdOclones1),
                 'clone2':(exp.oclones2, exp.stdOclones2),
                 'cloneAvr':(exp.avrOclones, exp.stdAvrOclones),
                 'read1':(exp.reads1o2, exp.stdReads1o2),
                 'read2':(exp.reads2o1, exp.stdReads2o1),
                 'readAvr':(exp.avrOreads, exp.stdAvrOreads)}
    return type2data[type][0], type2data[type][1]

def tab(f, exps, statsType, sampleOrder):
    if sampleOrder:
        exps = sortByOrder(exps, sampleOrder)
    else:
        exps = sorted( exps, key = lambda exp: (exp.exp, exp.size) )
    prevexp = ''
    for exp in exps:
        #nameitems = exp.name.split('-')
        #if len(nameitems) != 2:
        #    raise ValueError("Wrong filename format. Required: experimentName-samplingsize\n")
        #f.write( "%s & %s" %(nameitems[0].replace("_", " "), nameitems[1]) )
        if prevexp != '' and exp.exp != prevexp:
            f.write("\\hline\n")
        prevexp = exp.exp
        if exp.size == 0:
            f.write( "%s & %s & NA" %(exp.sample1, exp.sample2) )
        else:
            f.write( "%s & %s & %dk" %(exp.sample1, exp.sample2, exp.size) )
        
        data, stddata = getStats(exp, statsType)
        print data
        print stddata
        for i, p in enumerate( data ):
            if len( stddata ) > i:
                f.write( "& %.2f $\pm$ %.2f" % (p, stddata[i]) )
            else:
                f.write( "& %.2f" % p )
        f.write("\\\\\n")
    f.write("\\hline\n")

#def getOverlapReadsLatexTab(exps, file, sample, sampleOrder):
def getCutoffsLatexTab(exps, file, statsType, sampleOrder):
    if len(exps) == 0:
        return
    
    cutoffs = exps[0].cutoffs
    for exp in exps:
        for i, c in enumerate(exp.cutoffs):
            if c != cutoffs[i]:
                raise InconsistentCutoffsError("Input files don't have the same cutoffs.")
    
    f = open(file, 'w')
    iseqlib.writeDocumentStart(f)
    #Table:
    cutoffs = exps[0].cutoffs
    tabHeader(f, cutoffs)
    tab(f, exps, statsType, sampleOrder)

    label = ''
    captionStr = ''
    iseqlib.tableCloser(f, captionStr, label)
    iseqlib.writeDocumentEnd(f)

    f.close()
#====== END LATEX =====

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
        f.write( "%s\t%s" %(nameitems[0], nameitems[1]) )
        for i, p in enumerate( exp.avrOreads ):
            if len( exp.stdAvrOreads ) > i:
                f.write( "\t%.2f (%.2f)" % (p, exp.stdAvrOreads[i]) )
            else:
                f.write( "\t%.2f" % p )
        f.write("\n")

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
        if float(items[0]) >= 5.0:
            continue
        exp.addCutoffStats(items)
    f.close()
    return exp

#============= MODE 2: All statistics per cutoff =================
def tabHeader2(f):
    f.write("\\begin{table}\n")
    f.write("\\centering\n")
    f.write("\\scalebox{0.7}{%\n")
    f.write("\\begin{tabular}{r|r|r|r|r|r|r|r|r}\n")
    f.write("\\hline\n")
    #f.write("S1 & S2 & Size & Clones1 & Clones2 & \%O.Clones1 & \%O.Clones2 & \%O.Reads1 & \%O.Reads 2\\\\\n")
    f.write("S1 & S2 & Size & Clones1 & \%O.Clones1 & \%O.Reads1 & Clones2 & \%O.Clones2 & \%O.Reads2\\\\\n")
    f.write("\\hline\n")

def tab2(f, sampleToExps, cutoff, sampleOrder):
    samples = sampleOrder
    if not sampleOrder:
        samples = sorted( sampleToExps.keys() )

    for sample in samples:
        exps = sampleToExps[sample]
        for exp in exps:
            #Find index of the cutoff
            index = -1
            for i,c in enumerate(exp.cutoffs):
                if c == cutoff:
                    index = i
                    break
            if index == -1:
                continue

            samsize = str(exp.size)
            if samsize == '0':
                samsize = 'NA'
            if exp.sample1 == sample:
                #f.write("%s & %s & %s & %s & %s & %.2f & %.2f & %.2f & %.2f \\\\\n" %(exp.sample1, exp.sample2, samsize, iseqlib.prettyInt(int(exp.clones1[index])), iseqlib.prettyInt(int(exp.clones2[index])), exp.oclones1[index], exp.oclones2[index], exp.reads1o2[index], exp.reads2o1[index]) )
                f.write("%s & %s & %s & %s & %.2f & %.2f & %s & %.2f & %.2f \\\\\n" %(exp.sample1, exp.sample2, samsize, iseqlib.prettyInt(int(exp.clones1[index])), exp.oclones1[index], exp.reads1o2[index], iseqlib.prettyInt(int(exp.clones2[index])), exp.oclones2[index], exp.reads2o1[index]) )
            else:
                #f.write("%s & %s & %s & %s & %s & %.2f & %.2f & %.2f & %.2f \\\\\n" %(exp.sample2, exp.sample1, samsize, iseqlib.prettyInt(int(exp.clones2[index])), iseqlib.prettyInt(int(exp.clones1[index])), exp.oclones2[index], exp.oclones1[index], exp.reads2o1[index], exp.reads1o2[index]) )
                f.write("%s & %s & %s & %s & %.2f & %.2f & %s & %.2f & %.2f \n" %(exp.sample2, exp.sample1, samsize, iseqlib.prettyInt(int(exp.clones2[index])), exp.oclones2[index], exp.reads2o1[index], iseqlib.prettyInt(int(exp.clones1[index])), exp.oclones1[index], exp.reads1o2[index]) )
        f.write("\\hline\n")        

def getAllStatsLatexTab(sampleToExps, cutoff, outfile, sampleOrder):
    f = open(outfile, 'w') 
    iseqlib.writeDocumentStart(f)
    tabHeader2(f)
    tab2(f, sampleToExps, cutoff, sampleOrder)
    label = ''
    captionStr = ''
    iseqlib.tableCloser(f, captionStr, label)
    iseqlib.writeDocumentEnd(f)
    f.close()

#def textTab(f, sampleToExps, cutoff, sampleOrder, group2samples):
def textTab(f, sampleToExps, cutoff, sampleOrder):
    #f.write("#S1\tS2\tSize\tClones1\tClones2\t%O.Clones1\t%O.Clones2\t%O.Reads1\t%O.Reads2\n")
    f.write("#S1\tS2\tSize\tClones1\t%O.Clones1\t%O.Reads1\tClones2\t%O.Clones2\t%O.Reads2\n")
    samples = sampleOrder
    if not sampleOrder:
        samples = sorted( sampleToExps.keys() )

    for sample in samples:
        exps = sampleToExps[sample]
        for exp in exps:
            #Find index of the cutoff
            index = -1
            for i,c in enumerate(exp.cutoffs):
                if c == cutoff:
                    index = i
                    break
            if index == -1:
                continue

            samsize = str(exp.size)
            if samsize == '0':
                samsize = 'NA'
            if exp.sample1 == sample:
                f.write("%s\t%s\t%s\t%s\t%.2f\t%.2f\t%s\t%.2f\t%.2f \n" %(exp.sample1, exp.sample2, samsize, iseqlib.prettyInt(int(exp.clones1[index])), exp.oclones1[index], exp.reads1o2[index], iseqlib.prettyInt(int(exp.clones2[index])), exp.oclones2[index], exp.reads2o1[index]) )
                #f.write("%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f \n" %(exp.sample1, exp.sample2, samsize, iseqlib.prettyInt(int(exp.clones1[index])), iseqlib.prettyInt(int(exp.clones2[index])), exp.oclones1[index], exp.oclones2[index], exp.reads1o2[index], exp.reads2o1[index]) )
            else:
                f.write("%s\t%s\t%s\t%s\t%.2f\t%.2f\t%s\t%.2f\t%.2f \n" %(exp.sample2, exp.sample1, samsize, iseqlib.prettyInt(int(exp.clones2[index])), exp.oclones2[index], exp.reads2o1[index], iseqlib.prettyInt(int(exp.clones1[index])), exp.oclones1[index], exp.reads1o2[index]) )
                #f.write("%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f \n" %(exp.sample2, exp.sample1, samsize, iseqlib.prettyInt(int(exp.clones2[index])), iseqlib.prettyInt(int(exp.clones1[index])), exp.oclones2[index], exp.oclones1[index], exp.reads2o1[index], exp.reads1o2[index]) )
        f.write("#\n")        
    
#def getAllStatsLatexTabs(exps, outdir, cutoffs, sampleOrder, latex, group2samples):
def getAllStatsLatexTabs(exps, outdir, cutoffs, sampleOrder, latex):
    sampleToExps = {} #key = sampleName, vals = list of exps
    allFlag = False
    if 'all' in cutoffs:
        cutoffs = []
        allFlag = True
    for exp in exps:
        if exp.sample1 in sampleToExps:
            sampleToExps[ exp.sample1 ].append( (exp, exp.sample2) )
        else:
            sampleToExps[ exp.sample1 ] = [ (exp, exp.sample2) ]
        
        if exp.sample2 in sampleToExps:
            sampleToExps[ exp.sample2 ].append( (exp, exp.sample1) )
        else:
            sampleToExps[ exp.sample2 ] = [ (exp, exp.sample1) ]
        
        #Cutoffs:
        for c in exp.cutoffs:
            if allFlag and c not in cutoffs:
                cutoffs.append(c)

    for s in sampleToExps:
        sampleToExps[s] = sorted( sampleToExps[s], key=lambda e:e[1] )
        sampleToExps[s] = [ e[0] for e in sampleToExps[s] ]

    for c in cutoffs:
        if not isinstance(c, float):
            c = float(c)
        if latex:
            outfile = os.path.join(outdir, "overlap-%.3f.tex" %c)
            getAllStatsLatexTab(sampleToExps, c, outfile, sampleOrder)
        else:
            outfile = os.path.join(outdir, "overlap-%.3f.txt" %c)
            f = open(outfile, 'w')
            #textTab(f, sampleToExps, c, sampleOrder, group2samples)
            textTab(f, sampleToExps, c, sampleOrder)
            f.close()

#def checkOptions(options):
def readGroup2Samples(file): 
    group2samples = {}
    f = open(file, 'r')
    for line in f:
        items = line.strip().split()
        if len(items) < 2:
            continue
        group = items[0]
        sample = items[1]
        if group not in group2samples:
            group2samples[group] = [sample]
        else:
            group2samples[group].append(sample)
    f.close()
    return group2samples

def main():
    parser = iseqlib.initOptions()
    iseqlib.initPlotOptions(parser)
    parser.add_option('-i', '--indir', dest='indir')
    parser.add_option('-o', '--outdir', dest='outdir', default = '.')
    
    parser.add_option('-m', '--mode', dest='mode', default='1,2', type='string', help='Specify how you would like to aggregate the data. Should be a comma seperated list of any of the valid choices: [1, 2]. Mode 1 is one statistics of interest across different cutoffs. Mode 2 is different statistics for 1 cutoff.')

    parser.add_option('-s', '--stats_type', dest='stats_type', default='readAvr', help='Option for mode 1. Specify which overlapping statistics of interest to print out. Default = %default. Valid values are [clone1, clone2, cloneAvr, read1, read2, readAvr], where clone# is percentage of clones in sample # that passed each cutoff and are also in the other sample; cloneAvr is average of clone1 and clone2; read# is percentage of reads in clones of sample # that passed cutoff and also in the other sample; readAvr: average of read1 and read2')
    parser.add_option('-c', '--cutoffs', dest='cutoffs', default='all', help='Option for mode 2. Comma separated list of cutoffs of interest. Default=%default' )
    parser.add_option('-a', '--sampleOrder', dest='sampleOrder')
    parser.add_option('-l', '--latex', dest='latex', action='store_true', default=False)
    parser.add_option('-g', '--groupToSamples', dest='group2samples')

    options, args = parser.parse_args()
    iseqlib.checkPlotOptions( options, parser )

    if options.sampleOrder:
        options.sampleOrder = options.sampleOrder.split(',')
    indir = options.indir
    if not os.path.isdir(indir):
        raise ValueError("Input directory %s is not a directory\n" %indir)
    options.mode = options.mode.split(',')
    options.cutoffs = options.cutoffs.split(',')

    exps = []
    for file in os.listdir(indir):
        if os.path.isdir( os.path.join(indir, file) ):
            continue
        exp = readFile( os.path.join(indir, file) )
        exps.append(exp)
    
    if options.group2samples:
        options.group2samples = readGroup2Samples(options.group2samples)

    orfile = os.path.join(options.outdir, "overlapReads-%s.tex" %options.stats_type)
    #getOverlapReadsTab(exps, orfile)
    if '1' in options.mode:
        getCutoffsLatexTab(exps, orfile, options.stats_type, options.sampleOrder)
        #drawOverlapReads(exps, options)
    
    if '2' in options.mode:
        getAllStatsLatexTabs(exps, options.outdir, options.cutoffs, options.sampleOrder, options.latex)
        #getAllStatsLatexTabs(exps, options.outdir, options.cutoffs, options.sampleOrder, options.latex, options.group2samples)

if __name__ == '__main__':
    main()
