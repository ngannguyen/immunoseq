#!/usr/bin/env python2.6

#nknguyen at soe ucsc edu
#May 30 2012
#Draw scatter plot of each clones and their size for each pair of samples
#xaxis: Sample 1, yaxis: sample 2
#each dot represents a clonotype with same amino acid, and V and J genes
#

import os, re, sys, copy
from optparse import OptionParser

from immunoseq.lib.immunoseqLib import *
import matplotlib.cm as cm 
import matplotlib.pyplot as pyplot
from matplotlib.ticker import *
from matplotlib.font_manager import FontProperties
import matplotlib.backends.backend_pdf as pltBack

def getData(seqs1, seqs2):
    points = {} #key = (xval, yval) #val = counts (number of points with this x and y)
    for aa, v2j in seqs1.iteritems():
        for v, j2seq in v2j.iteritems():
            for j, seq in j2seq.iteritems():
                x = seq.count
                y = 0
                if aa in seqs2 and v in seqs2[aa] and j in seqs2[aa][v]:
                    y = seqs2[aa][v][j].count
                point = (x, y)
                if point not in points:
                    points[point] = 1
                else:
                    points[point] += 1

    for aa, v2j in seqs2.iteritems():
        for v, j2seq in v2j.iteritems():
            for j, seq in j2seq.iteritems():
                if aa in seqs1 and v in seqs1[aa] and j in seqs1[aa][v]:
                    continue
                x = 0
                y = seq.count
                point = (x, y)
                if point not in points:
                    points[point] = 1
                else:
                    points[point] += 1

    xdata = []
    ydata = []
    counts = []
    for point, count in points.iteritems():
        xdata.append(point[0])
        ydata.append(point[1])
        counts.append(count)
    return xdata, ydata, counts

def scatterPlot(pair, samples, options):
    isAbs = not options.rel
    sample1 = pair[0]
    sample2 = pair[1]
    seqs1 = samples[sample1].seqs
    seqs2 = samples[sample2].seqs

    #Plot:
    options.out = "%s-%s-scatter" %(sample1, sample2)
    if not isAbs:
        options.out += "-rel"
    fig, pdf = initImage(10.0, 10.0, options)

    data1, data2, counts = getData(seqs1, seqs2)
    
    if not isAbs:
        total1 = samples[sample1].total
        if total1 > 0:
            data1 = [ d*100.0/total1 for d in data1]
        total2 = samples[sample2].total
        if total2 > 0:
            data2 = [ d*100.0/total2 for d in data2]
    #else:
    nonzeros = []
    for i, d1 in enumerate(data1):
        if d1 > 0:
            nonzeros.append(d1)
        if data2[i] > 0:
            nonzeros.append(data2[i])
    minNonzeros = min(nonzeros)
    zeroReplace = float(minNonzeros)/10.0
    if not isAbs:
        zeroReplace = 1.0
        while minNonzeros < 1:
            zeroReplace /= 10.0
            minNonzeros *= 10

    for i in xrange( len(data1) ):
        if data1[i] == 0:
            data1[i] = zeroReplace
        if data2[i] == 0:
            data2[i] = zeroReplace

    axes = fig.add_axes( [0.08, 0.1, 0.87, 0.8] )
    
    xmax = max(data1)
    xmin = min(data1)
    ymax = max(data2)
    ymin = min(data2)
    xymax = max([xmax, ymax]) 
    xymin = min([xmin, ymin])
    range = xymax - xymin
    xymax += range*0.01
    #if not isAbs:
    #    xymin -= range*0.01
   
    #axes.plot( data1, data2, marker='o', markeredgecolor='b', markersize=5.0, linestyle='none' )
    maxcount = float(max(counts))
    #cm = matplotlib.cm.jet
    #counts = [ cm.jet( c/maxcount ) for c in counts]
    counts = [ cm.jet( c ) for c in counts]

    #colors = [colorConvert.to_rgb(float(c)/maxcount) for c in counts]
    cax = axes.scatter( data1, data2, color=counts, marker='o' )
    #fig.colorbar(cax, ticks=[min(counts), max(counts)], orientation='vertical')
    axes.plot( [xymin, xymax], [xymin, xymax], linestyle='-', linewidth=0.4, color="#848484" )

    #Manually draw grid lines because freaking keynote would not display the grid using axes.x/yaxis.grid
    for x in [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100]:
        axes.plot( [x, x], [xymin, 100], linestyle='-', linewidth=0.4, color='#848484')
    for y in [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100]:
        axes.plot( [xymin, 100], [y,y], linestyle='-', linewidth=0.4, color='#848484')


    #HACK
    #xymax = 1
    #END HACK
    axes.set_xlim(xymin, xymax)
    axes.set_ylim(xymin, xymax)
    
    #if isAbs:
    axes.set_xscale('log')
    axes.set_yscale('log')
    
    axes.set_title( 'Clone Size in %s versus %s' %(sample1, sample2) )
    axes.set_xlabel('%s read count' %sample1)
    axes.set_ylabel('%s read count' %sample2)
    if not isAbs:
        axes.set_xlabel('%s read frequency' %sample1)
        axes.set_ylabel('%s read frequency' %sample2)

    #axes.yaxis.grid(b=True, color="#848484", linestyle='-', linewidth=0.005)
    #axes.xaxis.grid(b=True, color="#848484", linestyle='-', linewidth=0.005)
    
    #fig.savefig( pdf, format='pdf' )
    #pdf.close()
    writeImage( fig, pdf, options)
    
    #Tab:
    #tabfile = "%s-%s-tab.tex"

def readPairsFile(file):
    pairs = []
    f = open(file, 'r')
    for line in f:
        items = line.strip().split()
        if len(items) != 2:
            sys.stderr.write("Wrong pair file format\n")
            sys.exit(1)
        pairs.append( items )
    return pairs

def main():
    usage = '%prog indir pairsFile'
    parser = OptionParser(usage = usage)
    initPlotOptions( parser )

    #parser.add_option('-p', '--percent', dest='percent', type='float', default=0.0, help='Minimum relative size (percentage of total reads) required. Default=%default')
    parser.add_option('-c', '--count', dest='count', type='int', default=1, help='Minimum number of reads required. Default=%default')
    parser.add_option('-r', '--relative', dest='rel', action='store_true', default=False, help='If specified, will draw the relative (frequency) scatter plot instead of reads count')
    options, args = parser.parse_args()
    
    checkPlotOptions( options, parser )
    
    indir = args[0]
    pairsFile = args[1]
    pairs = readPairsFile(pairsFile)
    sys.stderr.write("Done readding pairs\n")

    samples = readfiles(indir, options.count, 2)
    name2sample = {}
    for sample in samples:
        name2sample[sample.name] = sample

    for pair in pairs:
        scatterPlot(pair, name2sample, options)
        sys.stderr.write("Done plot of pair %s, %s\n" %(pair[0], pair[1]))

if __name__ == '__main__':
    main()














        

