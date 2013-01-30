#!/usr/bin/env python2.6

'''
04/23/2012 nknguyen soe ucsc edu

'''
import os, re, sys, random
from optparse import OptionParser

import immunoseq.lib.immunoseqLib as iseqlib
import matplotlib.pyplot as pyplot
from matplotlib.ticker import *
from matplotlib.font_manager import FontProperties


def readfile(file):
    freqs = []
    f = open(file, 'r')
    for line in f:
        if len(line) >0 and line[0] == '>':
            line = line.strip()
            items = line.split(';')
            freq = items[1].split('|')[0].lstrip('freq=')
            freq = float(freq)
            freqs.append(freq)
    f.close()
    return freqs

def readfiles(indir):
    ext = 'fa'
    files = iseqlib.getfiles(indir, ext)
    sample2freqs = {}
    for file in files:
        sample = file.rstrip(ext).rstrip('.')
        freqs = readfile( os.path.join(indir, file) )
        sample2freqs[sample] = freqs
    return sample2freqs

def binDataSamples(sample2freqs):
    sample2data = {}
    for sample, freqs in sample2freqs.iteritems():
        xdata, ydata = binData(freqs)
        sample2data[sample] = (xdata, ydata)
    return sample2data

def binData(freqs):
    xdata = [0, 0.001, 0.01, 0.1, 1.0]
    ydata = [0.0 for x in xdata]
    total = len(freqs)
    if total == 0:
        return xdata, ydata
    #Sort the frequencies
    freqs = sorted(freqs)
    freq = freqs.pop(0)
    for i in xrange( len(xdata) - 1 ):
        lowerX = xdata[i]
        higherX = xdata[i+1]
        if freq < lowerX:
            break

        while freq < higherX:
            ydata[i] += 1
            if len(freqs) > 0:
                freq = freqs.pop(0)
            else:
                break
    ydata[ -1 ] = len(freqs)
    if freq >= xdata[-1]:
        ydata[-1] += 1
    ydata = [ y*100.0/total for y in ydata] 
    return xdata, ydata

def drawClonesizeDist(outdir, options, sample2data):
    options.out = os.path.join( outdir, "cloneVsRead" )
    fig, pdf = iseqlib.initImage( 10.0, 8.0, options )
    axes = fig.add_axes( [0.12, 0.15, 0.85, 0.75] )

    lines = []
    sampleNames = sorted(sample2data.keys())
    #draw
    colors = iseqlib.getColors6()
    lightColors = iseqlib.getColors6light()
    markers = ['o', '*', 's', 'd', '^', 'p', 'v']
    c = 0

    xdata = []
    xticklabels = []
    for sample in sampleNames:
        #xticklabels, ydata = binData( sample2freqs[sample] )
        (xticklabels, ydata) = sample2data[sample]
        print sample
        print xdata
        print ydata
        xdata = range(0, len(xticklabels), 1)
        markersize = 10.0
        m = markers[c]
        if m == '*':
            markersize = 12.0
        elif m == 's':
            markersize = 8.0
        
        l = axes.plot(xdata, ydata, color=colors[c], marker=m, markeredgecolor=colors[c], markersize=markersize, linestyle='none')
        axes.plot(xdata, ydata, color=lightColors[c], linestyle='-', linewidth=0.2)
        lines.append(l)
        c += 1
    
    iseqlib.editSpine(axes)
    
    axes.set_title('Clone Size Distribution', size='xx-large')
    axes.set_xlabel('Clone size as percentage of total reads', size='large')
    axes.set_ylabel('Percentage of total clones', size='large')
    
    fontP = FontProperties()
    fontP.set_size('medium')
    legend = axes.legend( lines, sampleNames, numpoints=1, loc='best', ncol=1, prop=fontP)
    legend._drawFrame = False
    
    axes.xaxis.set_ticklabels( xticklabels )
    axes.xaxis.set_ticks( [x - 0.5 for x in xdata] )
    axes.set_xlim(-0.5, len(xdata) - 0.5)
    
    axes.set_ylim(-0.5, 50)
    for label in axes.get_xticklabels():
        label.set_fontsize('medium')
        label.set_rotation(45)
    for label in axes.get_yticklabels():
        label.set_fontsize('medium')
    
    axes.yaxis.grid(b=True, color="#BDBDBD", linestyle='-', linewidth=0.005)
    axes.xaxis.grid(b=True, color="#BDBDBD", linestyle='-', linewidth=0.005)

    iseqlib.writeImage(fig, pdf, options)

def drawClonesizeDisti2(outdir, options, sample2data):
    options.out = os.path.join( outdir, "cloneVsRead" )
    fig, pdf = iseqlib.initImage( 10.0, 8.0, options )
    #axes = fig.add_axes( [0.12, 0.15, 0.85, 0.75] )
    toprange = 0.95 - 0.61
    bottomrange = 0.22 - 0
    topaxes, axes = iseqlib.setAxes3(fig, toprange, bottomrange)

    lines = []
    sampleNames = sorted(sample2data.keys())
    #draw
    colors = iseqlib.getColors6()
    lightColors = iseqlib.getColors6light()
    markers = ['o', '*', 's', 'd', '^', 'p', 'v']
    c = 0

    xdata = []
    xticklabels = []
    for sample in sampleNames:
        #xticklabels, ydata = binData( sample2freqs[sample] )
        (xticklabels, ydata) = sample2data[sample]
        print sample
        print xdata
        print ydata
        xdata = range(0, len(xticklabels), 1)
        markersize = 10.0
        m = markers[c]
        if m == '*':
            markersize = 12.0
        elif m == 's':
            markersize = 8.0
        #Upper plot
        l = topaxes.plot(xdata, ydata, color=colors[c], marker=m, markeredgecolor=colors[c], markersize=markersize, linestyle='none')
        topaxes.plot(xdata, ydata, color=lightColors[c], linestyle='-', linewidth=0.2)
        lines.append(l)
        
        #Lower plot
        axes.plot(xdata, ydata, color=colors[c], marker=m, markeredgecolor=colors[c], markersize=markersize, linestyle='none')
        axes.plot(xdata, ydata, color=lightColors[c], linestyle='-', linewidth=0.2)
        c += 1
    
    iseqlib.editSpine2(topaxes, axes)
    
    topaxes.set_title('Clone Size Distribution', size='xx-large')
    axes.set_xlabel('Clone size as percentage of total reads', size='large')
    axes.set_ylabel('Percentage of total clones', size='large')
    
    fontP = FontProperties()
    fontP.set_size('medium')
    legend = topaxes.legend( lines, sampleNames, numpoints=1, loc='best', ncol=1, prop=fontP)
    legend._drawFrame = False
    
    #Draw the discontinue sign:
    iseqlib.drawDiscontinueSign(topaxes, axes, 61, 22)

    axes.xaxis.set_ticklabels( xticklabels )
    axes.xaxis.set_ticks( [x - 0.5 for x in xdata] )
    axes.set_xlim(-0.5, len(xdata) - 0.5)
    
    dummyxticklabels = [ "" for l in xticklabels ]
    topaxes.set_xticklabels(dummyxticklabels)
    topaxes.xaxis.set_ticks( [x - 0.5 for x in xdata] )
    topaxes.set_xlim(-0.5, len(xdata) - 0.5)
    topaxes.set_ylim(61, 95)
    
    axes.set_ylim(-0.5, 22)
    for label in axes.get_xticklabels():
        label.set_fontsize('medium')
        label.set_rotation(45)
    for label in axes.get_yticklabels():
        label.set_fontsize('medium')

    topaxes.yaxis.grid(b=True, color="#BDBDBD", linestyle='-', linewidth=0.005)
    topaxes.xaxis.grid(b=True, color="#BDBDBD", linestyle='-', linewidth=0.005)
    
    axes.yaxis.grid(b=True, color="#BDBDBD", linestyle='-', linewidth=0.005)
    axes.xaxis.grid(b=True, color="#BDBDBD", linestyle='-', linewidth=0.005)

    iseqlib.writeImage(fig, pdf, options)

def samplingSamples(sample2freqs, numSamplings):
    size = min( [ len(freqs) for freqs in sample2freqs.values() ] )
    sample2data = {}
    for sample, freqs in sample2freqs.iteritems(): #each sample
        xdata = []
        ydata = []
        for i in xrange(numSamplings):
            ifreqs = random.sample( freqs, size )
            ixdata, iydata = binData(ifreqs)
            if len(xdata) == 0:
                xdata = ixdata
                ydata = iydata
            else:
                for j, y in enumerate( iydata ):
                    ydata[j] += y
        for j in xrange( len(ydata) ):#average ydata
            ydata[j] /= numSamplings
        sample2data[sample] = (xdata, ydata)
    return sample2data

def main():
    usage=('Usage: %prog [options] inputDir outputDir')
    parser = OptionParser(usage=usage)
    parser.add_option('-s', '--sampling', dest='sampling', action='store_true', default=False, help='If specified, will randomly sampling sequences from each file with sampling size equals to size of the smallest sample')
    parser.add_option('--numSamplings', dest='numSamplings', type='int', default=100, help='Number of samplings. Default=%default')
    iseqlib.initPlotOptions(parser)
    options, args = parser.parse_args()
    iseqlib.checkPlotOptions(options, parser)
    
    indir = args[0]
    outdir = args[1]
    sample2freqs = readfiles(indir)
    if options.sampling:
        sample2data = samplingSamples(sample2freqs, options.numSamplings)
    else:
        sample2data = binDataSamples(sample2freqs) 

    drawClonesizeDist(outdir, options, sample2data)

if __name__ == "__main__":
    main()




