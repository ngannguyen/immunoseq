#!/usr/bin/env python2.6

import os, sys, re, gzip
import cPickle as pickle
from immunoseq.lib.immunoseqLib import * 
#from immunoseq.lib.overlapLib import *

import matplotlib.pyplot as pyplot
from matplotlib.ticker import *
from matplotlib.font_manager import FontProperties
import numpy as np

def drawPairOverlapData_hack(axes, group2avr, group2std, group2pairVec):
    markersize = 12.0 
    offset = 0.03
    xlabels = []
    group2properLabel = {'controls':'Control - Control', 'controls-patients':'Control - Patient', 'patients':'Patient - Patient'}
    #numpoints = len( group2avr.values()[0] )
    #xoffset = offset*(numpoints-1)/2
    miny = float('inf')
    maxy = 0

    lines = []
    linenames = []

    for index, group in enumerate( sorted(group2avr.keys()) ):
        ydata = group2avr[group]
        stddata = group2std[group]
        upperydata = [ y + stddata[yi] for yi, y in enumerate(ydata) ]
        lowerydata = [ y - stddata[yi] for yi, y in enumerate(ydata) ]
        
        numpoints = len(ydata)
        xoffset = offset*(numpoints-1)/2
        i = index + 1 - xoffset
        xlabels.append(group)
        xdata = [i + offset*k for k in xrange(numpoints)]
        #axes.errorbar(xdata, ydata, yerr=stddata, color='#377EB8', markeredgecolor = "#377EB8", markersize=markersize, fmt='o')

        #Separate into red and blue:
        ydataB27pos = [] #both samples in pair are B27+
        stddataB27pos = [] 
        xdataB27pos = []

        ydataB27neg = [] #at least one of the sample in the pair is B27-
        stddataB27neg = [] 
        xdataB27neg = []
        b27posSamples = ['as1D', 'as11D', 'as16D', 'asBD', 'as20D', 'as8D']
        pairs = group2pairVec[group]
        for pairIndex, pair in enumerate(pairs):
            samples = pair.split('-')
            if samples[0] in b27posSamples and samples[1] in b27posSamples:
                ydataB27pos.append( ydata[pairIndex] )
                stddataB27pos.append( stddata[pairIndex] )
                #xdataB27pos.append( xdata[pairIndex] )
            else:
                ydataB27neg.append( ydata[pairIndex] )
                stddataB27neg.append( stddata[pairIndex] )
                #xdataB27neg.append( xdata[pairIndex] )
        
        if len(ydataB27pos) > 0:
            xdataB27pos = xdata[ : len(ydataB27pos)]
            #l = axes.errorbar(xdataB27pos, ydataB27pos, yerr=stddataB27pos, color='#E31A1C', markeredgecolor = "#E31A1C", markersize=markersize, fmt='o')
            l = axes.errorbar(xdataB27pos, ydataB27pos, yerr=stddataB27pos, color='#353535', markeredgecolor = "#353535", markersize=markersize, fmt='o')
            #Average point:
            b27posMean = [ np.mean(ydataB27pos) ]
            b27posStd = [ np.std(ydataB27pos) ]
            #avrl = axes.errorbar( [np.mean(xdataB27pos)], b27posMean, yerr=b27posStd, color='#FE8E8F', markeredgecolor='#FE8E8F', markersize=markersize + 3.0, fmt='^' )
            avrl = axes.errorbar( [np.mean(xdataB27pos)], b27posMean, yerr=b27posStd, color='#656565', markeredgecolor='#656565', markersize=markersize + 3.0, fmt='^' )
            if 'B27+, B27+' not in linenames:
                lines.append(l[0])
                linenames.append('B27+, B27+')
            if 'B27+, B27+ Avr' not in linenames:
                lines.append(avrl[0])
                linenames.append('B27+, B27+ Avr')

        if len(ydataB27neg) > 0:
            xdataB27neg = xdata[ len(ydataB27pos): ]
            l = axes.errorbar(xdataB27neg, ydataB27neg, yerr=stddataB27neg, color='#959595', markeredgecolor = "#959595", markersize=markersize, fmt='o')
            #Average point:
            b27negMean = [ np.mean(ydataB27neg) ]
            b27negStd = [ np.std(ydataB27neg) ]
            avrl = axes.errorbar( [np.mean(xdataB27neg)], b27negMean, yerr=b27negStd, color='#BDBDBD', markeredgecolor='#BDBDBD', markersize=markersize + 3.0, fmt='^' )
            if 'B27+/-, B27-' not in linenames:
                lines.append(l[0])
                linenames.append('B27+/-, B27-')
            if 'B27+/-, B27- Avr' not in linenames:
                lines.append(avrl[0])
                linenames.append('B27+/-, B27- Avr')
        
        #update max, min:
        miny = min( [min(lowerydata), miny] )
        maxy = max( [max(upperydata), maxy] )

        #Average point:
        #groupMean = [ np.mean(ydata) ]
        #groupStd = [ np.std(ydata) ]
        #axes.errorbar( [np.mean(xdata)], groupMean, yerr=groupStd, color='#4DAF4A', markeredgecolor='#4DAF4A', markersize=markersize + 3.0, fmt='^' )
    
    legend = axes.legend(lines, linenames, numpoints=1, loc="best", ncol=1)
    legend.__drawFrame = False

    editSpine(axes)
    range = maxy - miny

    #Vertical lines that separate groups:
    for x in xrange(1, len(xlabels)):
        axes.plot( [x + 0.5, x + 0.5], [miny - range*0.05, maxy + range*0.05], color="#3F3F3F", linestyle='-', linewidth=0.5)
    axes.yaxis.grid(b=True, color="#3F3F3F", linestyle='-', linewidth=0.5)
    
    axes.set_xlim(0.5, len(group2avr) + 0.5)
    axes.set_ylim(miny -range*0.05, maxy + range*0.01)
    axes.xaxis.set_ticks( xrange(1, len(xlabels) + 1) )
    axes.xaxis.set_ticklabels( [ group2properLabel[xlabel] for xlabel in xlabels ] )
    
    for label in axes.get_xticklabels():
        label.set_fontsize ('x-large')
        label.set_fontweight ('bold')
    for label in axes.get_yticklabels():
        label.set_fontsize ('large')
        label.set_fontweight ('bold')

    axes.set_title("Pairwise Overlap", size='xx-large', weight='bold')
    axes.set_ylabel("Number of shared clones", size='x-large', weight='bold')

def drawPairOverlap_hack(outfile, group2avrVec, group2stdVec, group2pairVec):
    dpi = 300
    outformat = 'pdf'
    #outformat = 'eps'
    fig, pdf = initImage2(10.0, 10.0, outformat, outfile, dpi)
    axes = setAxes(fig)
    drawPairOverlapData_hack(axes, group2avrVec, group2stdVec, group2pairVec)
    writeImage2(fig, pdf, outformat, outfile, dpi)


########### MAIN ###############
pickleFile = sys.argv[1]
outfile = sys.argv[2]
(g2avr, g2std, g2p) = pickle.load( gzip.open(pickleFile, "rb") )
drawPairOverlap_hack(outfile, g2avr, g2std, g2p)

