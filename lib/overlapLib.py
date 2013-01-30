#nknguyen soe ucsc edu
#Wed Jul 18 16:56:23 PDT 2012

import os, sys, re, copy, time
import random as rnd
from immunoseq.lib.immunoseqLib import * 
from sonLib.bioio import system

import matplotlib.pyplot as pyplot
from matplotlib.ticker import *
from matplotlib.font_manager import FontProperties
from scipy.stats import ttest_ind
import numpy as np

###########DRAW PAIRWISE OVERLAP STATS ###################
#######HACK TO COLOR B27+ PAIRS RED #######
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
            avrl = axes.errorbar( [np.mean(xdataB27pos)], b27posMean, yerr=b27posStd, color='#BDBDBD', markeredgecolor='#BDBDBD', markersize=markersize + 3.0, fmt='^' )
            if 'B27+, B27+' not in linenames:
                lines.append(l[0])
                linenames.append('B27+, B27+')
            if 'B27+, B27+ Avr' not in linenames:
                lines.append(avrl[0])
                linenames.append('B27+, B27+ Avr')

        if len(ydataB27neg) > 0:
            xdataB27neg = xdata[ len(ydataB27pos): ]
            l = axes.errorbar(xdataB27neg, ydataB27neg, yerr=stddataB27neg, color='#848484', markeredgecolor = "#848484", markersize=markersize, fmt='o')
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
    fig, pdf = initImage2(10.0, 10.0, outformat, outfile, dpi)
    axes = setAxes(fig)
    drawPairOverlapData_hack(axes, group2avrVec, group2stdVec, group2pairVec)
    writeImage2(fig, pdf, outformat, outfile, dpi)
########### END HACK ###############

def drawPairOverlapData(axes, group2avr, group2std):
    markersize = 10.0 
    offset = 0.02
    xlabels = []
    #numpoints = len( group2avr.values()[0] )
    #xoffset = offset*(numpoints-1)/2
    miny = float('inf')
    maxy = 0
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
        axes.errorbar(xdata, ydata, yerr=stddata, color='#377EB8', markeredgecolor = "#377EB8", markersize=markersize, fmt='o')
        #axes.errorbar(xdata, ydata, yerr=stddata, color='#A6D7FE', markeredgecolor = "#A6D7FE", markersize=markersize, fmt='o')
        #axes.plot(xdata, ydata, color='#377EB8', marker='o', markeredgecolor = "#377EB8", markersize=markersize, linestyle='none')

        #update max, min:
        miny = min( [min(lowerydata), miny] )
        maxy = max( [max(upperydata), maxy] )

        #Average point:
        groupMean = [ np.mean(ydata) ]
        groupStd = [ np.std(ydata) ]
        #axes.errorbar( [np.mean(xdata)], groupMean, yerr=groupStd, color='#275880', markeredgecolor='#275880', markersize=markersize + 3.0, fmt='^' )
        axes.errorbar( [np.mean(xdata)], groupMean, yerr=groupStd, color='#E31A1C', markeredgecolor='#E31A1C', markersize=markersize + 3.0, fmt='^' )
    
    editSpine(axes)
    range = maxy - miny

    #Vertical lines that separate groups:
    for x in xrange(1, len(xlabels)):
        axes.plot( [x + 0.5, x + 0.5], [miny - range*0.05, maxy + range*0.05], color="#848484", linestyle='-', linewidth=0.005)
    axes.yaxis.grid(b=True, color="#848484", linestyle='-', linewidth=0.005)
    
    axes.set_xlim(0.5, len(group2avr) + 0.5)
    axes.set_ylim(miny -range*0.05, maxy + range*0.01)
    axes.xaxis.set_ticks( xrange(1, len(xlabels) + 1) )
    axes.xaxis.set_ticklabels( xlabels )
    
    for label in axes.get_xticklabels():
        label.set_fontsize ('x-large')
    
    axes.set_title("Pairwise Overlap", size='xx-large')
    axes.set_ylabel("Number of clones", size='x-large')

def drawPairOverlap(outfile, group2avrVec, group2stdVec):
    dpi = 300
    outformat = 'pdf'
    fig, pdf = initImage2(10.0, 10.0, outformat, outfile, dpi)
    axes = setAxes(fig)
    drawPairOverlapData(axes, group2avrVec, group2stdVec)
    writeImage2(fig, pdf, outformat, outfile, dpi)

###########################################

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

def getPairwiseOverlap( seqs1, seqs2, aa2v2j1, aa2v2j2, cutoffs, mode, discrete ):
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
        s1 = seqs1[k]
        #s2 = seqs2[k]
        s2 = findSeq(s1, aa2v2j2)
        if not s2: #not found in repertoire 2
            continue

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

    #printPairwiseOverlap(reads1, reads2, clones1, clones2, stats1, stats2, cutoffs, outfile)
    return reads1, reads2, clones1, clones2, stats1, stats2

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

def printPairwiseOverlapSequences(name1, name2, seqs1, seqs2, aa2v2j1, aa2v2j2, outdir, cutoffs, mode, discrete):
    cutoff2seqs1 = {}
    cutoff2seqs2 = {}
    #Initialize cutoff2seqs:
    for c in cutoffs:
        cutoff2seqs1[c] = []
        cutoff2seqs2[c] = []

    for k in seqs1:
        s1 = seqs1[k]
        if not hasSeq(s1, aa2v2j2):
            continue
        #s2 = seqs2[k]
        s2 = findSeq(s1, aa2v2j2)
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
#def fitPower(xdata, ydata):
#    #Try to fit a Y=aX^b to the observed data. i.e try to estimate a and b
def sample2colorFixed():
    colors = ["#E31A1C", "#FE8E8F", "#377EB8", "#A6D7FE", "#4DAF4A", "#B8FEB5", "#984EA3", "#F6BDFE",
              "#FF7F00", "#FEBF80", 
              "#525252", "#737373", "#969696", "#BDBDBD",
              #"#006837", "#31A354", "#78C679", "#C2E699",
              "#2B8CBE", "#7BCCC4", "#BAE4BC",
              "#225EA8", "#41B6C4", "#A1DAB4",
              "#253494",
              "#993404", "#D95F0E", 
              "#6A51A3", "#9E9AC8", "#CBC9E2"
             ]
    s2i = {#'as10R':0, 'as11R': 1, 'as12R': 2, 'as13R': 3, 'as1R': 4, 'as1D': 5, 'as8D': 6, 'as15D': 7,
           #'asBR': 8, 'asBD': 9,
           'adaptBSD1': 10, 'adaptBSD1cd8n': 11, 'adaptBSD2R': 12, 'adaptBSD3': 13,
           'adaptF28cd4': 14, 'adaptF57cd4': 15, 'adaptM35cd4': 16, 
           'adaptF28cd8': 17, 'adaptF57cd8': 18, 'adaptM35cd8': 19, 
           'adaptF28': 14, 'adaptF57': 15, 'adaptM35': 16, 
           'adaptRep': 20,
           'wangAll': 21, 'wangCyt': 22, 
           'warrenM1': 23, 'warrenM2': 24, 'warrenF': 25,
           'f28CD4memory': 0, 'f28CD4naive': 1, 'f28CD8memory': 2, 'f28CD8naive': 3, 
           'f57CD4memory': 4, 'f57CD4naive': 5, 'f57CD8memory': 6, 'f57CD8naive': 7, 
           'm35CD4memory': 8, 'm35CD4naive': 9, 'm35CD8memory': 10, 'm35CD8naive': 11,
           #'as11D':12, 'as16D':13, 'as1Ddraw2':21, 'as1Ddraw2notcd8':23,
           #'asBDdraw2':24, 'as20D':25,
           'as1D':0, 'as8D':1, 'as11D':2, 'as15D':3, 'as16D':4,
           'asBD':5, 'as20D':6,
           'adaptAS':7, 'adaptCD':8, 'adaptMA':9, 'adaptBSD1n':11,
           'as1Ddraw2':1, 'as1R':2, 'irep1D': 3, 'irep1R':4,
          }
    s2c = {}
    for s, i in s2i.iteritems():
        s2c[s] = colors[i]
    return s2c

def sample2color(names):
    s2cfixed = sample2colorFixed()
    sample2color = {}

    #colors = iseqlib.getColors0()
    colors = getColors0()
    for i, name in enumerate(names):
        if name in s2cfixed:
            sample2color[name] = s2cfixed[name]
        else:
            sample2color[name] = colors[i]
    return sample2color 

def drawOverlapPlotData(data, axes, labels, sam2color):
    #colors = iseqlib.getColors0()
    lines = []

    maxy = 0
    maxx = 0
    for i, d in enumerate(data):
        xdata = d[0]
        ydata = d[1]
        maxy = max([maxy, max(ydata)])
        maxx = max([maxx, max(xdata)])
        name = labels[i]
        #l = axes.plot(xdata, ydata, color=colors[i], linestyle='-', linewidth=3)
        l = axes.plot(xdata, ydata, color=sam2color[name], linestyle='-', linewidth=3)
        lines.append(l)
    #HACK:
    #axes.set_xlim(0, 1.5*(10**9))
    #axes.set_xlim(0, 3.5*(10**9))
    axes.set_ylim(0, max([maxy*0.5, 50]) )
    axes.set_xlim(0, maxx*0.5)

    legend = pyplot.legend(lines, labels, numpoints=1, loc='best')
    #axes.set_xlabel('Product of number of top unique sequences n_i1 x n_i2', size='x-large')
    axes.set_title("Number of shared clones over an increasing sampling space") 
    axes.set_xlabel('n1 x n2', size='x-large')
    axes.set_ylabel('Number of shared unique sequences', size='x-large')

def drawOverlapPlot(data, file, labels, options, sam2color):
    options.out = file
    fig, pdf = initImage(8.0, 10.0, options)
    axes = setAxes(fig)
    drawOverlapPlotData(data, axes, labels, sam2color)
    writeImage(fig, pdf, options)

def addSeq2aa2v2j(seq, aa2v2j):
    v2j = {}
    if seq.seq in aa2v2j:
        v2j = aa2v2j[seq.seq]
    else:
        aa2v2j[seq.seq] = v2j

    for v in seq.vs:
        j2seq = {}
        if v in v2j:
            j2seq = v2j[v]
        else:
            v2j[v] = j2seq

        for j in seq.js:
            if j not in j2seq:
                j2seq[j] = seq

def cmpseq(seq1, seq2):
    if seq1.seq == seq2.seq: #same amino acid sequence
        for v1 in seq1.vs:
            if v1 in seq2.vs:
                for j1 in seq1.js: 
                    if j1 in seq2.js: #at least one (v, j) pair of seq1 presents in seq2
                        return True
    return False

#def pairwiseOverlapPlot(seqs1, seqs2, outfile): 
def pairwiseOverlapPlot(seqs1, seqs2, isRandom): 
    #sorting seqs1
    seqs1list = [(header, seq) for header, seq in seqs1.iteritems()]
    seqs2list = [(header, seq) for header, seq in seqs2.iteritems()]
    if not isRandom: 
        seqs1list = sorted( seqs1list, key=lambda item:item[1].count, reverse=True )
        seqs2list = sorted( seqs2list, key=lambda item:item[1].count, reverse=True )
    else:
        rheaders1 = rnd.sample( seqs1.keys(), len(seqs1.keys()) )
        rheaders2 = rnd.sample( seqs2.keys(), len(seqs2.keys()) )
        seqs1list = [(header, seqs1[header]) for header in rheaders1]
        seqs2list = [(header, seqs2[header]) for header in rheaders2]

    top1 = {} #aa2v2j
    top2 = {}
    xdata = []
    ydata = []
    currx = 0
    curry = 0
    
    for i in xrange( max([len(seqs1list), len(seqs2list)]) ):
        if i >= len(seqs1list):
            (h1, s1) = ('', None)
        else:
            (h1, s1) = seqs1list[i]
        
        if i >= len(seqs2list):
            (h2, s2) = ('', None)
        else:
            (h2, s2) = seqs2list[i]

        if h1 != '' and h2 != '' and cmpseq(s1, s2): #s1 and s2 are the same sequence
            curry += 1
        else: #s1 and s2 are not the same sequence, check to see if s1 is already in top2, and s2 is already in top1
            if h1 != '' and hasSeq(s1, top2):
                curry += 1
            if h2 != '' and hasSeq(s2, top1):
                curry += 1
        
        #Update top1 and top2
        if h1 != '':
            addSeq2aa2v2j(s1, top1)
        if h2 != '':
            addSeq2aa2v2j(s2, top2)

        i1 = min([i + 1, len(seqs1list)]) 
        i2 = min([i + 1, len(seqs2list)])
        currx = i1*i2
        xdata.append( currx )
        ydata.append(curry)
    
    return xdata, ydata
################## END PAIRWISE OVERLAP PLOT (ROBINS et al.) #####################

#############################################################################
################## DRAW OVERLAP PLOT COMBINED (ALL) #########################
#############################################################################
def countOverlap(seqs1, aa2v2j2):
    numSharedSeqs = 0
    for h1, s1 in seqs1.iteritems():
        if hasSeq(s1, aa2v2j2):
            numSharedSeqs += 1
    return numSharedSeqs

def drawOverlapPlotAllData(sample2caseYdata, sample2controlYdata, axesList, group2samples, sample2group, group2color, samplesPerPlot):
    sampleNames = []
    for group in sorted(group2samples.keys()):
        sampleNames.extend( group2samples[group] )

    markers = ['^', 'o']
    markersize = 8.0
    lines = []
    offset = 0.02 #x offset between data points
    
    miny = float('inf')
    maxy = 0

    for i in xrange( len(axesList) ):
        xlabels = []
        axes = axesList[i]
        startIndex = i*samplesPerPlot
        endIndex = min( [startIndex + samplesPerPlot, len(sampleNames)] )
        for j in xrange(startIndex, endIndex):
            name = sampleNames[j]
            xlabels.append(name)
            caseYdata = sample2caseYdata[name]
            xdata = [ j - startIndex + offset*k for k in xrange( len(caseYdata) ) ] 
            color = group2color[ sample2group[name] ][0]
            darkcolor = group2color[ sample2group[name] ][1]
            l1 = axes.plot(xdata, caseYdata, color= darkcolor, marker=markers[0], markeredgecolor= darkcolor , markersize=markersize, linestyle='none')

            controlYdata = sample2controlYdata[name]
            #xdata = [ j - startIndex for y in controlYdata ] 
            xdata = [ j - startIndex + offset*(k + len(caseYdata)) for k in xrange( len(controlYdata) ) ] 
            l2 = axes.plot(xdata, controlYdata, color= color, marker=markers[1], markeredgecolor=color, markersize=markersize, linestyle='none')
            if j == 0:
                lines = [l1, l2]

            miny = min( [miny, min(caseYdata), min(controlYdata)] )
            maxy = max( [maxy, max(caseYdata), max(controlYdata)] )
        
        editSpine(axes)
        axes.xaxis.set_ticklabels( xlabels )
        numSamples = len(xlabels)
        xoffset = offset*(numSamples-2)/2
        range = maxy - miny
        
        #Draw vertical grid:
        for x in xrange(numSamples - 1):
            axes.plot( [x + xoffset + 0.5, x + xoffset + 0.5], [miny - range*0.05, maxy + range*0.05], color="#848484", linestyle='-', linewidth=0.005)

        axes.set_xlim(xoffset - 0.5, min([samplesPerPlot, len(sampleNames)]) + xoffset - 0.5 )
        axes.xaxis.set_ticks( [ x + xoffset for x in xrange(numSamples) ] )
        
        #HACK:
        #yticks = [ float(y)/(10**7) for y in xrange(2, 11, 2) ]
        #yticklabels = [ str(y) for y in xrange(2, 11, 2) ]
        #yticks = [ float(y)/(10**8) for y in xrange(2, 21, 2) ]
        #yticklabels = [ "%.2f" % (y*(10**7)) for y in yticks ]
        #axes.yaxis.set_ticks(yticks)
        #axes.yaxis.set_ticklabels( yticklabels )
        #for l in axes.get_yticklabels():
        #    l.set_fontsize('medium')

        axes.set_ylim( miny - range*0.01 , maxy + range*0.01 )
        #axes.set_ylim(-0.00000001 , 0.0000003)
        #axes.set_ylim(-0.0000001 , 0.000005)

        for label in axes.get_xticklabels():
            label.set_fontsize( 'medium' )
            label.set_rotation( 45 )
        #axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
        axes.yaxis.grid(b=True, color="#848484", linestyle='-', linewidth=0.005)
        
        if i == 0:
            axes.set_title("Shared sequences", size='xx-large')
            #HACK
            #legend = axes.legend(lines, ["RNA", "DNA"], numpoints=1, loc="best", ncol=1)
            legend = axes.legend(lines, ["AS", "Healthy"], numpoints=1, loc="best", ncol=1)
        #if i == len(axesList) /2:
        #    axes.set_ylabel("Ratio of shared sequences to sampling space (x 10^-7)", size='large')
        if i == len(axesList) - 1:
            axes.set_xlabel("Samples", size='large')
        
    legend.__drawFrame = False


def drawOverlapPlotAll(samples, sample2aa2v2j, cases, controls, options, group2samples, sample2group): 
    name2sample = {}
    for s in samples:
        name2sample[s.name] = s

    #Get data
    sample2caseYdata = {}
    sample2controlYdata = {}
    for s1 in samples:#each sample
        sample2caseYdata[s1.name] = []
        #aa2v2j1 = sample2aa2v2j[s1.name]

        for name in cases:       
            if s1.name == name:
                continue
            s2 = name2sample[name]
            aa2v2j2 = sample2aa2v2j[s2.name]
            y = countOverlap(s1.seqs, aa2v2j2)*1.0/( len(s1.seqs)*len(s2.seqs) )
            sample2caseYdata[s1.name].append( y )
        
        sample2controlYdata[s1.name] = []
        for name in controls:
            if s1.name == name:
                continue
            s2 = name2sample[name]
            aa2v2j2 = sample2aa2v2j[s2.name]
            y = countOverlap(s1.seqs, aa2v2j2)*1.0/( len(s1.seqs)*len(s2.seqs) )
            sample2controlYdata[s1.name].append( y )

    #Draw plot:
    options.out = os.path.join(options.outdir, "overlapAll")
    fig, pdf = initImage(10.0, 10.0, options)
    axesList = setAxes2(fig, len( sample2group.keys() ), options.samplesPerPlot)
    colors = getColors6()
    colorsDark = getColors6dark()
    group2color = {}
    for i, group in enumerate( sorted( group2samples.keys() ) ):
        group2color[group] = (colors[i + 1], colorsDark[i + 1])

    drawOverlapPlotAllData(sample2caseYdata, sample2controlYdata, axesList, group2samples, sample2group, group2color, options.samplesPerPlot)
    writeImage(fig, pdf, options)

################## END DRAW OVERLAP PLOT COMBINED (ALL) #########################

def countOverlap_3way(seqs1, aa2v2j2, aa2v2j3):
    numSharedSeqs = 0
    for h1, s1 in seqs1.iteritems():
        if hasSeq(s1, aa2v2j2) and hasSeq(s1, aa2v2j3):
            numSharedSeqs += 1
    return numSharedSeqs

def drawOverlapPlotAllData_3way(sample2casecaseYdata, sample2controlcontrolYdata, sample2casecontrolYdata, axesList, group2samples, sample2group, group2color, samplesPerPlot):
    sampleNames = []
    for group in sorted(group2samples.keys()):
        sampleNames.extend( group2samples[group] )

    markers = ['^', 'o', 'd']
    markersize = 8.0
    lines = []
    offset = 0.02 #x offset between data points
    
    miny = float('inf')
    maxy = 0

    for i in xrange( len(axesList) ):
        xlabels = []
        axes = axesList[i]
        startIndex = i*samplesPerPlot
        endIndex = min( [startIndex + samplesPerPlot, len(sampleNames)] )
        numx = 0
        for j in xrange(startIndex, endIndex):
            name = sampleNames[j]
            xlabels.append(name)
            
            lightcolor = group2color[ sample2group[name] ][0]
            color = group2color[ sample2group[name] ][1]
            darkcolor = group2color[ sample2group[name] ][2]
            
            #casecase
            casecaseYdata = sample2casecaseYdata[name]
            xdata = [ j - startIndex + offset*k for k in xrange( len(casecaseYdata) ) ] 
            l1 = axes.plot(xdata, casecaseYdata, color= darkcolor, marker=markers[0], markeredgecolor= darkcolor , markersize=markersize, linestyle='none')

            #casecontrol
            casecontrolYdata = sample2casecontrolYdata[name]
            xdata = [ j - startIndex + offset*(k + len(casecaseYdata)) for k in xrange( len(casecontrolYdata) ) ] 
            l2 = axes.plot(xdata, casecontrolYdata, color = color, marker=markers[1], markeredgecolor=color, markersize=markersize, linestyle='none')
            
            controlcontrolYdata = sample2controlcontrolYdata[name]
            #xdata = [ j - startIndex for y in controlYdata ] 
            xdata = [ j - startIndex + offset*( k + len(casecontrolYdata) + len(casecaseYdata) ) for k in xrange( len(controlcontrolYdata) ) ] 
            l3 = axes.plot(xdata, controlcontrolYdata, color= lightcolor, marker=markers[2], markeredgecolor= lightcolor, markersize=markersize, linestyle='none')
            if j == 0:
                lines = [l1, l2, l3]

            miny = min( [miny, min(casecaseYdata), min(casecontrolYdata), min(controlcontrolYdata)] )
            maxy = max( [maxy, max(casecaseYdata), max(casecontrolYdata), max(controlcontrolYdata)] )
            if numx == 0:
                numx = len(casecaseYdata) + len(casecontrolYdata) + len(controlcontrolYdata)
        
        editSpine(axes)
        axes.xaxis.set_ticklabels( xlabels )
        numSamples = len(xlabels)
        xoffset = offset*(numx-1)/2
        range = maxy - miny
        
        #Draw vertical grid:
        for x in xrange(numSamples - 1):
            axes.plot( [x + xoffset + 0.5, x + xoffset + 0.5], [miny - range*0.05, maxy + range*0.05], color="#848484", linestyle='-', linewidth=0.005)

        axes.set_xlim(xoffset - 0.5, min([samplesPerPlot, len(sampleNames)]) + xoffset - 0.5 )
        axes.xaxis.set_ticks( [ x + xoffset for x in xrange(numSamples) ] )
        
        #HACK:
        #yticks = [ float(y)/(10**7) for y in xrange(2, 11, 2) ]
        #yticklabels = [ str(y) for y in xrange(2, 11, 2) ]
        #yticks = [ float(y)/(10**8) for y in xrange(2, 21, 2) ]
        #yticklabels = [ "%.2f" % (y*(10**7)) for y in yticks ]
        #axes.yaxis.set_ticks(yticks)
        #axes.yaxis.set_ticklabels( yticklabels )
        #for l in axes.get_yticklabels():
        #    l.set_fontsize('medium')

        #axes.set_ylim( miny - range*0.01 , maxy + range*0.01 )

        for label in axes.get_xticklabels():
            label.set_fontsize( 'medium' )
            label.set_rotation( 45 )
        #axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
        axes.yaxis.grid(b=True, color="#848484", linestyle='-', linewidth=0.005)
        
        if i == 0:
            axes.set_title("Shared sequences", size='xx-large')
            #HACK
            #legend = axes.legend(lines, ["RNA", "DNA"], numpoints=1, loc="best", ncol=1)
            legend = axes.legend(lines, ["AS-AS", "AS-Healthy", "Healthy-Healthy"], numpoints=1, loc="best", ncol=1)
        #if i == len(axesList) /2:
        #    axes.set_ylabel("Ratio of shared sequences to sampling space (x 10^-7)", size='large')
        if i == len(axesList) - 1:
            axes.set_xlabel("Samples", size='large')
        
    legend.__drawFrame = False

def drawOverlapPlotAll_3way(samples, sample2aa2v2j, cases, controls, options, group2samples, sample2group): 
    name2sample = {}
    for s in samples:
        name2sample[s.name] = s
    
    #for each sample, there are going to be 3 types of data: case-case, case-control, control-control
    casecasePairs = []
    casecontrolPairs = []
    controlcontrolPairs = []
    for i in xrange( len(cases) -1 ):
        for j in xrange(i + 1, len(cases)):
            casecasePairs.append( [cases[i], cases[j]] )

    for i in xrange( len(controls) -1 ):
        for j in xrange(i + 1, len(controls)):
            controlcontrolPairs.append( [controls[i], controls[j]] )

    for case in cases:
        for control in controls:
            casecontrolPairs.append( [case, control] )

    #Get data
    normFactor = 10**10
    sample2casecaseYdata = {}
    sample2casecontrolYdata = {}
    sample2controlcontrolYdata = {}
    for s in samples:#each sample
        sample2casecaseYdata[s.name] = []
        for [name1, name2] in casecasePairs:
            if s.name == name1 or s.name == name2:
                continue
            s1 = name2sample[name1]
            s2 = name2sample[name2]
            aa2v2j1 = sample2aa2v2j[s1.name]
            aa2v2j2 = sample2aa2v2j[s2.name]
            y = countOverlap_3way(s.seqs, aa2v2j1, aa2v2j2)
            y = y*1.0*normFactor/( len(s.seqs)*len(s1.seqs)*len(s2.seqs) )
            sample2casecaseYdata[s.name].append( y )
        
        sample2controlcontrolYdata[s.name] = []
        for [name1, name2] in controlcontrolPairs:
            if s.name == name1 or s.name == name2:
                continue
            s1 = name2sample[name1]
            s2 = name2sample[name2]
            aa2v2j1 = sample2aa2v2j[s1.name]
            aa2v2j2 = sample2aa2v2j[s2.name]
            y = countOverlap_3way(s.seqs, aa2v2j1, aa2v2j2)
            y = y*1.0*normFactor/( len(s.seqs)*len(s1.seqs)*len(s2.seqs) )
            sample2controlcontrolYdata[s.name].append( y )
        
        sample2casecontrolYdata[s.name] = []
        for [name1, name2] in casecontrolPairs:
            if s.name == name1 or s.name == name2:
                continue
            s1 = name2sample[name1]
            s2 = name2sample[name2]
            aa2v2j1 = sample2aa2v2j[s1.name]
            aa2v2j2 = sample2aa2v2j[s2.name]
            y = countOverlap_3way(s.seqs, aa2v2j1, aa2v2j2)
            y = y*1.0*normFactor/( len(s.seqs)*len(s1.seqs)*len(s2.seqs) )
            sample2casecontrolYdata[s.name].append( y )
        
    #Draw plot:
    options.out = os.path.join(options.outdir, "overlapAll_3way")
    fig, pdf = initImage(10.0, 10.0, options)
    axesList = setAxes2(fig, len( sample2group.keys() ), options.samplesPerPlot)
    colors = getColors6()
    colorsDark = getColors6dark()
    colorsLight = getColors6light()
    group2color = {}
    for i, group in enumerate( sorted( group2samples.keys() ) ):
        group2color[group] = (colorsLight[i+1], colors[i + 1], colorsDark[i + 1])

    drawOverlapPlotAllData(sample2casecaseYdata, sample2controlcontrolYdata, axesList, group2samples, sample2group, group2color, options.samplesPerPlot)
    #drawOverlapPlotAllData_3way(sample2casecaseYdata, sample2controlcontrolYdata, sample2casecontrolYdata, axesList, group2samples, sample2group, group2color, options.samplesPerPlot)
    writeImage(fig, pdf, options)

################## END DRAW OVERLAP PLOT COMBINED (ALL) 3 WAYS #########################



################# GET FASTA FILES OF UNIQUE SEQUENCES, SHARED WITH SAME GROUP SEQUENCES, SHARED WITH OTHER GROUPS SEQUENCES, SHARED SEQUENCES ###################
def convertHeaderNt2aa(header):
    items = header.split('|')
    aa = nt2aa(items[0])
    return '|'.join( [ aa, items[1], items[2] ] )

def getSeq2sample2count(samples, nt2aa): #key = seq|vs|js, val = {sample: (count, freq)}
    seq2sample2count = {}
    for sample in samples:
        for header, seq in sample.seqs.iteritems():
            if nt2aa:
                header = convertHeaderNt2aa(header)

            if header not in seq2sample2count:
                seq2sample2count[header] = { sample.name: [ (seq.count, seq.freq) ] }
            else:
                if sample.name not in seq2sample2count[header]:
                    seq2sample2count[header][sample.name] = [ (seq.count, seq.freq) ]
                else:
                    seq2sample2count[header][sample.name].append( (seq.count, seq.freq) )
    return seq2sample2count

def uniqAndSharedSequences(outdir, sample, seq2sample2count, group2samples, group, nt2aa):
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
        if nt2aa:
            header = convertHeaderNt2aa(header)

        if header not in seq2sample2count:
            raise ValueError("Sequence %s of sample %s is not in seq2sample2count\n" %(header, sample.name))
        sample2count = seq2sample2count[header]
        fastaheader = "%s;freq=%f|%s|%s;size=%d" %(sample.name, seq.freq, ','.join(seq.vs), ','.join(seq.js), seq.count)
        samegroup = False
        diffgroup = False
        
        #Search to see if current sequence present in other samples of the same group and of different group(s)
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


################## TTEST ######################
#ttest to test the hypothesis: cases have more shared sequences than a case and a control or two controls.
#Calculate the below shared sequences stats for all posible pairs of (case-case) and compare with stats of all pairs of (case-control) and (control-control)
#For each pair, the stats is calculated as follow:
#Steps involved: 1/ randomly select N (ttestSamplingSize) uniq sequences from each sample
#                2/ calculate the number of shared sequences between the two samples
#                3/ Repeat (1) and (2) M times (M = ttestSamplingNum)
#                4/ Take the average of all the samplings
def getAvrSharedSeqs(sample1, sample2, samplingSize, numSamplings): 
    stime = time.clock()
    data = []
    
    maxTotal = max( [ sum([seq.count for seq in sample.seqs.values()])  for sample in [sample1, sample2] ] )
    
    for i in xrange(numSamplings): #each sampling:
        samplingtime = time.clock()

        s1 = sample1
        s2 = sample2
        sys.stderr.write("%s, %s\n" %(s1.name, s2.name))
        if samplingSize > 0:
            s1 = samplingSample_uniq(sample1, samplingSize, maxTotal)
            s2 = samplingSample_uniq(sample2, samplingSize, maxTotal)
        #seqs1 = sample1.seqs.values()
        #seqs2 = sample2.seqs.values()
        seqs1 = s1.seqs
        seqs2 = s2.seqs
        sys.stderr.write("Number of uniq sequences in sample %s : %d, %s: %d\n" %(s1.name, len(seqs1), s2.name, len(seqs2)) )
        shared = 0
        for seq1 in seqs1:
            if seq1 in seqs2:
                shared += 1
        data.append(shared)

        sys.stderr.write("One sampling in %f s.\n" %(time.clock() - samplingtime))
    sys.stderr.write("Calculate shared sequences for %d samplings took %f seconds.\n" %(numSamplings, time.clock() - stime))
    return np.mean(data)

def ttest( outfile, samples, cases, controls, samplingSize, numSamplings, mode ):
    #if mode == 1: group1 includes: case-case; group2 includes control-control, control-case
    #elif mode == 2: group2 includes: case-case and control-control, group2 includes control-case
    #

    if len(cases) < 2 or len(controls) == 0 :
        raise ValueError("Cannot perform ttest with < 2 cases or 0 control.\n")

    casePairs = [] #stats of case-case pairs
    controlPairs = [] #stats of case-control and control-control pairs
    caseNames = []
    controlNames = []

    #Convert samples two a hash for easy look up:
    name2sample = {} #key = name, val = Sample
    for s in samples:
        name2sample[s.name] = s

    #Case-case pairs:
    for i in xrange( len(cases) -1 ):
        case1 = cases[i]
        sample1 = name2sample[case1]
        for j in xrange(i+1, len(cases)):
            case2 = cases[j]
            #HACK
            #if sample2patient[case1] == sample2patient[case2]:
            #    continue
            sample2 = name2sample[case2]
            sys.stderr.write("Sample1: %s, sample2: %s\n" %(sample1.name, sample2.name))
            casePairs.append( getAvrSharedSeqs(sample1, sample2, samplingSize, numSamplings) )
            caseNames.append("%s-%s" %(case1, case2))

    #Control-control pairs:
    for i in xrange( len(controls) -1 ):
        control1 = controls[i]
        sample1 = name2sample[control1]
        for j in xrange(i+1, len(controls)):
            control2 = controls[j]
            #HACK
            #if sample2patient[control1] == sample2patient[control2]:
            #    continue
            sample2 = name2sample[control2]
            if mode == 1:
                controlPairs.append( getAvrSharedSeqs(sample1, sample2, samplingSize, numSamplings) )
                controlNames.append("%s-%s" %(control1, control2))
            elif mode == 2:
                casePairs.append( getAvrSharedSeqs(sample1, sample2, samplingSize, numSamplings) )
                caseNames.append("%s-%s" %(control1, control2))

    #Control-case pairs:
    for control in controls:
        sample1 = name2sample[control]
        for case in cases:
            sample2 = name2sample[case]
            controlPairs.append( getAvrSharedSeqs(sample1, sample2, samplingSize, numSamplings) )
            controlNames.append("%s-%s" %(control, case))

    #ttest:
    stime = time.clock()
    tval, pval = ttest_ind(casePairs, controlPairs)
    sys.stderr.write("ttest after gotten the data in %f s\n" %(time.clock() - stime))

    f = open(outfile, 'w')
    f.write("#Group 1: %s\n" %(','.join(caseNames) ) )
    f.write("#Group 2: %s\n" %(','.join(controlNames) ) )
    f.write("#Pval\tTval\tMean1\tStd1\tMean2\tStd2\n")
    f.write("%f\t%f\t%f\t%f\t%f\t%f\n" %(pval, tval, np.mean(casePairs), np.std(casePairs), np.mean(controlPairs), np.std(controlPairs) ))
    f.close()

################## END TTEST ##################

################ OLD MAIN, TO BE DELETED #############

def main():
    starttime = time.clock() 
    parser = initOptions()
    addOptions(parser)
    initPlotOptions(parser)
    options, args = parser.parse_args()
    checkOptions(parser, args, options)
    checkPlotOptions(options, parser)

    stime = time.clock() - starttime
    if options.verbose:
        sys.stderr.write("Reading and checking options in %d seconds.\n" %stime)
        sys.stderr.write("Reading in input files...\n")
    stime = time.clock()

    samples = readfiles(options.indir, options.minCount, 1) 
    sample2aa2v2j = getsample2aa2v2j(samples)
    #if options.nt2aa:
    #    transNt2aa(samples)
    
    if options.verbose:
        sys.stderr.write("Read input files in %f seconds.\n" %(time.clock() - stime))

    if options.sampling > 0:
        print options.sampling
        stime = time.clock()
        samples = sampling(samples, options.sampling, options.samplingUniq)
        if options.verbose:
            sys.stderr.write("Done sampling in %f seconds.\n" %(time.clock() - stime) )

    if options.overlapPlot:
        overlapPlotDir = os.path.join(options.outdir, "overlapPlots")
        system("mkdir -p %s" %(overlapPlotDir))

    if not options.noPairwiseOverlap:
        statdir = os.path.join(options.outdir, "overlapStats")
        system("mkdir -p %s" %statdir)
    
    ###########################################
    #Pairwise overlap statistics and sequences#
    ###########################################
    if options.verbose:
        sys.stderr.write("Pairwise overlap statistics and fasta, if specified:\n")
    pairwiseTime = time.clock()

    for i in xrange( len(samples) - 1 ):
        sample1 = samples[i]
        aa2v2j1 = sample2aa2v2j[ sample1.name ]
        for j in xrange(i+1, len(samples)):
            sample2 = samples[j]
            aa2v2j2 = sample2aa2v2j[ sample2.name ]
            pair = "%s_%s" %(sample1.name, sample2.name)
            if not options.noPairwiseOverlap:
                if options.verbose:
                    sys.stderr.write("Getting pairwise overlap statistics for samples %s and %s.\n" %(sample1.name, sample2.name))
                    stime = time.clock()
                
                outfile = os.path.join(statdir, "%s.txt" %pair)
                reads1, reads2, clones1, clones2, stats1, stats2 = getPairwiseOverlap(sample1.seqs, sample2.seqs, aa2v2j1, aa2v2j2, options.cutoffs, options.mode, options.discrete)
                printPairwiseOverlap(reads1, reads2, clones1, clones2, stats1, stats2, options.cutoffs, outfile)
                
                if options.verbose:
                    sys.stderr.write("Done in %f s.\n" %(time.clock() - stime))
            
            if options.fasta:
                if options.verbose:
                    sys.stderr.write("Print fasta of the pairwise overlap sequences between samples %s and %s\n" %(sample1.name, sample2.name))
                    stime = time.clock()

                fastaDir = os.path.join(options.outdir, "pairwiseSharedSeqs")
                system("mkdir -p %s" %fastaDir)
                printPairwiseOverlapSequences(sample1.name, sample2.name, sample1.seqs, sample2.seqs, aa2v2j1, aa2v2j2, fastaDir, options.cutoffs, options.mode, options.discrete)
                
                if options.verbose:
                    sys.stderr.write("Done in %f s.\n" %(time.clock() - stime))
    if options.verbose:
        sys.stderr.write("Done Pairwise Overlap Statistics and Fasta in %f s.\n" %(time.clock() - pairwiseTime))

    ############################
    #HARLAN ROBINS OVERLAP PLOT 
    #(Y axis = number of shared sequences, X axis = n1j x n2j where nij (i=1,2) is the number of sequences of sample i in the top jth clones. 
    ############################
    if options.verbose:
        sys.stderr.write("Starting overlapPlot if specified.\n")
    plottime = time.clock()
    
    if options.overlapPlot:
        sampleNames = sorted( [sample.name for sample in samples] )
        sam2color = sample2color(sampleNames)
        
        for i in xrange( len(samples) ):
            sample1 = samples[i]
            aa2v2j1 = sample2aa2v2j[sample1.name]
            data = []
            labels = []
            for j in xrange(len(samples) ):
                if i == j:
                    continue
                sample2 = samples[j]
                aa2v2j2 = sample2aa2v2j[sample2.name]

                #Robins overlap plot:
                stime = time.clock()
                xdata, ydata = pairwiseOverlapPlot(sample1.seqs, sample2.seqs, options.overlapPlotRandom)
                data.append( [xdata, ydata] )
                labels.append(sample2.name)

                if options.verbose:
                    sys.stderr.write("Calculated data for overlap plot of pair %s - %s  in %f seconds\n" %(sample1.name, sample2.name, time.clock() - stime))

                #ofile = os.path.join(overlapPlotDir, '%s-%s' %(sample1.name, sample2.name) )
                #drawOverlapPlot([ [xdata, ydata] ], ofile, [ "%s-%s" %(sample1.name, sample2.name)], options)
            ofile = os.path.join(overlapPlotDir, '%s' %(sample1.name) )
            drawOverlapPlot(data, ofile, labels, options, sam2color)

    if options.verbose:
        sys.stderr.write("Done OverlapPlot in %f s.\n" %(time.clock() - plottime))
    ############################
    ##DONE HARLAN ROBINS PLOTS##
    ############################

    #################################
    #====== IF GROUP2SAMPLEs =======#
    #################################
    uniqtime = time.clock()
    if options.verbose:
        sys.stderr.write("Getting uniq and shared sequences if requested...\n")

    #if options.group2samples:
    if options.group2samples:
        group2samples, sample2group = readGroup2samples(options.group2samples)
    if options.uniqAndSharedSeqs:
        seq2sample2count = getSeq2sample2count(samples, options.nt2aa) #key = seq, val = {sample: (count, freq)}
        for sample in samples:
            outdir = os.path.join(options.outdir, "uniqAndSharedSeqs", sample.name)
            system("mkdir -p %s" %outdir)

            stime = time.clock()
            uniqAndSharedSequences(outdir, sample, seq2sample2count, group2samples, sample2group[sample.name], options.nt2aa)
           
            if options.verbose:
                sys.stderr.write("uniq and shared sequences for sample %s in %f seconds.\n" %(sample.name, time.clock() - stime))
    
    if options.verbose:
        sys.stderr.write("Done uniq and shared sequences in %f seconds.\n" %(time.clock() - uniqtime))

    ###############################
    #=== TTEST ===
    ###############################
    if options.ttest or options.overlapPlotAll:
        sep = ','
        cases = readList(options.cases, sep)
        controls = readList(options.controls, sep)
    
    if options.ttest:
        ttesttime = time.clock()
        if options.verbose:
            sys.stderr.write("Ttest...\n")
        ttestOutfile = os.path.join(options.outdir, "ttest.txt")
        ttest( ttestOutfile, samples, cases, controls, options.ttestSamplingSize, options.ttestSamplingNum, options.ttestMode )
        
        if options.verbose:
            sys.stderr.write("Done in %f s\n" %(time.clock() - ttesttime))

    ###############################
    #==== COMBINE OVERLAP PLOT ===#
    ###############################
    if options.overlapPlotAll:
        drawOverlapPlotAll(samples, sample2aa2v2j, cases, controls, options, group2samples, sample2group)
        #drawOverlapPlotAll_3way(samples, sample2aa2v2j, cases, controls, options, group2samples, sample2group)


