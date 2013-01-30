#!/usr/bin/env python2.6

"""
01/16/2012 nknguyen at soe ucsc edu

Input: Input directory with sample fasta files containing the productive clones, with header has the format:
    >SomeID;size=##
        Output directory
Output: Plots of clonesize distribution: 
    a/ Proportion of total reads versus number of clones with size >= the appropriate proportion
    b/ Number of reads versus number of clones with size >= the appropriate reads
    c/ Table with columns as followed:
    TotalReads\tTotalClones\tNumberOfClonesWithSize >= 0.01\tNumberOfClonesWithSize >= 0.05\tNumberOfClonesWithSize >= 0.1
    Rows = Samples
"""

import os, sys, re
from optparse import OptionParser

#import libPlotting as libplot
import immunoseq.lib.immunoseqLib as iseqlib
import matplotlib.pyplot as pyplot
from matplotlib.ticker import *
from matplotlib.font_manager import FontProperties
#from numpy import array

class Sample:
    def __init__(self, name):
        self.name = name
        self.totalCount = 0
        self.totalClones = 0
        self.size2clones = {}
        #self.countArr = [1, 5, 10, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500]
        #self.countArr = xrange(0, 1001, 10)
        #self.countArr = xrange(0, 101, 1)
        self.countArr = [0, 10, 100, 1000, 10000]
        #self.countArr = [0, 100, 1000, 10000, 100000]
        self.percentArr = [0, 0.001, 0.01, 0.1, 1]
        #self.percentArr = [0, 0.1, 1]
        #HACK
        #self.percentArr = [0, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
        #self.percentArr = [0, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
        
        self.clonesPerCount = [ 0 for c in self.countArr ]
        self.readsPerCount = [0 for c in self.countArr] #Percentage of reads per category
        self.clonesPerPercent = [0.0 for p in self.percentArr]
        self.readsPerPercent = [0.0 for p in self.percentArr]
        
        #self.countArr = [1, 10, 100, 1000, 10000, 100000]
        #self.percentArr = [0, 0.001, 0.01, 0.1, 1, 10]
        #self.clonesPerCount = [0, 0, 0, 0, 0, 0]
        #self.readsPerCount = [0, 0, 0, 0, 0, 0] #Percentage of reads per category
        #self.clonesPerPercent = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        #self.readsPerPercent = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    def setUpperLimit(self, maxFreq):
        filteredSize2clones = {}
        total = sum( [ s*c for s,c in self.size2clones.iteritems() ] )
        for size, clones in self.size2clones.iteritems():
            freq = 0
            if total > 0:
                freq = size*100.0/total
            if freq <= maxFreq:
                filteredSize2clones[size] = clones
        self.size2clones = filteredSize2clones
    
    def setTotal(self):
        self.totalCount = sum( [ s*c for s,c in self.size2clones.iteritems() ] )
        self.totalClones = sum( [ c for s,c in self.size2clones.iteritems() ] )

    def setClonesPerCount(self, cumulative):
        for size, count in self.size2clones.iteritems():
            for i, bin in enumerate( self.countArr ):
                if size >= bin:
                    if cumulative or i == len(self.countArr) -1 or size < self.countArr[i+1]:
                        self.clonesPerCount[i] += count
                        self.readsPerCount[i] += count*size
                    
        for i, c in enumerate(self.readsPerCount):
            self.readsPerCount[i] = c*100.0/self.totalCount

    def setClonesPerPercent(self, cumulative):
        for size, count in self.size2clones.iteritems():
            percent = 0.0
            if self.totalCount > 0:
                percent = size*100.0/self.totalCount
            for i, bin in enumerate( self.percentArr ):
                if percent >= bin:
                    if cumulative or i == len(self.percentArr) - 1 or percent < self.percentArr[i+1]:
                        self.clonesPerPercent[i] += count
                        self.readsPerPercent[i] += count*size
        for i, c in enumerate(self.readsPerPercent):
            self.readsPerPercent[i] = c*100.0/self.totalCount

    def avr(self, denom):
        self.totalClones /= denom
        for i in xrange(len(self.countArr)):
            self.clonesPerCount[i] /= denom
        for i in xrange(len(self.percentArr)):
            self.clonesPerPercent[i] /= denom

def readfiles(indir, cumulative):
    samples = []
    for file in os.listdir(indir):
        if not os.path.exists(os.path.join(indir, file)) or (not re.search(".fa", file)):
            continue
        sampleName = os.path.basename(file).rstrip("fa").rstrip(".")
        s = Sample(sampleName)
        #sys.stderr.write("File %s\n" %file)
        f = open(os.path.join(indir, file), 'r')
        for line in f:
            if line[0] == '>':
                items = line.strip().split(';size=') 
                if len(items) < 2:
                    sys.stderr.write("Wrong header format\n")
                    sys.exit(1)
                count = int( items[ len(items) -1 ] )
                
                if count not in s.size2clones:
                    s.size2clones[ count ] = 1
                else:
                    s.size2clones[ count ] += 1
        f.close()
        #s.setUpperLimit(0.005) #HACK
        s.setTotal()
        s.setClonesPerCount(cumulative)
        s.setClonesPerPercent(cumulative)
        samples.append(s)
    return samples

def getAvrSamples(samples, group2samples, cumulative):
    avrSamples = []
    for group, names in group2samples.iteritems():
        groupAvrSample = Sample(group)
        for name in names:
            for s in samples:
                if s.name == name:
                    for size, clones in s.size2clones.iteritems():
                        if size not in groupAvrSample.size2clones:
                            groupAvrSample.size2clones[size] = clones
                        else:
                            groupAvrSample.size2clones[size] += clones
        #average
        #for size in groupAvrSample.size2clones:
        #    groupAvrSample.size2clones[size] /= len(names)
        groupAvrSample.setTotal()
        groupAvrSample.setClonesPerCount(cumulative)
        groupAvrSample.setClonesPerPercent(cumulative)
        
        #average
        groupAvrSample.avr( len(names) )

        avrSamples.append(groupAvrSample)
    return avrSamples

def setUsageDistAxes( fig, numSamples, samplesPerPlot ):
    numaxes = numSamples / samplesPerPlot 
    if numSamples % samplesPerPlot != 0:
        numaxes += 1
    
    axesList = []
    axleft = 0.13
    axright = 0.95
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

#AS color:
def getColors3():
    colors = [ "#00611C", "#0E8C3A", #green
               "#00008B", "#4372AA", #blue
               "#551A8B", "#72587F", # purple
               "#B6316C", "#CD6889", # cranberry
               "#8B4513", "#8B5742", # brown
               "#EE8833", "#EEC591",
               "#00CD00", "#00FF66",
               "#0147FA", "#0198E1",
               "#BC8F8F", "#C76114",
               "#AA6600", "#B5A642",
               "#008837", "#A6DBA0", #green
               "#0571B0", "#92C5DE"] #blue
    return colors

def getColors6():
    #colors = ["#E31A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#1B9E77", "#FFFF33", "#A65628", "#CE1256"] #red, blue, green, purple, orange, green-ish, yellow, brown, pink 
    #GRAY SCALE
    colors = ["#353535", "#959595", "#4DAF4A", "#984EA3", "#FF7F00", "#1B9E77", "#FFFF33", "#A65628", "#CE1256"] #red, blue, green, purple, orange, green-ish, yellow, brown, pink 
    return colors

def getColors6light():
    #colors = ["#FE8E8F", "#A6D7FE", "#B8FEB5", "#F6BDFE", "#FEBF80", "#95FEDF", "#FFFFB3", "#D8885A", "#D7B5D8"] #red, blue, green, purple, orange, green-ish
    #GRAY SCALE
    colors = ["#656565", "#BDBDBD", "#B8FEB5", "#F6BDFE", "#FEBF80", "#95FEDF", "#FFFFB3", "#D8885A", "#D7B5D8"] #red, blue, green, purple, orange, green-ish
    return colors

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
           "f57CD8naive":3, "m35CD8naive":3,
           "adaptAS8Na":0, "adaptAS8Na2":0, "adaptAS8Nb":0, "adaptAS8Nb2":0, "adaptAS8Nc":0,
           "adaptCD8Na":1, "adaptCD8Nb":1, "adaptCD8n2":1,
           "adaptMA8Na":2, "adaptMA8Na2":2, "adaptMA8Nb":2, "adaptMA8Nb2":2, "adaptMA8Nc2":2,
           "SEQ-A1":5, "SEQ-A1b":5, "SEQ-A2":5, "SEQ-B":5,
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
           "f57CD8naive":0, "m35CD8naive":1,
           "adaptAS8Na":0, "adaptAS8Na2":1, "adaptAS8Nb":2, "adaptAS8Nb2":3, "adaptAS8Nc":4,
           "adaptCD8Na":0, "adaptCD8Nb":1, "adaptCD8n2":2,
           "adaptMA8Na":0, "adaptMA8Na2":1, "adaptMA8Nb":2, "adaptMA8Nb2":3, "adaptMA8Nc2":4,
           "SEQ-A1":0, "SEQ-A1b":1, "SEQ-A2":2, "SEQ-B":3,
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

def drawCloneSizeData( axesList, samples, samplesPerPlot, options, isAbs, yaxisPcReads, yaxisPcClones, cumulative ):
    if len( samples ) <= 0:
        return
    colors = getColors6()
    s2c, s2m, s2cLight = sample2color( colors )
    #markers = ['o', '^']
    #markers=['o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', '^', '^']

    #c = -1
    #textsize = 'x-large'
    textsize = 'medium'
    fontP = FontProperties()
    fontP.set_size(textsize)
    linesDict = {}
    labelsDict = {}

    if isAbs:
        xtickLabels = samples[0].countArr
    else:
        xtickLabels = samples[0].percentArr
    
    #get x location
    xdata = range( 0, len(xtickLabels), 1 )
    #offset = 0.01
    offset = 0.0
    maxy = 0
    for i in range( len(axesList) ):
        lines = []
        sampleNames = []
        axes = axesList[i]
        #HACK
        #if not yaxisPcReads and not yaxisPcClones:
        if not yaxisPcReads:
            axes.set_yscale('log')
        #axes.set_yscale('log')
        startIndex = i*samplesPerPlot
        endIndex = min( [startIndex + samplesPerPlot, len(samples)] )
        for j in range( startIndex, endIndex ):
            sample = samples[j]
            #sampleNames.append( "%s" % (sample.name))
            sampleNames.append( "%s" % (iseqlib.properName(sample.name)))
            ydata = sample.clonesPerCount
            if not yaxisPcReads and not isAbs:
                ydata = sample.clonesPerPercent
            elif yaxisPcReads and isAbs:
                ydata = sample.readsPerCount
            elif yaxisPcReads and not isAbs:
                ydata = sample.readsPerPercent
            
            if yaxisPcClones:
                ydata = [ 100*y/sample.totalClones for y in ydata ]
            maxy = max( [maxy, max(ydata)] )

            #if (i*len(axesList) + j) %2 == 0:
            #    c += 1
            #l = axes.plot( xdata, ydata, color=colors[c], marker=markers[(i*len(axesList) + j)%2], markersize=6.0, linestyle='none' )
            #axes.plot( xdata, ydata, color=colors[c], linestyle='-', linewidth=0.01 )
            currxdata = [x + offset*(j - startIndex) for x in xdata]
            if isAbs:
                markersize = 8.0
                #markersize = 10.0
            else:
                markersize = 8.0
                #markersize = 10.0
            if s2m[sample.name] == '*':
                markersize += 2.0
            elif s2m[sample.name] == 's':
                markersize -= 2.0

            print sample.name
            print ydata

            #l = axes.plot( xdata, ydata, color=s2c[sample.name], marker=s2m[sample.name], markeredgecolor=s2c[sample.name], markersize=6.0, linestyle='none' )
            l = axes.plot( currxdata, ydata, color=s2c[sample.name], marker=s2m[sample.name], markeredgecolor=s2c[sample.name], markersize=markersize, linestyle='none' )
            #l = axes.plot( currxdata, ydata, color="#377EB8", marker='o', markeredgecolor=s2c[sample.name], markersize=markersize, linestyle='none' )
            #axes.plot( xdata, ydata, color=s2cLight[sample.name], linestyle='-', linewidth=0.3 )
            axes.plot( currxdata, ydata, color=s2cLight[sample.name], linestyle='-', linewidth=0.5 )
            #axes.plot( currxdata, ydata, color="#377EB8", linestyle='-', linewidth=0.7 )
            lines.append( l )
        
        #if yaxisPcReads or yaxisPcClones:
        #    axes.plot( [-0.5, len(xtickLabels)], [90, 90], linestyle='-', linewidth=0.2, color="#CCCCCC" )
        linesDict[i] = lines
        labelsDict[i] = sampleNames
    
        if cumulative:
            #axes.set_title( 'Cumulative Distribution of Clones', size="xx-large", weight='bold' )
            #axes.set_title( 'Cumulative Distribution of Clones', size="large", weight='bold' )
            axes.set_title( 'A. Cumulative Distribution of Clones', size="large", weight='bold' )
            if yaxisPcReads:
                #axes.set_title( 'Cumulative Distribution of Sequences', size="xx-large", weight='bold')
                axes.set_title( 'B. Cumulative Distribution of Sequences', size="large", weight='bold')
        else:
            axes.set_title( 'Clone Size Distribution', size="xx-large" )

    for i in range( len(axesList) ):
        axes = axesList[ i ]
        #libplot.editSpine( axes )
        iseqlib.editSpine( axes )
        axes.set_xlabel('Clone size (number of sequences)', size = textsize)
        if not isAbs:
            axes.set_xlabel('Clone size (% of total sequences)', size = textsize, weight='bold')
        axes.set_ylabel('Frequency (number of clones)', size=textsize)
        if yaxisPcReads:
            axes.set_ylabel('Frequency (% of total sequences)', size=textsize, weight='bold')
        if yaxisPcClones:
            axes.set_ylabel('Frequency (% of total clones)', size = textsize, weight='bold')
        #Legend
        #legend = pyplot.legend( linesDict[i], labelsDict[i], numpoints=1, prop=fontP, loc="best" )
        #legend = axes.legend( linesDict[i], labelsDict[i], numpoints = 1, "upper right", ncol=3)
        if yaxisPcClones and not yaxisPcReads:
            legend = axes.legend( linesDict[i], labelsDict[i], numpoints = 1, loc="best", ncol=1, prop=fontP)
            #for t in legend.get_texts():
            #    t.set_fontsize(textsize)
            legend._drawFrame = False

        #HACK
        #if isAbs:
        if False:
            xtickdata = xrange(10, 101, 10)
            #grid:
            for x in xtickdata:
                axes.plot([x,x], [0.0001, maxy], color='#3F3F3F', linestyle='-', linewidth=0.1)
            axes.xaxis.set_ticks(xtickdata)
            minx = 0
            maxx = 101
            axes.set_xlim(minx, maxx)
            axes.xaxis.set_ticklabels([str(x) for x in xtickdata] , size='medium')
        else:
            axes.xaxis.set_ticklabels( xtickLabels, size='medium' )
            axes.xaxis.set_ticks([x + 0.5*offset*(samplesPerPlot-1) - 0.5 for x in xdata] )
            minx = -0.5
            maxx = len(xtickLabels) - 0.5
            
            xticks = [x + 0.5*offset*(samplesPerPlot-1) - 0.5 for x in xdata]
            for xi, x in enumerate(xticks):
                if xi == 0:
                    continue
                axes.plot([x,x], [0.0001, maxy], color='#3F3F3F', linestyle='-', linewidth=0.1)
            axes.set_xlim(minx, maxx)

        #END HACK

        #axes.xaxis.set_ticks( xdata )
        #axes.xaxis.set_ticks([x + 0.5*offset*(samplesPerPlot-1) for x in xdata] )
        #axes.set_xlim(-0.5, len(xtickLabels))

        #numTicks = 20
        #yticks = [ float(t)/numTicks for t in range(numTicks +1) ]
        #ytickLabels = []
        #for y in yticks:
        #    ytickLabels.append( "%d" %(y*100) )
        
        if not yaxisPcReads and not yaxisPcClones:
            #yticks = [100, 1000, 5000, 10000, 15000, 20000]
            #ytickLabels = ["100", "1k", "5k", "10k", "15k", "20k"]
            yticks = [1, 10, 100, 1000, 10000, 50000, 100000, 150000]
            #ytickLabels = ["1", "10", "100", "1k", "10k", "50k", "100k", "150k"]
            ytickLabels = ["1", "10", "100", "1k", "10k", "50k", "100k", "150k"]
            axes.yaxis.set_ticklabels( ytickLabels )
            axes.yaxis.set_ticks( yticks )
            for y in yticks:
                if y in [1, 10, 100, 1000, 10000, 100000]:
                #if y in [100, 1000, 5000, 10000, 15000, 20000]:
                    axes.plot([minx, maxx], [y, y], color='#3F3F3F', linestyle='-', linewidth=0.1)
                else:
                    axes.plot([minx, maxx], [y, y], color='#3F3F3F', linestyle='-.', linewidth=0.1)
            axes.set_ylim(0.8, maxy + 15000)
        elif yaxisPcClones:
        #elif yaxisPcClones or yaxisPcReads:
            #yticks = [0.001, 0.01, 0.1, 1, 10, 25, 50, 75, 100]
            yticks = [0.001, 0.01, 0.1, 1, 10, 100]
            #yticks = xrange(0, 101, 20)
            ytickLabels = [str(y) for y in yticks]
            axes.yaxis.set_ticks(yticks)
            axes.yaxis.set_ticklabels( ytickLabels )
            for y in yticks:
                if y in [0.001, 0.01, 0.1, 1, 10, 100]:
                #if y in xrange(0, 101, 20):
                    axes.plot([minx, maxx], [y, y], color='#3F3F3F', linestyle='-', linewidth=0.1)
                else:
                    axes.plot([minx, maxx], [y, y], color='#3F3F3F', linestyle='-.', linewidth=0.1)
            axes.set_ylim(0.0005, maxy + 25)
            #axes.set_ylim(0.001, maxy + 25)
        else:
            yticks = xrange(0, 101, 20)
            ytickLabels = [str(y) for y in yticks]
            axes.yaxis.set_ticks(yticks)
            axes.yaxis.set_ticklabels( ytickLabels )
            for y in yticks:
                axes.plot([minx, maxx], [y, y], color='#3F3F3F', linestyle='-', linewidth=0.1)
            axes.set_ylim(-1, maxy + 1)
        #axes.set_ylim(-1, maxy + 1)
        #if isAbs: #HACK
        #    axes.set_ylim(-0.1, 15)

        #axes.set_ylim(-0.1, 50) #HACK
        for label in axes.get_xticklabels():
            #label.set_fontsize( 'large' )
            label.set_fontsize( 'small' )
            #label.set_rotation( 45 )

        for label in axes.get_yticklabels():
            #label.set_fontsize( 'large' )
            label.set_fontsize( 'small' )
        
        #axes.yaxis.grid(b=True, color="#BCBCBC", linestyle='-', linewidth=0.005)
        #axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    return 

def drawCloneSizeDist( samples, options, isAbs, yaxisPcReads, yaxisPcClones, cumulative ):
    options.out = os.path.join( options.outdir, "cloneSizeDist" )
    if cumulative:
        options.out += "-cumulative"
    if not isAbs:
        options.out = os.path.join( options.outdir, "cloneSizeDist-Rel" )
    if yaxisPcReads:
        options.out += "-pcReads"
    if yaxisPcClones:
        options.out += "-pcClones"
    #fig, pdf = libplot.initImage( 10.0, 8.0, options )
    fig, pdf = iseqlib.initImage( 10.0, 8.0, options )
    axesList = setUsageDistAxes( fig, len(samples), options.samplesPerPlot )
    drawCloneSizeData( axesList, samples, options.samplesPerPlot, options, isAbs, yaxisPcReads, yaxisPcClones, cumulative )
    #libplot.writeImage( fig, pdf, options )
    iseqlib.writeImage( fig, pdf, options )

#============ Top Clones Versus Read Proportion ===============
def getClonesVsReads( sample ):
    #HACK:
    xdata = range(1, 101, 1) 
    xlabels = [str(x)  for x in xdata ] 
    
    #xdata = [1, 5, 10, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 500, 1000]
    #xlabels = ['1', '5', '10', '25', '50', '75', '100', '125', '150', '175', '200', '225', '250', '275', '300', '350', '400', '500', '1000' ]
    #xdata = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    #xlabels = [str(x) for x in xdata]
    readPercentage = {}
    for x in xdata:
        readPercentage[x] = 0
    #xlabels = ['1', '5', '10', '50', '100', '500', '1000', '5000', 'all']
    #xdata = [1, 5, 10, 50, 100, 500, 1000, 5000, 10000]
    #readPercentage = {1:0, 5:0, 10:0, 50:0, 100:0, 500:0, 1000:0, 5000:0, 'all':0}
    cumulClones = 0
    for i,size in enumerate( sorted( sample.size2clones.keys(), reverse=True ) ):
        #HACK:
        #if i == 0:
        #    continue
        #END HACK
        numClones = sample.size2clones[size]
        cumulClones += numClones
        for k in readPercentage:
            if k == 'all':
                readPercentage[k] += size*numClones
            elif cumulClones <= k:
                readPercentage[k] += size*numClones
            elif cumulClones - numClones < k:
                readPercentage[k] += size*( k - (cumulClones -numClones) )
    ydata = []
    for k in sorted(readPercentage.keys()):
        ydata.append( readPercentage[k]*100.0/sample.totalCount )
    return xlabels, xdata, ydata 

def getClonesVsReadsRel( sample ):
    readPercentage = {1:0, 5:0, 10:0, 25:0, 50:0, 75:0, 100:0}
    cumulClones = 0
    for size in sorted( sample.size2clones.keys(), reverse=True ):
        numClones = sample.size2clones[size]
        cumulClones += numClones
        for k in readPercentage:
            clones = float(k)*sample.totalClones/100.0
            if k == 100:
                readPercentage[k] += size*numClones
            elif cumulClones <= clones:
                readPercentage[k] += size*numClones
            elif cumulClones - numClones < clones:
                readPercentage[k] += size*( k - (cumulClones -numClones) )
    ydata = []
    xdata = []
    for k in sorted(readPercentage.keys()):
        xdata.append(k)
        ydata.append( readPercentage[k]*100.0/sample.totalCount )
    xlabels = [ str(x) for x in xdata]
    return xlabels, xdata, ydata 


def drawCloneVsReadData( axesList, samples, samplesPerPlot, options, isAbs ):
    if len( samples ) <= 0:
        return
    colors = getColors6()
    s2c, s2m, s2cLight = sample2color( colors )
    #textsize = 'x-large'
    textsize = 'medium'
    fontP = FontProperties()
    fontP.set_size(textsize)
    linesDict = {}
    labelsDict = {}
    xticklabels = [] 
    
    #get x location
    #xdata = range( 0, len(xtickLabels), 1 )
    maxy = 0
    #markersize = 10.0
    markersize = 8.0
    name2line = {} 
    for i in range( len(axesList) ):
        lines = []
        sampleNames = []
        axes = axesList[i]
        #HACK
        #if isAbs:
        #    axes.set_xscale('log')
        startIndex = i*samplesPerPlot
        endIndex = min( [startIndex + samplesPerPlot, len(samples)] )
        for j in range( startIndex, endIndex ):
            sample = samples[j]
            if isAbs:
                xtickLabels, xdata, ydata = getClonesVsReads( sample )
            else:
                xtickLabels, xdata, ydata = getClonesVsReadsRel( sample )
           
            #HACK discrete:
            discreteydata = [ydata[0]]
            if not options.cumulative:
                for k in xrange(1, len(ydata)):
                    discreteydata.append( ydata[k] - ydata[k-1] )
                ydata = discreteydata


            #sampleNames.append( "%s-%d" % (sample.name, sample.totalCount))
            sampleNames.append( "%s" % (sample.name))
            maxy = max( [maxy, max(ydata)] )

            l = axes.plot( xdata, ydata, color=s2c[sample.name], marker=s2m[sample.name], markeredgecolor=s2c[sample.name], markersize=markersize, linestyle='none' )
            axes.plot( xdata, ydata, color=s2cLight[sample.name], linestyle='-', linewidth=0.2 )
            lines.append( l )
            name2line[sample.name] = l 
        
        #if yaxisPcReads or yaxisPcClones:
        #    axes.plot( [-0.5, len(xtickLabels)], [90, 90], linestyle='-', linewidth=0.3, color="#CCCCCC" )
        linesDict[i] = lines
        labelsDict[i] = sampleNames
        if options.cumulative:
            #axes.set_title( 'Cumulative Distribution of 50 Largest Clones', size='xx-large', weight='bold' )
            axes.set_title( 'C. Cumulative Distribution of 50 Largest Clones', size='large', weight='bold' )
        else:
            axes.set_title( 'Distribution of 50 Largest Clones', size='xx-large' )

    if not options.cumulative:
        axes.set_yscale('log')

    for i in range( len(axesList) ):
        axes = axesList[ i ]
        #libplot.editSpine( axes )
        iseqlib.editSpine( axes )
        axes.set_xlabel('Number of top clones', size = textsize, weight='bold')
        if not options.cumulative:
            axes.set_xlabel('Clone rank', size = textsize, weight='bold')

        if not isAbs:
            axes.set_xlabel('Percentage of total clones', size = textsize, weight='bold')
        axes.set_ylabel('Frequency (% of total sequences)', size = textsize, weight= 'bold')
        #legend = pyplot.legend( linesDict[i], labelsDict[i], numpoints=1, prop=fontP, loc="best" )
        
        ##HACK
        #labelorder = ['as15D', 'as20D', 'asBD', 'asBDdraw2', 'as8D', 'as16D', 'as1D', 'as11D']
        #currlabels = []
        #currlines = []
        #for label in labelorder:
        #    if label in labelsDict[i]:
        #        currlabels.append( iseqlib.properName(label) )
        #        currlines.append(name2line[label])
        ##legend = pyplot.legend( currlines, currlabels, numpoints=1, prop=fontP, loc="upper right" )
        #legend = pyplot.legend( currlines, currlabels, numpoints=1, prop=fontP, loc="best" )
        ##END HACK
        #legend._drawFrame = False

        axes.xaxis.set_ticklabels( xtickLabels )
        axes.xaxis.set_ticks( xdata )
        
        #axes.set_xlim(-0.5, len(xtickLabels))
        
        for label in axes.get_xticklabels():
            #label.set_fontsize( 'large' )
            label.set_fontsize( 'small' )
            #label.set_rotation( 90 )
        for label in axes.get_yticklabels():
            #label.set_fontsize( 'large' )
            label.set_fontsize( 'small' )
        
        #HACK: yaxis tick labels:
        #ytickdata = [0.1, 1, 2.5, 5, 7.5, 10, 15, 20]
        #for y in ytickdata:
        #    if y == 20:
        #        continue
        #    if y in [0.1, 1, 10, 100]:
        #        #axes.plot([0,50.5], [y, y], color='#838383', linestyle='-', linewidth=0.005)
        #        axes.plot([0,50.5], [y, y], color='#BCBCBC', linestyle='-', linewidth=0.005)
        #    else:
        #        axes.plot([0,50.5], [y, y], color='#BCBCBC', linestyle='-.', linewidth=0.005)
        #axes.yaxis.set_ticks(ytickdata)
        #axes.set_ylim(0, 21)
        #axes.yaxis.set_ticklabels([str(y) for y in ytickdata] , size='medium')
        
        #x ticks:
        xticks = [1]
        xticks.extend( xrange(5, 51, 5) )
        axes.xaxis.set_ticks( xticks )
        axes.xaxis.set_ticklabels( [str(x) for x in xticks] , size = 'large')
        #x tick lines:
        for x in xdata:
            if x in xrange(0, 51, 5):
                axes.plot([x,x], [0.02, maxy+5], color='#3F3F3F', linestyle='-', linewidth=0.05)
            #else:
            #    axes.plot([x,x], [0.02, maxy+5], color='#BCBCBC', linestyle='-.', linewidth=0.005)
        
        #END HACK
        if isAbs:
            #HACK
            axes.set_xlim(0, 50.5)
            #axes.set_xlim(0, 525)
            #axes.set_xlim(0, 10500) #before hack
        else:
            axes.set_xlim(0, 101)
            #axes.set_xlim(-1, 51)

        #axes.set_ylim(0.02, maxy + 5)
        axes.set_ylim(0, maxy + 2)
        #HACK
        #axes.set_ylim(-1, maxy + 2)

        axes.yaxis.grid(b=True, color="#3F3F3F", linestyle='-', linewidth=0.05)
        #axes.xaxis.grid(b=True, color="#3F3F3F", linestyle='-', linewidth=0.05)
    
    return 

def drawCloneVsRead( samples, options, isAbs ):
    options.out = os.path.join( options.outdir, "cloneVsRead" )
    if not isAbs:
        options.out = os.path.join( options.outdir, "cloneVsRead-Rel" )
    #fig, pdf = libplot.initImage( 10.0, 8.0, options )
    fig, pdf = iseqlib.initImage( 10.0, 8.0, options )
    axesList = setUsageDistAxes( fig, len(samples), options.samplesPerPlot )
    drawCloneVsReadData( axesList, samples, options.samplesPerPlot, options, isAbs )
    #libplot.writeImage( fig, pdf, options )
    iseqlib.writeImage( fig, pdf, options )

def drawCombine( samples, options ):
    options.out = os.path.join(options.outdir, "combine")
    fig, pdf = iseqlib.initImage( 8.0, 12.0, options )
    #fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1) 
    fig.subplots_adjust(hspace=.5)

    #drawCloneSizeDist( samples, options, isAbs, yaxisPcReads, yaxisPcClones, cumulative )
    ax1 = fig.add_subplot(311)
    drawCloneSizeData( [ax1], samples, options.samplesPerPlot, options, False, False, True, True )
    
    ax2 = fig.add_subplot(312)
    drawCloneSizeData( [ax2], samples, options.samplesPerPlot, options, False, True, False, True )
    
    ax3 = fig.add_subplot(313)
    #drawCloneVsRead( samples, options, True )
    drawCloneVsReadData( [ax3], samples, options.samplesPerPlot, options, True )

    #Write to output file  
    #import matplotlib.backends.backend_pdf as pltBack
    #pdf = None
    #if options.outFormat == 'pdf' or options.outFormat == 'all':
    #    pdf = pltBack.PdfPages( options.out + '.pdf' )
    iseqlib.writeImage( fig, pdf, options )

def readGroup2Samples(file):
    group2samples = {}
    f = open(file, 'r')
    for line in f:
        items = line.strip().split()
        group2samples[items[0]] = items[1].split(',')
    f.close()
    return group2samples

def checkOptions( args, options, parser ):
    if len(args) < 1:
        parser.error("No input directory was specified!\n")
    if not os.path.isdir( args[0] ):
        parser.error("Input directory '%s' does not exist.\n" % args[0])
    if not os.path.isdir( options.outdir ):
        parser.error("Output directory '%s' does not exist.\n" %options.outdir)

def initOptions( parser ):
    parser.add_option('-o', '--outdir', dest='outdir', default='.', help="Output directory. Default = current direcotory")
    parser.add_option('--samplesPerPlot', dest='samplesPerPlot', type='int', default=25)
    parser.add_option('-c', '--cumulative', dest='cumulative', action='store_true', default=False, help="If specified, will plot the cumululative distribution instead of discrete.Default=%default")
    parser.add_option('-g', '--group', dest='group', help='Group2samples file (Format: group<space>commaSeparatedListOfSamples. If specified, only draw average of each group.')
    #parser.add_option('--isAbs', dest='isAbs', action='store_true', default=False, help='If specified, cloneSize is displayed as read counts, otherwise, as percentage of total read counts.')
    #parser.add_option('--yaxisPcReads', dest='yaxisPcReads', action='store_true', default=False, help='If specified, yaxis is Percentage of Total Reads. Default=Number of Clones')
    #parser.add_option('--yaxisPcClones', dest='yaxisPcClones', action='store_true', default=False, help='If specified, yaxis is Percentage of Total Clones. Default=Number of Clones')
    
def main():
    usage = ('Usage: %prog [options] inputDirectory\n')

    parser = OptionParser(usage = usage)
    initOptions( parser )
    #libplot.initOptions( parser )
    iseqlib.initPlotOptions( parser )
    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    #libplot.checkOptions( options, parser )
    iseqlib.checkPlotOptions( options, parser )

    samples = readfiles( args[0], options.cumulative )
    if options.group:
        group2samples = readGroup2Samples(options.group) 
        samples = getAvrSamples(samples, group2samples, options.cumulative)
    
    #drawCloneSizeDist( samples, options, True, False, False, options.cumulative )
    #drawCloneSizeDist( samples, options, False, False, False, options.cumulative )
    #drawCloneSizeDist( samples, options, True, False, True, options.cumulative )
    #drawCloneSizeDist( samples, options, False, False, True, options.cumulative )
    #drawCloneSizeDist( samples, options, True, True, False, options.cumulative )
    #drawCloneSizeDist( samples, options, False, True, False, options.cumulative )
    
    #drawCloneSizeDist( samples, options, options.isAbs, options.yaxisPcReads, options.yaxisPcClones, options.cumulative )
    
    #if not options.group:
    #    drawCloneVsRead( samples, options, True )
    #    drawCloneVsRead( samples, options, False )

    drawCombine( samples, options )

if __name__ == "__main__":
    main()
