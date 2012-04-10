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
        self.countArr = [1, 5, 10, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500]
        #self.percentArr = [0, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
        #HACK
        self.percentArr = [0, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
        
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
        sampleName = os.path.basename(file).rstrip(".fa")
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
    axleft = 0.1
    axright = 0.97
    axwidth = axright - axleft
    axbottom = 0.13
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
    colors = ["#E31A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#1B9E77", "#FFFF33", "#A65628", "#CE1256"] #red, blue, green, purple, orange, green-ish, yellow, brown, pink 
    return colors

def getColors6light():
    colors = ["#FE8E8F", "#A6D7FE", "#B8FEB5", "#F6BDFE", "#FEBF80", "#95FEDF", "#FFFFB3", "#D8885A", "#D7B5D8"] #red, blue, green, purple, orange, green-ish
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
           "as8D":0, "as10R":0, "as11R":0, "as12R":0, "as13R":0, "as15D":0, "as1DR":0,
           "diffGroupShared":0, "sameGroupShared":1, "shared":2, "uniq":3
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
    
    markers = ['o', '*', 's', 'd', '^', 'p', 'v']
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
           "as8D":1, "as10R":2, "as11R":3, "as12R":4, "as13R":5, "as15D":6, "as1DR":0,
           "diffGroupShared":0, "sameGroupShared":1, "shared":2, "uniq":3
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
        #    axes.set_yscale('log')
        startIndex = i*samplesPerPlot
        endIndex = min( [startIndex + samplesPerPlot, len(samples)] )
        for j in range( startIndex, endIndex ):
            sample = samples[j]
            sampleNames.append( "%s-%d" % (sample.name, sample.totalCount))
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
            markersize = 8.0
            if s2m[sample.name] == '*':
                markersize = 10.0
            elif s2m[sample.name] == 's':
                markersize = 6.0
            #l = axes.plot( xdata, ydata, color=s2c[sample.name], marker=s2m[sample.name], markeredgecolor=s2c[sample.name], markersize=6.0, linestyle='none' )
            l = axes.plot( currxdata, ydata, color=s2c[sample.name], marker=s2m[sample.name], markeredgecolor=s2c[sample.name], markersize=markersize, linestyle='none' )
            #axes.plot( xdata, ydata, color=s2cLight[sample.name], linestyle='-', linewidth=0.3 )
            axes.plot( currxdata, ydata, color=s2cLight[sample.name], linestyle='-', linewidth=0.2 )
            lines.append( l )
        
        if yaxisPcReads or yaxisPcClones:
            axes.plot( [-0.5, len(xtickLabels)], [90, 90], linestyle='-', linewidth=0.2, color="#CCCCCC" )
        linesDict[i] = lines
        labelsDict[i] = sampleNames
    
        if cumulative:
            axes.set_title( 'Cumulative Clone Size Distribution', size="xx-large" )
        else:
            axes.set_title( 'Clone Size Distribution', size="xx-large" )

    for i in range( len(axesList) ):
        axes = axesList[ i ]
        #libplot.editSpine( axes )
        iseqlib.editSpine( axes )
        axes.set_xlabel('Clone size as number of reads', size = textsize)
        if not isAbs:
            axes.set_xlabel('Clone size as percentage of total reads', size = textsize)
        axes.set_ylabel('Number of clones', size=textsize)
        if yaxisPcReads:
            axes.set_ylabel('Percentage of total reads', size=textsize)
        if yaxisPcClones:
            axes.set_ylabel('Percentage of total clones', size = textsize)
        #Legend
        #legend = pyplot.legend( linesDict[i], labelsDict[i], numpoints=1, prop=fontP, loc="best" )
        #legend = axes.legend( linesDict[i], labelsDict[i], numpoints = 1, "upper right", ncol=3)
        legend = axes.legend( linesDict[i], labelsDict[i], numpoints = 1, loc="best", ncol=1, prop=fontP)
        #for t in legend.get_texts():
        #    t.set_fontsize(textsize)
        legend._drawFrame = False

        axes.xaxis.set_ticklabels( xtickLabels )
        #axes.xaxis.set_ticks( xdata )
        #axes.xaxis.set_ticks([x + 0.5*offset*(samplesPerPlot-1) for x in xdata] )
        #axes.set_xlim(-0.5, len(xtickLabels))
        axes.xaxis.set_ticks([x + 0.5*offset*(samplesPerPlot-1) - 0.5 for x in xdata] )
        axes.set_xlim(-0.5, len(xtickLabels) - 0.5)

        #numTicks = 20
        #yticks = [ float(t)/numTicks for t in range(numTicks +1) ]
        #ytickLabels = []
        #for y in yticks:
        #    ytickLabels.append( "%d" %(y*100) )
        
        if not yaxisPcReads and not yaxisPcClones:
            yticks = [1, 10, 100, 1000,5000, 10000, 50000, 100000, 150000]
            ytickLabels = ["1", "10", "100", "1,000", "5,000", "10,000", "50,000", "100,000", "150,000"]
            axes.yaxis.set_ticklabels( ytickLabels )
            axes.yaxis.set_ticks( yticks )
        #axes.set_ylim(-0.1, maxy + 0.1)
        axes.set_ylim(-0.1, 5) #HACK
        for label in axes.get_xticklabels():
            label.set_fontsize( textsize )
            label.set_rotation( 45 )
        for label in axes.get_yticklabels():
            label.set_fontsize( textsize )
        axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
        axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
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
    xdata = range(1, 525) 
    xlabels = [str(x)  for x in xdata ] 
    #xlabels = ['1', '5', '10', '25', '50', '75', '100', '125', '150', '175', '200', '225', '250', '275', '300', '350', '400', '500' ]
    #xdata = [1, 5, 10, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 500]
    readPercentage = {}
    for x in xdata:
        readPercentage[x] = 0
    #xlabels = ['1', '5', '10', '50', '100', '500', '1000', '5000', 'all']
    #xdata = [1, 5, 10, 50, 100, 500, 1000, 5000, 10000]
    #readPercentage = {1:0, 5:0, 10:0, 50:0, 100:0, 500:0, 1000:0, 5000:0, 'all':0}
    cumulClones = 0
    for i,size in enumerate( sorted( sample.size2clones.keys(), reverse=True ) ):
        #HACK:
        if i == 0:
            continue
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
    textsize = 'medium'
    fontP = FontProperties()
    fontP.set_size(textsize)
    linesDict = {}
    labelsDict = {}
    xticklabels = [] 
    
    #get x location
    #xdata = range( 0, len(xtickLabels), 1 )
    maxy = 0
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
            
            #sampleNames.append( "%s-%d" % (sample.name, sample.totalCount))
            sampleNames.append( "%s" % (sample.name))
            maxy = max( [maxy, max(ydata)] )

            l = axes.plot( xdata, ydata, color=s2c[sample.name], marker=s2m[sample.name], markeredgecolor=s2c[sample.name], markersize=8.0, linestyle='none' )
            axes.plot( xdata, ydata, color=s2cLight[sample.name], linestyle='-', linewidth=0.2 )
            lines.append( l )
        
        #if yaxisPcReads or yaxisPcClones:
        #    axes.plot( [-0.5, len(xtickLabels)], [90, 90], linestyle='-', linewidth=0.3, color="#CCCCCC" )
        linesDict[i] = lines
        labelsDict[i] = sampleNames
        axes.set_title( 'Clonotypes versus Reads' )

    for i in range( len(axesList) ):
        axes = axesList[ i ]
        #libplot.editSpine( axes )
        iseqlib.editSpine( axes )
        axes.set_xlabel('Number of clonotypes', size = textsize)
        if not isAbs:
            axes.set_xlabel('Percentage of total clonotypes', size = textsize)
        axes.set_ylabel('Percentage of total reads', size=textsize)
        legend = pyplot.legend( linesDict[i], labelsDict[i], numpoints=1, prop=fontP, loc="best" )
        legend._drawFrame = False

        axes.xaxis.set_ticklabels( xtickLabels )
        axes.xaxis.set_ticks( xdata )
        #axes.set_xlim(-0.5, len(xtickLabels))
        if isAbs:
            #HACK
            axes.set_xlim(0, 75)
            #axes.set_xlim(0, 525)
            #axes.set_xlim(0, 10500) #before hack
        else:
            axes.set_xlim(0, 101)

        #axes.set_ylim(-0.05, maxy + 0.05)
        #HACK
        axes.set_ylim(0, maxy + 0.05)
        for label in axes.get_xticklabels():
            label.set_fontsize( textsize )
            label.set_rotation( 45 )
        for label in axes.get_yticklabels():
            label.set_fontsize( textsize )
        axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
        axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
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
    
    drawCloneSizeDist( samples, options, True, False, False, options.cumulative )
    drawCloneSizeDist( samples, options, False, False, False, options.cumulative )
    drawCloneSizeDist( samples, options, True, False, True, options.cumulative )
    drawCloneSizeDist( samples, options, False, False, True, options.cumulative )
    drawCloneSizeDist( samples, options, True, True, False, options.cumulative )
    drawCloneSizeDist( samples, options, False, True, False, options.cumulative )
    
    #drawCloneSizeDist( samples, options, options.isAbs, options.yaxisPcReads, options.yaxisPcClones, options.cumulative )
    
    if not options.group:
        drawCloneVsRead( samples, options, True )
        drawCloneVsRead( samples, options, False )

if __name__ == "__main__":
    main()
