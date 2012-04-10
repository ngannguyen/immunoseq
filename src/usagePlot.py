#!/usr/bin/env python2.6

"""
08/30/2011 nknguyen at soe ucsc edu

Input: location of vjUsage matrices where rows = Vs and columns= Js.
Each cell X_ij contains the (clonotype/read) count that have the combination V_i and J_j.
Each file corresponds to one sample.

Output: j usage plot, v usage plot (x = genes, y = usage,  different lines = samples)
        vjUsage heatmap for each sample
        vjUsage heatmaps for all samples in one file
"""

import os, sys, re
from optparse import OptionParser

#import libPlotting as libplot
import matplotlib.pyplot as pyplot
from matplotlib.ticker import *
from matplotlib.font_manager import FontProperties
#from numpy import array
import immunoseq.lib.immunoseqLib as immunoseqLib
from scipy.stats.stats import pearsonr
from scipy.stats import ttest_ind
import numpy as np
import mdp
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#from mpl_toolkits.mplot3d import Axes3D

class Sample:
    def __init__(self, name):
        self.name = name
        self.vjusage = [] #list of list (2D matrix: rows = Vs, Columns = Js)
        self.vcounts = []
        self.jcounts = []
        self.totalCount = 0.0

    def setUsage(self, matrix, rownames, colnames):
        self.vjusage = matrix
        self.vgenes = rownames
        self.jgenes = colnames
        if len(self.vgenes) != len(self.vjusage):
            sys.stderr.write("Row labels vector has different length from the number of rows.\n")
            sys.exit(1)
        for row in self.vjusage:
            if len(row) != len(self.jgenes):
                sys.stderr.write("Column labels vector has different length from the number of columns.\n")
                sys.stderr.write("Column labels: %s\n" %('\t'.join(self.jgenes)))
                sys.stderr.write("Sample: %s\n" %self.name)
                sys.exit(1)
    
    def setVcounts( self ):
        self.vcounts = [ sum(v) for v in self.vjusage ]
    
    def setJcounts( self ):
        self.jcounts = [ sum([v[j] for v in self.vjusage]) for j in range( len(self.jgenes) ) ]

    def setTotalCount( self ):
        self.totalCount = sum( [sum(v) for v in self.vjusage] )

    def getVfreqs( self ):
        if self.totalCount == 0:
            return []
        return [ c/self.totalCount for c in self.vcounts ]
    
    def getVfreqs2( self, genes, abs ):
        if self.totalCount == 0:
            return []
        freqs = []
        for g in genes:
            if g not in self.vgenes:
                raise KeyError("Gene %s should be in vgenes, but were not found of %s. Vgenes: %s\n" %(g, self.name, ",".join(self.vgenes)))
            for i, currgene in enumerate(self.vgenes):
                if g == currgene:
                    if abs:
                        freqs.append(self.vcounts[i])
                    else:
                        freqs.append( self.vcounts[i]/self.totalCount )
        return freqs

    def getJfreqs( self ):
        if self.totalCount == 0:
            return []
        return [ c/self.totalCount for c in self.jcounts ]

    def getJfreqs2( self, genes, abs ):
        if self.totalCount == 0:
            return []
        freqs = []
        for g in genes:
            if g not in self.jgenes:
                raise KeyError("Gene %s should be in jgenes, but were not found of %s. Vgenes: %s\n" %(g, self.name, ",".join(self.jgenes)))
            for i, currgene in enumerate(self.jgenes):
                if g == currgene:
                    if abs:
                        freqs.append(self.jcounts[i])
                    else:
                        freqs.append( self.jcounts[i]/self.totalCount )
        return freqs
    
    def normalizeVJusage( self ):
        m = []
        for r in self.vjusage:
            row = [ float(c)/self.totalCount for c in r ]
            m.append( row )
        return m

def getGeneIndices(genes, subgenes):
    indices = []
    for g in subgenes:
        for i, gene in enumerate(genes):
            if g == gene:
                indices.append(i)
                break
    return indices

def intersectVJ(samples, type2intersectGenes):
    ivs = type2intersectGenes['v']
    ijs = type2intersectGenes['j']
    for sample in samples:
        sample.intersectVJusage = intersectVJsample( sample.vjusage, sample.vgenes, sample.jgenes, ivs, ijs)
        sample.normIntersectVJusage = normalizeVJusage( sample.intersectVJusage )
    return

def intersectVJsample(vjusage, vs, js, ivs, ijs):
    #vs: current V genes (or vjusage rownames), js: current J genes (or vjusage colnames), ivs: itersect v genes, ijs: intersect j genes
    vindices = getGeneIndices(vs, ivs)
    jindices = getGeneIndices(js, ijs)
    m = []
    for vindex in vindices: #rows:
        row = vjusage[vindex]
        newrow = []
        for jindex in jindices:
            newrow.append( row[jindex] )
        m.append(newrow)
    return m

def getTotalCount( vjusage ):
    return sum( [sum(v) for v in vjusage] )

def normalizeVJusage(vjusage):
    m = []
    totalCount = getTotalCount(vjusage)
    if totalCount == 0:
        return vjusage

    for r in vjusage:
        row = [float(c)/totalCount for c in r]
        m.append(row)
    return m

def setUsageDistAxes( fig, numSamples, samplesPerPlot ):
    numaxes = numSamples / samplesPerPlot 
    if numSamples % samplesPerPlot != 0:
        numaxes += 1
    
    axesList = []
    axleft = 0.1
    axright = 0.95
    axwidth = axright - axleft
    axbottom = 0.15
    axtop = 0.9
    axheight = axtop - axbottom
    margin = 0.1

    h = ( axheight - ( margin * (numaxes - 1) ) )/float( numaxes )

    bottom = axtop - h
    for i in range( numaxes ):#add axes from top down
        axesList.append( fig.add_axes( [axleft, bottom, axwidth, h] ) )
        bottom = bottom - (margin + h)
    
    return axesList

def getGeneNumber(genename):
    items = re.split('\D+', genename)
    nums = []
    for i in items:
        if i != '':
            nums.append( int(i) )
    return nums
    
def getIntersect(samples):
    type2genes = {} 
    types = ['v', 'j']
    for type in types:
        genes = getIntersectGenes(samples, type)
        genes = sorted( genes, key=lambda g: getGeneNumber(g) )
        type2genes[type] = genes
        #type2genes[type] = sorted( getIntersectGenes(samples, type) )
    return type2genes

def getIntersectGenes(samples, genetype):
    if genetype not in ['v', 'j']:
        raise ValueError("genetype must be either 'v' or 'j'\n")
    if len(samples) == 0:
        raise ValueError("len(samples) = 0! (zero sample!)")
    genes = samples[0].vgenes
    filterGenes = []
    if genetype == 'j':
        genes = samples[0].jgenes

    for i in xrange(1, len(samples)):#gene must be in all other samples, otherwise remove from list
        s = samples[i]
        for j,g in enumerate(genes):
            if (genetype == 'j' and g not in s.jgenes) or (genetype == 'v' and g not in s.vgenes):
                if g not in filterGenes:
                    filterGenes.append(g)
    intersectGenes = []
    for g in genes:
        if g not in filterGenes:
            intersectGenes.append(g)
    return intersectGenes

#def sortGenes(genes, sample, genetype):
def sortGenes(genes, samples, genetype):
    if genetype not in ['v', 'j']:
        raise KeyError("genetype must be either 'v' or 'j'\n")
    sortedgenes = []
    gc = []
    if genetype == 'v':
        for g in genes:
            c = 0
            for sample in samples:
                for i, sg in enumerate(sample.vgenes):
                    if sg == g:
                        c += sample.vcounts[i]
                        break
            gc.append( (g, c) )
    else:
        for g in genes:
            c = 0
            for sample in samples:
                for i, sg in enumerate(sample.jgenes):
                    if sg == g:
                        c += sample.jcounts[i]
                        break
            gc.append( (g, c) )
    #if genetype == 'v':
    #    for i,g in enumerate(sample.vgenes):
    #        gc.append( (g, sample.vcounts[i]) )
    #else:
    #    for i,g in enumerate(sample.jgenes):
    #        gc.append( (g, sample.jcounts[i]) )
    sortedgc = sorted(gc, key=lambda item:item[1], reverse=True)
    sortedgenes = [item[0] for item in sortedgc]
    #for item in sortedgc:
    #    if item[0] in genes:
    #        sortedgenes.append( item[0] )
    return sortedgenes

def getStdData(std, sample, genetype, genes):
    stdData = getGeneUsage(std, genetype, genes, True)
    ydata = getGeneUsage(sample, genetype, genes, True)
    total = sum(ydata)
    if len(stdData) == 0 or total == 0:
        return [], []
    higherYdata = [ float(y + stdData[i])/total for i, y in enumerate(ydata) ]
    lowerYdata = [ float(y - stdData[i])/total for i, y in enumerate(ydata) ]

    return higherYdata, lowerYdata

def getGeneUsage( sample, genetype, genes, abs ):
    if genetype == "v":
        return sample.getVfreqs2(genes, abs)
    elif genetype == "j":
        return sample.getJfreqs2(genes, abs)
    else:
        sys.stderr.write( "Genetype '%s' is unrecognized. Genetype can only be 'v' or 'j'.\n" %genetype )
        sys.exit(1)

#def getGeneUsage( sample, genetype ):
#    if genetype == "v":
#        return sample.getVfreqs()
#    elif genetype == "j":
#        return sample.getJfreqs()
#    else:
#        sys.stderr.write( "Genetype '%s' is unrecognized. Genetype can only be 'v' or 'j'.\n" %genetype )
#        sys.exit(1)

#AS color:
def getColors3():
    colors = [ "#00611C", "#0E8C3A", #green
               "#00CD00", "#00FF66", # bright green
               "#551A8B", "#72587F", # purple
               "#B6316C", "#CD6889", # cranberry
               "#FFFF33", "#FFFFB3", #yellow
               "#EE8833", "#EEC591", # orange
               "#00008B", "#4372AA", #blue
               "#0147FA", "#0198E1", # bright blue
               "#984EA3", "#F6BDFE", #violet
               "#8B4513", "#8B5742", # brown
               "#BC8F8F", "#C76114", #
               "#AA6600", "#B5A642"] #
    return colors

def drawUsageData( axesList, samples, stds, options, genetype, intersectGenes ):
    if len( samples ) <= 0:
        return
    samplesPerPlot = options.samplesPerPlot 
    
    #HACK COLORS
    colors2 = [ "#E31A1C", "#FE8E8F", "#E31A1C", "#FE8E8F","#E31A1C", "#FE8E8F","#E31A1C", "#FE8E8F", "#984EA3", "#F6BDFE","#377EB8", "#A6D7FE", "#377EB8", "#A6D7FE","#377EB8", "#A6D7FE"]
    markers2=['o', 'o', 'o', 'o', 'd', '^', '^', '^', '^', 'o', 'o','o','o']
    colors1 = ["#E31A1C", "#FE8E8F","#984EA3", "#F6BDFE","#377EB8", "#A6D7FE"]
    markers1 =['o', 'd', '^']
    
    #colors = ["#252525", "#252525", "#636363", "#636363", "#969696", "#969696"] #sequential gray
    #colors = ["#1B9E77", "#7FC97F", "#D95F02", "#FDC086", "#7570B3", "#BEAED4"] #color brewer qualitative
    #colors = ["#0868AC", "#0868AC", "#43A2CA", "#43A2CA", "#7BCCC4", "#7BCCC4"] #sequential blue
    colors3 = getColors3()
    markers3=['o', 'o', 'o', 'o', 'o', 'o', '^', '^', 'o', 'o','o','o']
    
    if len(samples) == 3:
        colors = colors1
        markers = markers1
    elif len(samples) == 2:
        colors = ["#E31A1C", "#FE8E8F","#377EB8", "#A6D7FE"]
        markers =['o', '^']
    elif len(samples) == 8 :
        colors = colors2
        markers = markers2
    elif len(samples) == 7:
        colors = [ "#E31A1C", "#FE8E8F", "#E31A1C", "#FE8E8F","#E31A1C", "#FE8E8F","#E31A1C", "#FE8E8F", "#377EB8", "#A6D7FE", "#377EB8", "#A6D7FE","#377EB8", "#A6D7FE"]
        markers = ['o', 'o', 'o', 'o', '^', '^', '^']
    elif len(samples) == 11:
        colors = [ "#E31A1C", "#FE8E8F", "#E31A1C", "#FE8E8F","#E31A1C", "#FE8E8F","#E31A1C", "#FE8E8F", "#E31A1C", "#FE8E8F", "#E31A1C", "#FE8E8F", "#E31A1C", "#FE8E8F", "#377EB8", "#A6D7FE", "#377EB8", "#A6D7FE","#377EB8", "#A6D7FE", "#377EB8", "#A6D7FE"]
        markers = ['o', 'o', 'o', 'o', 'o', 'o', 'o', '^', '^', '^', '^']
    else:
        colors = colors3
        markers = markers3
    #END HACK COLORS

    c = -2
    textsize = 'medium'
    fontP = FontProperties()
    fontP.set_size(textsize)
    linesDict = {}
    labelsDict = {}

    #xtickLabels = getIntersectGenes(samples, genetype)
    xtickLabels = intersectGenes
    #xtickLabels = sortGenes(xtickLabels, samples, genetype)
    
    #if genetype == 'v':
        #xtickLabels = samples[0].vgenes
    #else:
    #    xtickLabels = samples[0].jgenes
    
    #get x location
    xdata = range( 0, len(xtickLabels), 1 )
    maxy = 0
    for i in range( len(axesList) ):
        lines = []
        sampleNames = []
        axes = axesList[i]
        startIndex = i*samplesPerPlot
        endIndex = min( [startIndex + samplesPerPlot, len(samples)] )
        for j in range( startIndex, endIndex ):
            sample = samples[j]
            #sampleNames.append( "%s-%d" % (sample.name, sample.totalCount))
            sampleNames.append( "%s" % (sample.name))
            ydata = getGeneUsage( sample, genetype, xtickLabels, False )
            maxy = max( [maxy, max(ydata)] )

            c += 2
            l = axes.plot( xdata, ydata, color=colors[c], marker=markers[c/2], markersize=8.0, markeredgecolor=colors[c], linestyle='none' )
            axes.plot( xdata, ydata, color=colors[c +1], linestyle='-', linewidth=0.01 )
            lines.append( l )

            #Draw standard deviation if the data is available:
            if sample.name in options.avr2std and options.avr2std[sample.name] in stds:
                std = stds[ options.avr2std[sample.name] ]
                higherYdata, lowerYdata = getStdData( std, sample, genetype, xtickLabels )
                if len(higherYdata) == 0 or len(lowerYdata) == 0:
                    continue
                for i in xrange( len(ydata) ):
                    axes.plot( [xdata[i], xdata[i]], [higherYdata[i], lowerYdata[i]], linestyle='-', marker='_', markersize=8.0, linewidth=0.02, color=colors[c+1] )
                maxy = max([maxy, max(higherYdata)])

        linesDict[i] = lines
        labelsDict[i] = sampleNames
        axes.set_title( '%s Usage Distribution' % genetype.upper(), size='xx-large' )

    for i in range( len(axesList) ):
        axes = axesList[ i ]
        immunoseqLib.editSpine( axes )
        axes.set_xlabel('Genes', size = 'x-large')
        axes.set_ylabel('Frequency', size='x-large')

        #Legend
        if not options.std:
            legend = pyplot.legend( linesDict[i], labelsDict[i], numpoints=1, prop=fontP, loc="best" )
            #legend = axes.legend( linesDict[i], labelsDict[i], numpoints = 1, "upper right", ncol=3)
            #for t in legend.get_texts():
            #    t.set_fontsize('xx-small')
            legend._drawFrame = False

        axes.xaxis.set_ticklabels( xtickLabels )
        axes.xaxis.set_ticks( xdata )
        axes.set_xlim(-0.5, len(xtickLabels) + 0.5)

        numTicks = 20
        yticks = [ float(t)/numTicks for t in range(numTicks +1) ]
        ytickLabels = []
        for y in yticks:
            ytickLabels.append( "%d" %(y*100) )
        axes.yaxis.set_ticklabels( ytickLabels )
        axes.yaxis.set_ticks( yticks )
        axes.set_ylim(-0.01, maxy + 0.01)
        #axes.set_ylim(-0.01, 0.25)
        for label in axes.get_xticklabels():
            label.set_fontsize( textsize )
            label.set_rotation( 90 )
        for label in axes.get_yticklabels():
            label.set_fontsize( textsize )
        axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
        axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    return 

def drawUsageDist( samples, stds, options, genetype, intersectGenes ):
    #Draw V usage distribution of each sample onto the same pdf
    options.out = os.path.join( options.outdir, "%sUsage" %genetype)
    fig, pdf = immunoseqLib.initImage( 10.0, 12.0, options )
    axesList = setUsageDistAxes( fig, len(samples), options.samplesPerPlot )
    drawUsageData( axesList, samples, stds, options, genetype, intersectGenes )
    immunoseqLib.writeImage( fig, pdf, options )

################# TWO-SAMPLE STUDENT'S T TEST TO SEE IF THERE IS BIAS AMONG DIFFERENT GROUPS ######################
#In this section, for each genetype (V, J, or VJ), we will calculate the Pearson Correlation of the relative V usage (frequencies instead of raw counts) of all pairs of samples that are WITHIN the same technology and all pairs of samples from two DIFFERENT technologies. And do a t test to see whether to accept for reject the null hypothesis (there is no bias among the technologies, i.e the mean PearsonCorrelations of pairs with same tech is the same with the mean PearsonCorrelation of pairs with different techs.
def readGroup2samples(file):
    group2samples = {}
    f = open(file, 'r')
    for line in f:
        items = line.strip().split()
        group2samples[items[0]] = items[1].split(',')
    f.close()
    return group2samples

def getCorrVorJ(sample1, sample2, genetype, genes):
    abs = False
    vec1 = getGeneUsage(sample1, genetype, genes, abs)
    vec2 = getGeneUsage(sample2, genetype, genes, abs)
    corr, pval = pearsonr(vec1, vec2)
    return corr

def getVJusage(sample):
    vec = []
    for r in  sample.normIntersectVJusage :
        vec.extend(r)
    return vec

def getCorrVJ(sample1, sample2):
    vec1 = getVJusage(sample1)
    vec2 = getVJusage(sample2)
    corr, pval = pearsonr(vec1, vec2)
    return corr

def getCorr(sample1, sample2, genetype, type2intersectGenes):
    if genetype == 'vj':
        return getCorrVJ(sample1, sample2)
    else:
        return getCorrVorJ(sample1, sample2, genetype, type2intersectGenes[genetype])

def getSampleByName(samples, name):
    for s in samples:
        if s.name == name:
            return s
    return None

def pairTtest(samples, group2samples, genetype, type2intersectGenes):
    #genes = getIntersectGenes(samples, genetype)
    samePairs = []
    diffPairs = []
    
    #Same pairs:
    for g, gsamples in group2samples.iteritems():
        for i in xrange( len(gsamples) - 1 ):
            s1 = getSampleByName(samples, gsamples[i])
            if not s1:
                raise KeyError("Could not found sample %s in the sample list\n" %gsamples[i])
            for j in xrange( i+1, len(gsamples) ):
                s2 = getSampleByName(samples, gsamples[j])
                if not s2:
                    raise KeyError("Could not found sample %s in the sample list\n" %gsamples[j])
                samePairs.append( getCorr(s1, s2, genetype, type2intersectGenes) )
    #Diff pairs:
    groups = group2samples.keys()
    for i in xrange( len(groups) - 1 ):
        for gs1 in group2samples[ groups[i] ]:
            s1 = getSampleByName(samples, gs1)
            for j in xrange( i+1, len(groups) ):
                for gs2 in group2samples[ groups[j] ]:
                    s2 = getSampleByName(samples, gs2)
                    diffPairs.append( getCorr(s1, s2, genetype, type2intersectGenes) )

    #t-test:
    tval, pval = ttest_ind(samePairs, diffPairs)
    return tval, pval, np.mean(samePairs), np.std(samePairs), np.mean(diffPairs), np.std(diffPairs)

def ttests(samples, group2samples, outfile, type2intersectGenes):
    #group2samples = readGroup2samples(group2samplesFile)
    genetypes = ['v', 'j', 'vj']
    f = open(outfile, 'w')
    f.write("Genetype\tPval\tSameMean\tSameStd\tDiffMean\tDiffStd\n")
    for genetype in genetypes:
        tval, pval, sameMean, sameStd, diffMean, diffStd = pairTtest(samples, group2samples, genetype, type2intersectGenes)
        f.write("%s\t%f\t%f\t%f\t%f\t%f\n" %(genetype, pval, sameMean, sameStd, diffMean, diffStd) )
    f.close()

############# DRAW VJ HEATMAPS #####################
def setCompareAxes( fig ):
    """
    Set axes for the VJheatmap Plot. There are 3 subplots total:
    Plot 0: Heatmap plot, where rows = Vs, cols = Js
    Plot 1: Right at the top of plot 0, shows J usage
    Plot 2: To the right of plot 1, shows V usage
    """
    axesList = []
    axleft = 0.1
    axright = 0.98
    axwidth = axright - axleft
    axbottom = 0.08
    axtop = 0.95
    axheight = axtop - axbottom
    margin = 0.07 #space between plots

    plot0wh = 0.7 #plot1 width = plot1 height
    plot1h = axheight - (plot0wh + margin) 
    plot2w = plot1h

    axesList.append( fig.add_axes( [axleft, axbottom, plot0wh, plot0wh] ) ) #Plot 0
    axesList.append( fig.add_axes( [axleft, axbottom + plot0wh + margin, plot0wh, plot1h] ) ) #Plot 1
    axesList.append( fig.add_axes( [axleft + plot0wh + margin, axbottom, plot2w, plot0wh] ) ) #Plot 2
    return axesList

#def drawVJdata( axesList, sample, options ):
def drawVJdata( fig, axes, sample, vgenes, jgenes, options ):
    #Draw heatmap:
    #data = sample.normalizeVJusage()
    data = sample.normIntersectVJusage
    hmaxes = axes.imshow( data, interpolation='nearest' )
    #Colorbar:
    #from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    #axins = inset_axes(axes, width="5%", height="10%", loc=3, bbox_to_anchor=(1.05, 0, 1, 1), bbox_transform=axes.transAxes, borderpad=0)
    #matplotlib.pyplot.colorbar(hmaxes, cax=axins, ticks=[0, 0.5, 1])
    cbar = fig.colorbar(hmaxes, shrink=0.3)
    #cbar.ax.set_yticklabels(['0', '0.5', '1'])

    immunoseqLib.editSpine( axes )
    axes.set_title( "VJ usage of sample %s" %sample.name )
    xticks = [ x for x in range( len(jgenes) ) ]
    axes.xaxis.set_ticks( xticks )
    axes.xaxis.set_ticklabels( jgenes )
    
    yticks = [ y for y in range( len(vgenes) ) ]
    axes.yaxis.set_ticks( yticks )
    axes.yaxis.set_ticklabels( vgenes )
    
    textsize = 'x-small'
    for label in axes.get_xticklabels():
        label.set_fontsize( textsize )
        label.set_rotation( 70 )
    for label in axes.get_yticklabels():
        label.set_fontsize( textsize )

    #Draw Jusage histogram:
    #p1axes = axesList[1]
    #p1axes.bar( range( len(sample.jcounts) ), sample.jcounts )

    #Draw Vusage histogram:
    #p2axes = axesList[2]
    #bottoms = range( len(sample.vgenes) )
    #p2axes.barh( bottoms, sample.vcounts )

def drawVJ( sample, vgenes, jgenes, options ):
    options.out = os.path.join( options.outdir, "%s-vjUsage" %sample.name )
    fig, pdf = immunoseqLib.initImage( 8.0, 10.0, options )
    axes = immunoseqLib.setAxes(fig)
    #axesList = setCompareAxes( fig )

    drawVJdata( fig, axes, sample, vgenes, jgenes, options )
    immunoseqLib.writeImage( fig, pdf, options )
#################### DONE DRAWING VJ HEATMAPS ##################

################### Principal Component Analysis PCA ###########
def getSample2group(group2samples):
    s2g = {}
    for g, samples in group2samples.iteritems():
        for s in samples:
            s2g[s] = g
    return s2g

def preparePcaMatrix(samples, group2samples, type2intersectGenes, genetype):
    m = []
    rownames = [] #sample name of each row
    abs = False 
    for g, names in group2samples.iteritems():
        for n in names:
            rownames.append( (n, g) )
            sample = getSampleByName(samples, n)
            if genetype == 'vj':
                rowdata = getVJusage(sample)
            else:
                rowdata = getGeneUsage( sample, genetype, type2intersectGenes[genetype], abs )
            m.append(rowdata)
    return np.array(m), rownames

def getPcaTransformedMatrix(samples, group2samples, type2intersectGenes, genetype):
    m, rownames = preparePcaMatrix(samples, group2samples, type2intersectGenes, genetype)
    transformedM = mdp.pca(m, output_dim=4)
    for i,r in enumerate(transformedM):
        sys.stdout.write("%s\t%s\t%s\n" %(rownames[i][0], rownames[i][1], "\t".join( [str(c) for c in r] )))

################### END OF Principal Component Analysis PCA ###########

def readfiles( indir, noheader ):
    samples = []
    stds = {} #key = name, val = Sample

    for file in os.listdir( indir ):
        if os.path.isdir(file):
            continue
        #name = os.path.basename( file )
        #name = file.split('-')[0]
        filenameItems = file.split('-')
        name = '-'.join( filenameItems[:len(filenameItems) -1] )
        if name == '':
            name = file.split('.')[0]
        name = name.rstrip('_uniq').rstrip('_abs')
        sample = Sample( name )
        fh = open( os.path.join(indir, file), 'r' )
        m = []
        for line in fh.readlines():
            if re.search('/OR', line):
                continue
            m.append( line.strip().split('\t') )

        if not noheader and len(m) > 0:
            colnames = m.pop(0)
            rownames = []
            for i in range( len(m) ):
                rownames.append( m[i].pop(0) )
        else:
            rownames = [ i for i in range( len(m) ) ]
            colnames = []
            if len(rownames) > 0:
                rownames = [ j for j in range( len(m[0]) ) ]
        
        #Convert value to float:
        for i in range( len(m) ):
            m[i] = [ float(c) for c in m[i] ]

        sample.setUsage(m, rownames, colnames)
        sample.setVcounts()
        sample.setJcounts()
        sample.setTotalCount()
        if re.match('std', sample.name):
            stds[sample.name] = sample
        else:
            samples.append( sample )

    return samples, stds

def checkSamples( samples ):
    if len(samples) == 0:
        return
    vgenes = samples[0].vgenes
    jgenes = samples[0].jgenes
    for s in samples:
        if vgenes != s.vgenes:
            sys.stderr.write("Sample %s has different number of Vs compared with sample %s\n" %(s.name, samples[0].name) )
            sys.exit(1)
        elif jgenes != s.jgenes:
            sys.stderr.write("Sample %s has different number of Js compared with sample %s\n" %(s.name, samples[0].name) )
            sys.exit(1)

def avr2std(file): 
    avr2std = {}
    f = open(file, 'r')
    for line in f:
        items = line.strip().split()
        avr2std[items[0]] = items[1]
    f.close()
    return avr2std

def checkOptions( args, options, parser ):
    if not options.indir or not os.path.isdir( options.indir ):
        parser.error("Input directory is required. Please specify a valid input directory.\n")
    if not os.path.isdir( options.outdir ):
        parser.error("Output directory '%s' does not exist.\n" %options.outdir)
    if options.sampleOrder:
        options.sampleOrder = options.sampleOrder.split(',')
    
    options.avr2std = {}
    if options.std:
        options.avr2std = avr2std(options.std)

def initOptions( parser ):
    parser.add_option('-i', '--indir', dest='indir', help="Required argument: Input directory")
    parser.add_option('-o', '--outdir', dest='outdir', default='.', help="Output directory. Default = current directory")
    parser.add_option('--noheader', dest='noheader', action="store_true", default=False, help="If NOT specified (default), it is assumed that the first row has the matrix column names (J genes), and first column has the row names (vgenes).")
    parser.add_option('--samplesPerPlot', dest='samplesPerPlot', type='int', default=15)
    parser.add_option('--sampleOrder', dest='sampleOrder', type='string', help='Comma seperated list of samples in the order of interest')
    parser.add_option('--std', dest='std', help='Required if wish to add standard deviation for the average usage. File must have the format: <averageUsageFile> <stdFile>. Default=None')
    parser.add_option('--group2samples', dest='group2samples', help='Format: <groupname> <Comma separated list of samples of this group>')
    parser.add_option('--ttest', dest='ttest', action = 'store_true', default=False, help='If specified, will perform ttest')
    parser.add_option('--vj', dest='vj', action = 'store_true', default=False, help='If specified, will draw VJ usage heatmap for each sample')
    #parser.add_option('--vjPlotPerPage', dest='vjPlotPerPage', type='int', default=8)

def main():
    parser = immunoseqLib.initOptions()
    initOptions( parser )
    immunoseqLib.initPlotOptions( parser )
    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    immunoseqLib.checkPlotOptions( options, parser )

    samples, stds = readfiles( options.indir, options.noheader )
    type2intersectGenes = getIntersect(samples)
    group2samples = readGroup2samples(options.group2samples)
    
    #Intersect VJ:
    if options.ttest or options.vj:
        intersectVJ(samples, type2intersectGenes) #added intersectVJusage, normIntersectVJusage

    #Pca:
    #getPcaTransformedMatrix(samples, group2samples, type2intersectGenes, 'v')
    #return

    #Student's t test:
    if options.ttest:
        outfile = os.path.join(options.outdir, "ttest.txt")
        ttests(samples, group2samples, outfile, type2intersectGenes)
    
    #HACK: Sort samples:
    #sampleOrder = 'SBC1,SBC7,SBC3,SBC4,SBC5,SBC6,SBC2,SBC8,Patient-01-R,Patient-01-D,Patient-10,Patient-11,Patient-12,Patient-13,Patient-B-R,Patient-B-D,Patient-8-D,Patient-15-D'
    if options.sampleOrder:
        orders = options.sampleOrder
        sortedSamples = []
        for s in orders:
            for sample in samples:
                if sample.name == s:
                    sortedSamples.append(sample)
        for s in samples:
            if s not in sortedSamples:
                sortedSamples.append(s)
    else:
        sortedSamples = samples

    #checkSamples( samples )
    drawUsageDist( sortedSamples, stds, options, 'v', type2intersectGenes['v'] )
    drawUsageDist( sortedSamples, stds, options, 'j', type2intersectGenes['j'] )
    
    if options.vj:
        for sample in samples:
            drawVJ( sample, type2intersectGenes['v'], type2intersectGenes['j'], options )

if __name__ == "__main__":
    main()
