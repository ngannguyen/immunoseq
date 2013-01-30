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
import matplotlib as mpl
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
                sys.stderr.write("Column labels vector has different length (%d) from the number of columns (%d).\n" %(len(self.jgenes), len(row)))
                sys.stderr.write("Column labels: %s\n" %('\t'.join(self.jgenes)))
                sys.stderr.write("Sample: %s\n" %self.name)
                print row
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
    axtop = 0.92
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

def filterGenes(genes, keepGenes):
    #Get the intersection of genes & keepGenes:
    intersectGenes = []
    for g in keepGenes:
        if g in genes:
            intersectGenes.append(g)
    return intersectGenes

def getSelectedGenes(type2genes, vs, js):
    newtype2genes = {}
    types = {'v': vs, 'j':js}
    for type, selectedGenes in types.iteritems():
        if selectedGenes:
            genes = filterGenes(type2genes[type], selectedGenes)
            newtype2genes[type] = genes
        else:
            newtype2genes[type] = type2genes[type]
    return newtype2genes

def getIntersect(samples):
    type2genes = {} 
    types=['v', 'j']
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

def getStdData(std, sample, genetype, genes, abs):
    stdData = getGeneUsage(std, genetype, genes, abs)
    ydata = getGeneUsage(sample, genetype, genes, abs)
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
               "#AA6600", "#B5A642",
               "#E31A1C", "#FE8E8F",
               "#984EA3", "#F6BDFE",
               ] #
    return colors

def drawUsageData( axesList, samples, stds, options, genetype, intersectGenes ):
    if len( samples ) <= 0:
        return
    samplesPerPlot = options.samplesPerPlot 
    
    #HACK COLORS
    #colors = []
    #if options.sampleColor:
    #    colors = options.sampleColor
    #if options.sampleMarker:
    #    markers = options.sampleMarker

    #if options.sampleOrder and not options.sampleColor:
    #    colors = [ "#E31A1C", "#FE8E8F"]*len(options.sampleOrder) + ["#377EB8", "#A6D7FE"]*(len(samples) - len(options.sampleOrder))
    #    markers = ["o"]*len(options.sampleOrder) + ["^"]*(len(samples) - len(options.sampleOrder))
    #if not options.groupDrawInfo and len(colors) < len(samples):
    #    raise ValueError("Too many samples (%d), not enough colors (%d). Please specify the options sampleColor and sampleMarker to fix the problem." %(len(samples), len(colors)))
    #END HACK COLORS

    c = -2
    textsize = 'large'
    fontP = FontProperties()
    fontP.set_size(textsize)
    linesDict = {}
    labelsDict = {}

    xtickLabels = intersectGenes
    
    #get x location
    xdata = range( 0, len(xtickLabels), 1 )
    maxy = 0

    #HACK COLORS:
    s2group = {'as1D': 'B27+,AS', 'as8D': 'B27+,AS', 'as11D': 'B27+,AS', 'as16D': 'B27+,AS',
               'as15D': 'B27-,AS', 
               'asBD': 'B27+,Healthy', 'as20D': 'B27+,Healthy',
               'adaptMA': 'B27-,Healthy', 'adaptM35': 'B27-,Healthy', 'adaptF57': 'B27-,Healthy', 'adaptAS': 'B27-,Healthy', 'adaptCD': 'B27-,Healthy', 'adaptF28':'B27-,Healthy',
               'irep1D': 'B27+,AS', 'as10R': 'B27+,AS', 'as11R': 'B27+,AS', 'as12R': 'B27+,AS', 'as13R': 'B27+,AS', 'as1R': 'B27+,AS', 'irepBD': 'B27+,Healthy', 'asBR': 'B27+,Healthy',
               'sameGroupShared': 'AS', 'diffGroupShared': 'Healthy', 'uniq': 'Unique',
               'b27pos': 'B27+', 'b27neg':'B27-'}

    #B27 +-
    #g2c = {'B27+,AS': '#E31A1C', 'B27+,Healthy': '#E31A1C', 'B27-,AS':'#377EB8', 'B27-,Healthy':'#377EB8', 'B27+':'#E31A1C', 'B27-': '#377EB8'}
    #g2lc = {'B27+,AS': '#FE8E8F', 'B27+,Healthy': '#FE8E8F', 'B27-,AS':'#A6D7FE', 'B27-,Healthy':'#A6D7FE', 'B27+': '#FE8E8F', 'B27-': '#A6D7FE'}
    #g2m = {'B27+,AS': 'o', 'B27+,Healthy': '^', 'B27-,AS':'o', 'B27-,Healthy':'^', 'B27+':'o', 'B27-':'^'}
    
    #AS vs Healthy
    #g2c = {'B27+,AS': '#E31A1C', 'B27+,Healthy': '#9E1114', 'AS':'#377EB8', 'Healthy':'#275880', 'Unique': '#4DAF4A'}
    #g2lc = {'B27+,AS': '#FE8E8F', 'B27+,Healthy': '#FE8E8F', 'AS':'#A6D7FE', 'Healthy':'#A6D7FE', 'Unique': '#B8FEB5'}
    #g2m = {'B27+,AS': 'o', 'B27+,Healthy': '^', 'AS':'o', 'Healthy':'^', 'Unique': 's'}
    
    #AS vs Healthy, no color (gray scale)
    g2c = {'B27+,AS': '#353535', 'B27+,Healthy': '#848484', 'AS':'#353535', 'Healthy':'#848484', 'Unique': '#4DAF4A'}
    g2lc = {'B27+,AS': '#BDBDBD', 'B27+,Healthy': '#BDBDBD', 'AS':'#BDBDBD', 'Healthy':'#BDBDBD', 'Unique': '#B8FEB5'}
    g2m = {'B27+,AS': 'o', 'B27+,Healthy': '^', 'AS':'o', 'Healthy':'^', 'Unique': 's'}
    #END HACK COLORS

    if options.groupDrawInfo:
        s2group = options.s2g
        g2c = options.g2c
        g2lc = options.g2lc
        g2m = options.g2m
    
    for i in range( len(axesList) ):
        lines = []
        sampleNames = []
        axes = axesList[i]
        startIndex = i*samplesPerPlot
        endIndex = min( [startIndex + samplesPerPlot, len(samples)] )
        groups = []
        for j in range( startIndex, endIndex ):
            sample = samples[j]
            sampleNames.append( "%s" % (sample.name))
            ydata = getGeneUsage( sample, genetype, xtickLabels, options.abs )
            maxy = max( [maxy, max(ydata)] )

            #c += 2
            #l = axes.plot( xdata, ydata, color=colors[c], marker=markers[c/2], markersize=8.0, markeredgecolor=colors[c], linestyle='none' )
            #axes.plot( xdata, ydata, color=colors[c +1], linestyle='-', linewidth=0.01 )
            #lines.append( l )
            
            group = s2group[sample.name]
            l = axes.plot( xdata, ydata, color=g2c[group], marker=g2m[group], markersize=8.0, markeredgecolor=g2c[group], linestyle='none' )
            axes.plot( xdata, ydata, color=g2lc[group], linestyle='-', linewidth=0.01 )
            if group not in groups:
                lines.append( l )
                groups.append(group)
            
            #Draw standard deviation if the data is available:
            if sample.name in options.avr2std and options.avr2std[sample.name] in stds:
                std = stds[ options.avr2std[sample.name] ]
                higherYdata, lowerYdata = getStdData( std, sample, genetype, xtickLabels, options.abs )
                if len(higherYdata) == 0 or len(lowerYdata) == 0:
                    continue
                for i in xrange( len(ydata) ):
                    #axes.plot( [xdata[i], xdata[i]], [higherYdata[i], lowerYdata[i]], linestyle='-', marker='_', markersize=8.0, linewidth=0.02, color=colors[c+1] )
                    axes.plot( [xdata[i], xdata[i]], [higherYdata[i], lowerYdata[i]], linestyle='-', marker='_', markersize=8.0, linewidth=0.02, color=g2lc[group] )
                maxy = max([maxy, max(higherYdata)])

        linesDict[i] = lines
        #labelsDict[i] = sampleNames
        labelsDict[i] = groups
        #axes.set_title( 'TRB%s Usage Distribution' % genetype.upper(), size='xx-large', weight='bold' )

    #HACK
    for i, xlabel in enumerate(xtickLabels):
        if xlabel == 'TRBV6-5/TRBV6-6':
            xtickLabels[i] = 'TRBV6-5/6-6'
    xtickLabels = [xlabel.lstrip('TRB') for xlabel in xtickLabels]
    for i, xlabel in enumerate(xtickLabels):
        items = xlabel.split('|')
        items = [item.lstrip('TRB') for item in items]
        xtickLabels[i] = '|'.join(items)
    #END HACK

    for i in range( len(axesList) ):
        axes = axesList[ i ]
        immunoseqLib.editSpine( axes )
        axes.set_xlabel('Gene', size = 'xx-large', weight='bold')
        axes.set_ylabel('Frequency (% of total sequences)', size='xx-large', weight='bold')

        #Legend
        if not options.std:
            #HACK
            cat2index = {}
            for j, label in enumerate(labelsDict[i]):
                cat2index[label] = j
            labels = cat2index.keys()
            lines = [linesDict[i][cat2index[label]] for label in labels]
            legend = pyplot.legend( lines, labels, numpoints=1, prop=fontP, loc="best" )
            #END HACK

            #legend = pyplot.legend( linesDict[i], labelsDict[i], numpoints=1, prop=fontP, loc="best" )
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
        #axes.set_ylim(-0.01, 0.22)
        for label in axes.get_xticklabels():
            label.set_fontsize( 'small' )
            label.set_fontweight( 'bold' )
            label.set_rotation( 75 )
        for label in axes.get_yticklabels():
            label.set_fontsize( textsize )
            label.set_fontweight( 'bold' )
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
    if not file:
        return group2samples
    f = open(file, 'r')
    for line in f:
        items = line.strip().split()
        group2samples[items[0]] = items[1].split(',')
    f.close()
    return group2samples

#========== ttest for each gene (targetGroup has higher/smaller usage than otherGroup) ============
def getGene2countVorJ(sample, genetype, genes, abs):
    gene2count = {}
    if genetype == 'v':
        for i, gene in enumerate(sample.vgenes):
            if gene in genes:
                gene2count[gene] = sample.vcounts[i]
                if not abs:
                    gene2count[gene] /= sample.totalCount
    else: 
        for i, gene in enumerate(sample.jgenes):
            if gene in genes:
                gene2count[gene] = sample.jcounts[i]
                if not abs:
                    gene2count[gene] /= sample.totalCount
    if len(gene2count.keys()) != len(genes):
        raise ValueError("Some gene in the gene list %s are not present in sample genes %s\n" %(','.join(genes), ','.join(sample.vgenes)))
    return gene2count

def getGene2countVJ(sample, selectedvjs, vs, js, abs):
    gene2count = {}
    if abs:
        vjusage = sample.intersectVJusage
    else:
        vjusage = sample.normIntersectVJusage
    for i, r in  enumerate( vjusage ) :
        for j, c in enumerate(r):
            vname = vs[i]
            jname = js[j]
            vjname = "%s|%s" %(vname, jname)
            if not selectedvjs or vjname in selectedvjs:
                gene2count[vjname] = c
    return gene2count

def pairTtestSingleGene(samples, group2samples, genetype, type2intersectGenes, selectedvjs, targetGroup, abs):
    '''ttest for each single gene
    '''
    #genes = getIntersectGenes(samples, genetype)
    gene2results = {} #key = gene, val = [tval, pval, meanGroup1, stdGroup1, meanGroup2, stdGroup2]
    gene2counts1 = {} #key = gene, val = list of counts of targetGroup's samples
    gene2counts2 = {} #key = gene, val = list of counts of other samples (non-targetGroup)
    
    targetSamples = group2samples[targetGroup]
    
    for sample in samples:
        if genetype == 'vj':
            gene2count = getGene2countVJ(sample, selectedvjs, type2intersectGenes['v'], type2intersectGenes['j'], abs)
        else:
            gene2count = getGene2countVorJ(sample, genetype, type2intersectGenes[genetype], abs)
        if sample.name in targetSamples:
            for gene, count in gene2count.iteritems():
                if gene not in gene2counts1:
                    gene2counts1[gene] = [count]
                else:
                    gene2counts1[gene].append(count)
        else:
            for gene, count in gene2count.iteritems():
                if gene not in gene2counts2:
                    gene2counts2[gene] = [count]
                else:
                    gene2counts2[gene].append(count)

    #t-test:
    for gene, vec1 in gene2counts1.iteritems():
        vec2 = gene2counts2[gene]
        tval, pval = ttest_ind(vec1, vec2)
        gene2results[gene] = [tval, pval, np.mean(vec1), np.std(vec1), np.mean(vec2), np.std(vec2)]
    return gene2results

#=========== ttests of correlation of gene usage among samples (hypothesis: same group has higher correlation) =================
def getCorrVorJ(sample1, sample2, genetype, genes, abs):
    #abs = False
    vec1 = getGeneUsage(sample1, genetype, genes, abs)
    vec2 = getGeneUsage(sample2, genetype, genes, abs)
    corr, pval = pearsonr(vec1, vec2)
    return corr

def getVJusage(sample, selectedvjs, vs, js, abs):
    #vs = rownames, js = column names of sample.normIntersectVJusage
    vec = []
    if abs:
        vjusage = sample.intersectVJusage
    else:
        vjusage = sample.normIntersectVJusage
    #for i, r in  enumerate( sample.normIntersectVJusage ) :
    for i, r in  enumerate( vjusage ) :
        for j, c in enumerate(r):
            vname = vs[i]
            jname = js[j]
            vjname = "%s|%s" %(vname, jname)
            if not selectedvjs or vjname in selectedvjs:
                vec.append(c)
    return vec

def getCorrVJ(sample1, sample2, type2intersectGenes, selectedvjs, abs):
    vs = type2intersectGenes['v']
    js = type2intersectGenes['j']
    vec1 = getVJusage(sample1, selectedvjs, vs, js, abs)
    vec2 = getVJusage(sample2, selectedvjs, vs, js, abs)
    corr, pval = pearsonr(vec1, vec2)
    return corr

def getCorr(sample1, sample2, genetype, type2intersectGenes, selectedvjs, abs):
    if genetype == 'vj':
        return getCorrVJ(sample1, sample2, type2intersectGenes, selectedvjs, abs)
    else:
        return getCorrVorJ(sample1, sample2, genetype, type2intersectGenes[genetype], abs)

def getSampleByName(samples, name):
    for s in samples:
        if s.name == name:
            return s
    return None

def pairTtest(samples, group2samples, genetype, type2intersectGenes, selectedvjs, targetGroup, abs):
    #genes = getIntersectGenes(samples, genetype)
    samePairs = [] #all pairs of samples in the targetGroup
    diffPairs = [] #all pairs of samples in the other groups (all intra-group pairs, no inter-group)
    
    #Same pairs:
    for g, gsamples in group2samples.iteritems():
        for i in xrange( len(gsamples) - 1 ):
            s1 = getSampleByName(samples, gsamples[i])
            if not s1:
                print [s.name for s in samples]
                raise KeyError("Could not found sample %s in the sample list\n" %gsamples[i])
            for j in xrange( i+1, len(gsamples) ):
                s2 = getSampleByName(samples, gsamples[j])
                if not s2:
                    print [s.name for s in samples]
                    raise KeyError("Could not found sample %s in the sample list\n" %gsamples[j])
                if targetGroup and targetGroup == g:
                    samePairs.append( getCorr(s1, s2, genetype, type2intersectGenes, selectedvjs, abs) )
                else:
                    diffPairs.append( getCorr(s1, s2, genetype, type2intersectGenes, selectedvjs, abs) )

    #t-test:
    tval, pval = ttest_ind(samePairs, diffPairs)
    return tval, pval, np.mean(samePairs), np.std(samePairs), np.mean(diffPairs), np.std(diffPairs)

def pairTtest_intraVsInterGroups(samples, group2samples, genetype, type2intersectGenes, selectedvjs, targetGroup, abs):
    #genes = getIntersectGenes(samples, genetype)
    samePairs = []
    diffPairs = []
    
    #Same pairs:
    for g, gsamples in group2samples.iteritems():
        issame = True
        if targetGroup and g != targetGroup: #if targetGroup is specified, samePairs includes pairs of the target group only
            issame = False
        for i in xrange( len(gsamples) - 1 ):
            s1 = getSampleByName(samples, gsamples[i])
            if not s1:
                print [s.name for s in samples]
                raise KeyError("Could not found sample %s in the sample list\n" %gsamples[i])
            for j in xrange( i+1, len(gsamples) ):
                s2 = getSampleByName(samples, gsamples[j])
                if not s2:
                    print [s.name for s in samples]
                    raise KeyError("Could not found sample %s in the sample list\n" %gsamples[j])
                if issame:
                    samePairs.append( getCorr(s1, s2, genetype, type2intersectGenes, selectedvjs, abs) )
                #else:
                #    diffPairs.append( getCorr(s1, s2, genetype, type2intersectGenes, selectedvjs, abs) )
    #Diff pairs:
    groups = group2samples.keys()
    for i in xrange( len(groups) - 1 ):
        for gs1 in group2samples[ groups[i] ]:
            s1 = getSampleByName(samples, gs1)
            for j in xrange( i+1, len(groups) ):
                for gs2 in group2samples[ groups[j] ]:
                    s2 = getSampleByName(samples, gs2)
                    diffPairs.append( getCorr(s1, s2, genetype, type2intersectGenes, selectedvjs, abs) )

    #t-test:
    tval, pval = ttest_ind(samePairs, diffPairs)
    return tval, pval, np.mean(samePairs), np.std(samePairs), np.mean(diffPairs), np.std(diffPairs)

def ttests(samples, group2samples, outfile, type2intersectGenes, selectedvjs, targetGroup, abs, pvalcutoff):
    #group2samples = readGroup2samples(group2samplesFile)
    genetypes = ['v', 'j', 'vj']
    f = open(outfile, 'w')
    f.write("Genetype\tPval\tSameMean\tSameStd\tDiffMean\tDiffStd\n")
    for genetype in genetypes:
        tval, pval, sameMean, sameStd, diffMean, diffStd = pairTtest(samples, group2samples, genetype, type2intersectGenes, selectedvjs, targetGroup, abs)
        f.write("%s\t%f\t%f\t%f\t%f\t%f\n" %(genetype, pval, sameMean, sameStd, diffMean, diffStd) )

    #Separate genes:
    if targetGroup:
        for genetype in genetypes:
            f.write("#\n#%s\n" %genetype)
            gene2results = pairTtestSingleGene(samples, group2samples, genetype, type2intersectGenes, selectedvjs, targetGroup, abs)
            
            pfile = os.path.join( os.path.split(outfile)[0], "%s-pvals.txt" %genetype) #also print all pvalues of (separate) ttests of all genes:
            pfh = open( pfile, 'w' )
            for gene, results in gene2results.iteritems():
                pval = results[1]
                if pval <= pvalcutoff:
                    f.write("%s\t%s\n" %(gene, '\t'.join([str(c) for c in results[1:]])) )
                #Hack:
                #pfh.write("%s\t%f\n" %(gene, pval))
                maxMean = max([ results[2], results[4] ])
                #if (maxMean <= 1 and maxMean >= 0.01) or maxMean >= 1000 : 
                if (maxMean <= 1 and maxMean >= 0.0001) or maxMean >= 1000 : 
                    pfh.write("%s\t%s\n" %(gene, '\t'.join([str(c) for c in results[1:]])) )
            pfh.close()
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

def getMinMaxVJusage(samples, abs):
    minvj = float('inf')
    maxvj = 0
    for sample in samples:
        if abs:
            data = sample.intersectVJusage
        else:
            data = sample.normIntersectVJusage
        for r in data:
            minvj = min([min(r), minvj])
            maxvj = max([max(r), maxvj])
    return minvj, maxvj

#def drawVJdata( axesList, sample, options ):
def drawVJdata( fig, axes, sample, vgenes, jgenes, options, minvj, maxvj ):
    #Draw heatmap:
    #data = sample.normalizeVJusage()
    if options.abs:
        data = sample.intersectVJusage
    else:
        data = sample.normIntersectVJusage

    #Normalize data to the range minvj-maxvj:
    if options.heatmapNoScale:
        hmaxes = axes.imshow( data, interpolation='nearest' )
    else:
        norm = mpl.colors.Normalize(vmin=minvj, vmax=maxvj)
        cmap = mpl.cm.get_cmap('rainbow', 20)
        hmaxes = axes.imshow( data, interpolation='nearest', cmap=cmap, norm=norm )
        
    #Colorbar:
    #from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    #axins = inset_axes(axes, width="5%", height="10%", loc=3, bbox_to_anchor=(1.05, 0, 1, 1), bbox_transform=axes.transAxes, borderpad=0)
    #matplotlib.pyplot.colorbar(hmaxes, cax=axins, ticks=[0, 0.5, 1])
    cbar = fig.colorbar(hmaxes, shrink=0.3)
    #cbar.ax.set_yticklabels(['0', '0.5', '1'])

    immunoseqLib.editSpine( axes )
    axes.set_title( "VJ usage of sample %s" % immunoseqLib.properName(sample.name) )
    xticks = [ x for x in range( len(jgenes) ) ]
    axes.xaxis.set_ticks( xticks )
    for i, xlabel in enumerate(jgenes):
        items = xlabel.split('|')
        items = [item.lstrip('TRB') for item in items]
        jgenes[i] = '|'.join(items)
    axes.xaxis.set_ticklabels( jgenes )
    
    yticks = [ y for y in range( len(vgenes) ) ]
    axes.yaxis.set_ticks( yticks )
    for i, ylabel in enumerate(vgenes):
        items = ylabel.split('|')
        items = [item.lstrip('TRB') for item in items]
        vgenes[i] = '|'.join(items)
    axes.yaxis.set_ticklabels( vgenes )
    
    textsize = 'x-small'
    for label in axes.get_xticklabels():
        label.set_fontsize( textsize )
        label.set_rotation( 80 )
    for label in axes.get_yticklabels():
        label.set_fontsize( textsize )

    #Draw Jusage histogram:
    #p1axes = axesList[1]
    #p1axes.bar( range( len(sample.jcounts) ), sample.jcounts )

    #Draw Vusage histogram:
    #p2axes = axesList[2]
    #bottoms = range( len(sample.vgenes) )
    #p2axes.barh( bottoms, sample.vcounts )

def drawVJ( sample, vgenes, jgenes, options, minvj, maxvj ):
    options.out = os.path.join( options.outdir, "%s-vjUsage" %sample.name )
    fig, pdf = immunoseqLib.initImage( 8.0, 10.0, options )
    axes = immunoseqLib.setAxes(fig)
    #axesList = setCompareAxes( fig )

    drawVJdata( fig, axes, sample, vgenes, jgenes, options, minvj, maxvj )
    immunoseqLib.writeImage( fig, pdf, options )
#################### DONE DRAWING VJ HEATMAPS ##################

################### Principal Component Analysis PCA ###########
def getSample2group(group2samples):
    s2g = {}
    for g, samples in group2samples.iteritems():
        for s in samples:
            s2g[s] = g
    return s2g

def drawPcaData(axes, rownames, transformedM, options):
    group2xdata = {}
    group2ydata = {}

    for i,r in enumerate(transformedM):
        sample = rownames[i][0]
        group = rownames[i][1]
        x = r[0]
        y = r[1]
        if group not in group2xdata:
            group2xdata[group] = [x]
            group2ydata[group] = [y]
        else:
            group2xdata[group].append(x)
            group2ydata[group].append(y)

    groups = sorted(group2xdata.keys())
    if options.groupDrawInfo:
        g2c = options.g2c
        g2m = options.g2m
        colors = [g2c[g] for g in groups]
        markers = [g2m[g] for g in groups]
    else:
        colors = immunoseqLib.getColors6()
        markers = ['^', 'o', 'd', 'p', 'v', '*', 's']
    
    lines = []
    minx = float('inf')
    maxx = float('-inf')
    miny = float('inf')
    maxy = float('-inf')
    for i, group in enumerate(groups):
        if i >= len(colors):
            raise ValueError("drawPcaData: Need more color!")
        color = colors[i]
        marker = markers[i]
        xdata = group2xdata[group]
        ydata = group2ydata[group]
        
        minx = min(minx, min(xdata))
        miny = min(miny, min(ydata))
        maxx = max(maxx, max(xdata))
        maxy = max(maxy, max(ydata))

        l = axes.plot(xdata, ydata, color=color, marker=marker, markersize=15.0, markeredgecolor=color, linestyle='none' )
        lines.append(l)
    
    for label in axes.get_xticklabels():
        label.set_fontsize( 'large' )
        label.set_fontweight( 'bold' )
    for label in axes.get_yticklabels():
        label.set_fontsize( 'large' )
        label.set_fontweight( 'bold' )

    rangex = maxx - minx
    axes.set_xlim(minx - rangex*0.1, maxx + rangex*0.1)
    rangey = maxy - miny
    axes.set_ylim(miny - rangey*0.1, maxy + rangey*0.1 )
    
     
    linenames = []
    for g in groups:
        linenames.append( immunoseqLib.properName(g) )
    #legend = pyplot.legend( lines, sorted(group2xdata.keys()), numpoints=1, loc="lower left" )
    legend = pyplot.legend( lines, linenames, numpoints=1, loc="best" )
    #legend._drawFrame = False
    #axes.set_title( 'PCA', size='xx-large', weight='bold' )
    axes.set_xlabel('PC1', size='xx-large', weight='bold')
    axes.set_ylabel('PC2', size='xx-large', weight='bold')
    axes.yaxis.grid(b=True, color="#3F3F3F", linestyle='-', linewidth=0.5)
    axes.xaxis.grid(b=True, color="#3F3F3F", linestyle='-', linewidth=0.5)

def drawPca(rownames, transformedM, outfile, options):
    #Draw V usage distribution of each sample onto the same pdf
    options.out = outfile
    fig, pdf = immunoseqLib.initImage( 12.5, 8.0, options )
    axes = immunoseqLib.setAxes(fig)
    drawPcaData( axes, rownames, transformedM, options )
    immunoseqLib.writeImage( fig, pdf, options )

def preparePcaMatrix(samples, group2samples, type2intersectGenes, selectedvjs, genetype, abs):
    m = []
    rownames = [] #sample name of each row
    #abs = False 
    for g, names in group2samples.iteritems():
        for n in names:
            rownames.append( (n, g) )
            sample = getSampleByName(samples, n)
            if genetype == 'vj':
                rowdata = getVJusage(sample, selectedvjs, type2intersectGenes['v'], type2intersectGenes['j'], abs)
            else:
                rowdata = getGeneUsage( sample, genetype, type2intersectGenes[genetype], abs )
            m.append(rowdata)
    return np.array(m), rownames

def getPcaTransformedMatrix(samples, group2samples, type2intersectGenes, selectedvjs, genetype, abs, outfile, options):
    m, rownames = preparePcaMatrix(samples, group2samples, type2intersectGenes, selectedvjs, genetype, abs)
    transformedM = mdp.pca(m, output_dim=4)
    
    #Write to text file
    f = open("%s.txt" %outfile, 'w')
    for i,r in enumerate(transformedM):
        f.write("%s\t%s\t%s\n" %(rownames[i][0], rownames[i][1], "\t".join( [str(c) for c in r] )))
    f.close()

    #Draw plot:
    drawPca(rownames, transformedM, outfile, options)

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
        #name = name.rstrip('_uniq').rstrip('_abs')
        name = name.split('_')[0]
        sample = Sample( name )
        fh = open( os.path.join(indir, file), 'r' )
        m = []
        for line in fh.readlines():
            if re.search('/OR', line) or re.search('undefined', line):
                continue
            if len(line.strip()) == 0:
                continue
            m.append( line.strip().split('\t') )

        if not noheader and len(m) > 0:
            colnames = m.pop(0)
            rownames = []
            #for i,row in enumerate(m):
            i = 0
            while i < len(m):
                row = m[i]
                if len(row) != len(colnames) + 1:
                    m.remove(m[i])
                    continue

                rownames.append( m[i].pop(0) )
                i += 1
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

def getGroupDrawInfo(file):
    # Format: <groupname> <comma,sep,sample,list> <groupcolor> <groupLightColor> <groupmarker>
    s2g = {}
    g2c = {}
    g2lc = {}
    g2m = {}
    f = open(file, 'r')
    for line in f:
        if len(line) == 0 or line[0] == '#':
            continue
        items = line.strip().split()
        if len(items) < 5:
            raise ValueError("Wrong groupDrawInfo format. Expected 5 fields, only has %d fields. Line %s\n" %(len(items), line))
    
        group = items[0]
        samples = items[1].split(',')
        color = items[2]
        lightcolor = items[3]
        marker = items[4]
        for s in samples:
            s2g[s] = group
        g2c[group] = color
        g2lc[group] = lightcolor
        g2m[group] = marker
    f.close()
    return s2g, g2c, g2lc, g2m

def checkOptions( args, options, parser ):
    if not options.indir or not os.path.isdir( options.indir ):
        parser.error("Input directory is required. Please specify a valid input directory.\n")
    if not os.path.isdir( options.outdir ):
        parser.error("Output directory '%s' does not exist.\n" %options.outdir)
    if options.sampleOrder:
        options.sampleOrder = options.sampleOrder.split(',')
    if options.sampleColor:
        options.sampleColor = options.sampleColor.split(',')
    if options.sampleMarker:
        options.sampleMarker = options.sampleMarker.split(',')
    if options.groupDrawInfo:
        options.s2g , options.g2c, options.g2lc, options.g2m = getGroupDrawInfo(options.groupDrawInfo)

    options.avr2std = {}
    if options.std:
        options.avr2std = avr2std(options.std)
    sep = ','
    if options.vs:
        options.vs = immunoseqLib.readList(options.vs, sep)
    if options.js:
        options.js = immunoseqLib.readList(options.js, sep)
    if options.vjs:
        options.vjs = immunoseqLib.readList(options.vjs, sep)


def initOptions( parser ):
    parser.add_option('-i', '--indir', dest='indir', help="Required argument: Input directory")
    parser.add_option('-o', '--outdir', dest='outdir', default='.', help="Output directory. Default = current directory")
    parser.add_option('--noheader', dest='noheader', action="store_true", default=False, help="If NOT specified (default), it is assumed that the first row has the matrix column names (J genes), and first column has the row names (vgenes).")
    parser.add_option('-a', '--abs', dest='abs', default=False, action='store_true', help='If specified, will plot the absolute counts instead of pecentage')
    parser.add_option('--samplesPerPlot', dest='samplesPerPlot', type='int', default=15)
    parser.add_option('--std', dest='std', help='Required if wish to add standard deviation for the average usage. File must have the format: <averageUsageFile> <stdFile>. Default=None')
    parser.add_option('--group2samples', dest='group2samples', help='Format: <groupname> <Comma separated list of samples of this group>')
    parser.add_option('--ttest', dest='ttest', action = 'store_true', default=False, help='If specified, will perform ttest')
    parser.add_option('--ttestTargetGroup', dest='ttestTargetGroup', help='If specified, will compare samples belong to this group versus all other samples')
    parser.add_option('--vj', dest='vj', action = 'store_true', default=False, help='If specified, will draw VJ usage heatmap for each sample')
    parser.add_option('--vs', dest='vs', help='File containing a list of V genes of interest. If specified, only look at usage of these genes')
    parser.add_option('--js', dest='js', help='File containing a list of J genes of interest. If specified, only look at usage of these genes')
    parser.add_option('--vjs', dest='vjs', help='File containing a list of VJ gene combinations of interest. If specified, only look at usage of these genes')
    parser.add_option('-p', '--pval', dest='pval', default=1.0, type='float', help='pvalue cutoff. Only report ttests with pval <= this value')
    parser.add_option('--heatmapNoScale', dest='heatmapNoScale', action='store_true', default=False, help='Normally, the VJ heatmap plots are scale to the min and the max values of all samples. If this option is specified, heatmap of each sample is scale to its own range')
    parser.add_option('--noplot', dest='noplot', action='store_true', default=False, help='If specified, do not draw geneusage plot')
    parser.add_option('--pca', dest='pca', action='store_true', default=False, help='If specified, print out the transform pca matrices')
    parser.add_option('--sampleOrder', dest='sampleOrder', type='string', help='Optional. Comma seperated list of samples in the order of interest')
    parser.add_option('--sampleColor', dest='sampleColor', type='string', help='Optional. Comma seperated list of sample color in the order of sampleOrder. Default=None')
    parser.add_option('--sampleMarker', dest='sampleMarker', type='string', help='Optional. Comma seperated list of sample Marker in the order of sampleOrder. Default=None')
    parser.add_option('--groupDrawInfo', dest='groupDrawInfo', help='Optional. File with plot drawing instructions. Format: <groupname> <comma,sep,sample,list> <groupcolor> <groupLightColor> <groupmarker>')
    #parser.add_option('--vjPlotPerPage', dest='vjPlotPerPage', type='int', default=8)

def main():
    parser = immunoseqLib.initOptions()
    initOptions( parser )
    immunoseqLib.initPlotOptions( parser )
    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    immunoseqLib.checkPlotOptions( options, parser )

    samples, stds = readfiles( options.indir, options.noheader )
    type2intersectGenes = getIntersect(samples) #Get genes that are present in all samples
    type2selectedGenes = getSelectedGenes(type2intersectGenes, options.vs, options.js) #Select only genes of interest that are present in all samples
    
    group2samples = readGroup2samples(options.group2samples)
    
    #Intersect VJ:
    if options.ttest or options.vj or options.pca:
        #intersectVJ(samples, type2intersectGenes, options.vjs) #added intersectVJusage, normIntersectVJusage
        intersectVJ(samples, type2selectedGenes) #added intersectVJusage, normIntersectVJusage

    #Pca:
    if options.pca:
        vPcaOut = os.path.join(options.outdir, "vPca")
        getPcaTransformedMatrix(samples, group2samples, type2intersectGenes, options.vjs, 'v', options.abs, vPcaOut, options)
        vjPcaOut = os.path.join(options.outdir, "vjPca")
        getPcaTransformedMatrix(samples, group2samples, type2intersectGenes, options.vjs, 'vj', options.abs, vjPcaOut, options)
    #return

    #Student's t test:
    if options.ttest:
        outfile = os.path.join(options.outdir, "ttest.txt")
        ttests(samples, group2samples, outfile, type2selectedGenes, options.vjs, options.ttestTargetGroup, options.abs, options.pval)
    
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
    if not options.noplot:
        drawUsageDist( sortedSamples, stds, options, 'v', type2selectedGenes['v'] )
        drawUsageDist( sortedSamples, stds, options, 'j', type2selectedGenes['j'] )
    
    if options.vj:
        minvj, maxvj = getMinMaxVJusage(samples, options.abs)
        for sample in samples:
            drawVJ( sample, type2selectedGenes['v'], type2selectedGenes['j'], options, minvj, maxvj )

if __name__ == "__main__":
    main()
