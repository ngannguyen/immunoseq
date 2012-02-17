#!/usr/bin/env python2.6

"""
08/30/2011 nknguyen at soe ucsc edu

Input: location of vjUsage matrices where rows = Vs and columns= Js.
Each cell X_ij contains the (clonotype/read) count that have the combination V_i and J_j.
Each file corresponds to one sample.

Output: 3D histogram of Vusage, Jusage (x = genes, y = counts, z = samples)
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
    
    def getVfreqs2( self, genes ):
        if self.totalCount == 0:
            return []
        freqs = []
        for g in genes:
            if g not in self.vgenes:
                raise KeyError("Gene %s should be in vgenes, but were not found of %s. Vgenes: %s\n" %(g, self.name, ",".join(self.vgenes)))
            for i, currgene in enumerate(self.vgenes):
                if g == currgene:
                    freqs.append( self.vcounts[i]/self.totalCount )
        return freqs

    def getJfreqs( self ):
        if self.totalCount == 0:
            return []
        return [ c/self.totalCount for c in self.jcounts ]

    def getJfreqs2( self, genes ):
        if self.totalCount == 0:
            return []
        freqs = []
        for g in genes:
            if g not in self.jgenes:
                raise KeyError("Gene %s should be in jgenes, but were not found of %s. Vgenes: %s\n" %(g, self.name, ",".join(self.jgenes)))
            for i, currgene in enumerate(self.jgenes):
                if g == currgene:
                    freqs.append( self.jcounts[i]/self.totalCount )
        return freqs
     
    def normalizeVJusage( self ):
        m = []
        for r in self.vjusage:
            row = [ float(c)/self.totalCount for c in r ]
            m.append( row )
        return m

def setUsageDistAxes( fig, numSamples, samplesPerPlot ):
    numaxes = numSamples / samplesPerPlot 
    if numSamples % samplesPerPlot != 0:
        numaxes += 1
    
    axesList = []
    axleft = 0.08
    axright = 0.95
    axwidth = axright - axleft
    axbottom = 0.1
    axtop = 0.9
    axheight = axtop - axbottom
    margin = 0.1

    h = ( axheight - ( margin * (numaxes - 1) ) )/float( numaxes )

    bottom = axtop - h
    for i in range( numaxes ):#add axes from top down
        axesList.append( fig.add_axes( [axleft, bottom, axwidth, h] ) )
        bottom = bottom - (margin + h)
    
    return axesList

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

def sortGenes(genes, sample, genetype):
    if genetype not in ['v', 'j']:
        raise KeyError("genetype must be either 'v' or 'j'\n")
    sortedgenes = []
    gc = []
    if genetype == 'v':
        for i,g in enumerate(sample.vgenes):
            gc.append( (g, sample.vcounts[i]) )
    else:
        for i,g in enumerate(sample.jgenes):
            gc.append( (g, sample.jcounts[i]) )
    sortedgc = sorted(gc, key=lambda item:item[1], reverse=True)
    for item in sortedgc:
        if item[0] in genes:
            sortedgenes.append( item[0] )
    return sortedgenes

def getGeneUsage( sample, genetype, genes ):
    if genetype == "v":
        return sample.getVfreqs2(genes)
    elif genetype == "j":
        return sample.getJfreqs2(genes)
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

def drawUsageData( axesList, samples, samplesPerPlot, options, genetype ):
    if len( samples ) <= 0:
        return
    colors = getColors3()
    markers=['o', 'o', 'o', 'o', 'o', 'o', '^', '^', 'o', 'o','o','o']

    c = -2
    textsize = 'small'
    fontP = FontProperties()
    fontP.set_size(textsize)
    linesDict = {}
    labelsDict = {}

    xtickLabels = getIntersectGenes(samples, genetype)
    xtickLabels = sortGenes(xtickLabels, samples[0], genetype)
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
            sampleNames.append( "%s-%d" % (sample.name, sample.totalCount))
            ydata = getGeneUsage( sample, genetype, xtickLabels )
            maxy = max( [maxy, max(ydata)] )

            c += 2
            l = axes.plot( xdata, ydata, color=colors[c], marker=markers[c/2], markersize=6.0, markeredgecolor=colors[c], linestyle='none' )
            axes.plot( xdata, ydata, color=colors[c +1], linestyle='-', linewidth=0.01 )
            lines.append( l )

        linesDict[i] = lines
        labelsDict[i] = sampleNames
        axes.set_title( '%s Usage Distribution' % genetype.upper() )

    for i in range( len(axesList) ):
        axes = axesList[ i ]
        immunoseqLib.editSpine( axes )
        axes.set_xlabel('Genes', size = textsize)
        axes.set_ylabel('Percentage of reads', size=textsize)

        #Legend
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
        axes.set_ylim(-0.05, maxy + 0.05)
        #axes.set_ylim(-0.01, 0.25)
        for label in axes.get_xticklabels():
            label.set_fontsize( textsize )
            label.set_rotation( 90 )
        for label in axes.get_yticklabels():
            label.set_fontsize( textsize )
        axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
        axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    return 

def drawUsageDist( samples, options, genetype ):
    #Draw V usage distribution of each sample onto the same pdf
    options.out = os.path.join( options.outdir, "%sUsage" %genetype)
    fig, pdf = immunoseqLib.initImage( 10.0, 12.0, options )
    axesList = setUsageDistAxes( fig, len(samples), options.samplesPerPlot )
    drawUsageData( axesList, samples, options.samplesPerPlot, options, genetype )
    immunoseqLib.writeImage( fig, pdf, options )

#=======DRAW VJ, INCOMPLETE==========
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

def drawVJdata( axesList, sample, options ):
    #Draw heatmap:
    p0axes = axesList[0]
    data = sample.normalizeVJusage()
    hmaxes = p0axes.imshow( data, interpolation='nearest' )
   
    immunoseqLib.editSpine( p0axes )
    p0axes.set_title( "VJ usage of sample %s" %sample.name )
    xticks = [ x + 0.5 for x in range( len(sample.jgenes) ) ]
    p0axes.xaxis.set_ticks( xticks )
    p0axes.xaxis.set_ticklabels( sample.jgenes )
    
    yticks = [ x + 0.5 for y in range( len(sample.vgenes) ) ]
    p0axes.yaxis.set_ticks( yticks )
    p0axes.yaxis.set_ticklabels( sample.vgenes )
    
    textsize = 'x-small'
    for label in p0axes.get_xticklabels():
        label.set_fontsize( textsize )
        label.set_rotation( 70 )
    for label in p0axes.get_yticklabels():
        label.set_fontsize( textsize )

    #Draw Jusage histogram:
    #p1axes = axesList[1]
    #p1axes.bar( range( len(sample.jcounts) ), sample.jcounts )

    #Draw Vusage histogram:
    #p2axes = axesList[2]
    #bottoms = range( len(sample.vgenes) )
    #p2axes.barh( bottoms, sample.vcounts )




def drawVJ( sample, options ):
    options.out = os.path.join( options.outdir, "%s-vjUsage" %sample.name )
    fig, pdf = immunoseqLib.initImage( 8.0, 10.0, options )
    axesList = setCompareAxes( fig )

    drawVJdata( axesList, sample, options )
    immunoseqLib.writeImage( fig, pdf, options )


def readfiles( indir, noheader ):
    samples = []
    for file in os.listdir( indir ):
        if os.path.isdir(file):
            continue
        #name = os.path.basename( file )
        #name = file.split('-')[0]
        filenameItems = file.split('-')
        name = '-'.join( filenameItems[:len(filenameItems) -1] )
        if name == '':
            name = file.split('.')[0]
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
        samples.append( sample )

    return samples

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


def checkOptions( args, options, parser ):
    if not options.indir or not os.path.isdir( options.indir ):
        parser.error("Input directory is required. Please specify a valid input directory.\n")
    if not os.path.isdir( options.outdir ):
        parser.error("Output directory '%s' does not exist.\n" %options.outdir)
    if options.sampleOrder:
        options.sampleOrder = options.sampleOrder.split(',')

def initOptions( parser ):
    parser.add_option('-i', '--indir', dest='indir', help="Required argument: Input directory")
    parser.add_option('-o', '--outdir', dest='outdir', default='.', help="Output directory. Default = current directory")
    parser.add_option('--noheader', dest='noheader', action="store_true", default=False, help="If NOT specified (default), it is assumed that the first row has the matrix column names (J genes), and first column has the row names (vgenes).")
    parser.add_option('--samplesPerPlot', dest='samplesPerPlot', type='int', default=10)
    parser.add_option('--sampleOrder', dest='sampleOrder', type='string', help='Comma seperated list of samples in the order of interest')
    #parser.add_option('--vjPlotPerPage', dest='vjPlotPerPage', type='int', default=8)

def main():
    parser = immunoseqLib.initOptions()
    initOptions( parser )
    immunoseqLib.initPlotOptions( parser )
    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    immunoseqLib.checkPlotOptions( options, parser )

    samples = readfiles( options.indir, options.noheader )
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
    drawUsageDist( sortedSamples, options, 'v' )
    drawUsageDist( sortedSamples, options, 'j' )
    #for sample in samples:
    #    drawVJ( sample, options )

if __name__ == "__main__":
    main()
