#!/usr/bin/env python2.6

'''
Mon Jul  9 13:08:05 PDT 2012
nknguyen soe ucsc edu
'''
import os, re, sys
import immunoseq.lib.immunoseqLib as iseqlib
import matplotlib.pyplot as pyplot
import matplotlib.pylab as pylab
from matplotlib.ticker import *
from matplotlib.font_manager import FontProperties

####  colors ####
def sample2colorFixed():
    colors = ["#E31A1C", "#FE8E8F", "#377EB8", "#A6D7FE", "#4DAF4A", "#B8FEB5", "#984EA3", "#F6BDFE",
              "#FF7F00", "#FEBF80", 
              "#525252", "#737373", "#969696", "#BDBDBD",
              "#2B8CBE", "#7BCCC4", "#BAE4BC",
              "#225EA8", "#41B6C4", "#A1DAB4",
              "#253494",
              "#993404", "#D95F0E", 
              "#6A51A3", "#9E9AC8", "#CBC9E2",
              "#80FE80"
             ]
    s2i = {
           'adapt10R':0, 'adapt11R': 2, 'adapt12R': 4, 'adapt13R': 6, 'adapt1R': 8, 'adaptBR': 10,
           'adapt8D': 0, 'adapt11D': 2, 'adapt15D': 4, 'adapt16D': 6, 'adapt1D': 8, 'adaptBD': 10, 'adapt20D': 23, 'adapt1Ddraw2': 9, 'adaptBDdraw2': 11, 
           'adapt1Ddraw2notcd8': 22, 
           'as11D':0, 'as16D':2, 'as1D':4, 'asBD':26, 'as20D':18, 'as8D':8, 'as15D':6, 'adaptF57':10, 'adaptMA':12,
           #'as10R':0, 'as11R': 1, 'as12R': 2, 'as13R': 3, 'as1R': 4, 'as1D': 5, 'as8D': 6, 'as15D': 7,
           #'asBR': 8, 'asBD': 9,
           #'as11D':12, 'as16D':13, 'as1Ddraw2':21, 'as1Ddraw2notcd8':23,
           #'asBDdraw2':24, 'as20D':25,
           #'adaptBSD1': 10, 'adaptBSD1cd8n': 11, 'adaptBSD2R': 12, 'adaptBSD3': 13,
           #'adaptF28cd4': 14, 'adaptF57cd4': 15, 'adaptM35cd4': 16, 
           #'adaptF28cd8': 17, 'adaptF57cd8': 18, 'adaptM35cd8': 19, 
           #'adaptF28': 14, 'adaptF57': 15, 'adaptM35': 16, 
           #'adaptRep': 20,
           #'wangAll': 21, 'wangCyt': 22, 
           #'warrenM1': 23, 'warrenM2': 24, 'warrenF': 25,
           #'f28CD4memory': 0, 'f28CD4naive': 1, 'f28CD8memory': 2, 'f28CD8naive': 3, 
           #'f57CD4memory': 4, 'f57CD4naive': 5, 'f57CD8memory': 6, 'f57CD8naive': 7, 
           #'m35CD4memory': 8, 'm35CD4naive': 9, 'm35CD8memory': 10, 'm35CD8naive': 11,
          }
    s2c = {}
    for s, i in s2i.iteritems():
        s2c[s] = colors[i]
    return s2c

def sample2color(names):
    s2cfixed = sample2colorFixed()
    sample2color = {}

    #colors = iseqlib.getColors0()
    colors = iseqlib.getColors0()
    for i, name in enumerate(names):
        if name in s2cfixed:
            sample2color[name] = s2cfixed[name]
        else:
            sample2color[name] = colors[i]
    return sample2color 

def drawAll(options, outdir, rowname2cells, index2colname):
    if options.infile == '-':
        outname = 'all'
    else:
        outname = os.path.basename(options.infile).split('.')[0]
    options.out = os.path.join(outdir, outname)
    fig, pdf = iseqlib.initImage(10.0, 8.0, options)
    axes = fig.add_axes( [0.12, 0.15, 0.85, 0.75] )

    lines = []
    #rownames = sorted( rowname2cells.keys() )
    rownames = ['as11D', 'as16D', 'as1D', 'asBD', 'as20D', 'as15D', 'as8D']
    name2color = sample2color(rownames)
    
    #xdata = sorted( [ int(x) for x in colname2index.keys() ] )
    xmax = 0.0
    xmin = float('inf')
    #ymax = 0
    markersize = 12.0
    xindices = []
    #for rowname in rownames:
    r = 0
    while r < len(rownames):
        rowname = rownames[r]
        row = rowname2cells[rowname]
        means = []
        stds = []
        xdata = []
        for i, m in enumerate( row ):
            if i % 2 == 0 and m != 'NA' and m != '' and m != '-':
                colname = index2colname[i]
                try:
                    colname = int(colname)
                    xdata.append( colname )
                except:
                    xdata.append(i/2)
                    if i not in xindices:
                        xindices.append(i)
                means.append( float(m) )

            elif i%2 == 1 and m != 'NA' and m != '' and m != '-':
                stds.append( float(m) )
        xmax = max([xmax, max(xdata)])
        xmin = min([xmin, min(xdata)])
         
        #HACK
        #if rowname == 'uniqClones':
        #if rowname == 'mountford' or rowname == 'horn':
        #if rowname != 'horn':
        #if rowname == 'manhattan' or rowname == 'euclidean' or rowname == 'binomial':
        exceptions = ['manhattan', 'euclidean', 'binomial', 'kulczynski', 'canberra', 'jaccard']
        if rowname in exceptions:
            rownames.remove(rowname)
            continue
        #means = [m/min(means) for m in means]
        #stds = [0.0 for s in stds]
        #END HACK
        color = name2color[rowname]
        axes.errorbar(xdata, means, yerr=stds, color=color, markeredgecolor=color, markersize=markersize, fmt='.')
        line = axes.plot(xdata, means, color=color, linestyle='-', linewidth=4.0)
        lines.append(line)
        r += 1
    
    #axes.set_xscale('log')
    #axes.set_yscale('log')
    iseqlib.editSpine(axes)
    #axes.set_title("%s index across different sampling sizes" %outname, size='xx-large')
    #axes.set_xlabel("Sampling size (number of reads)", size='large' )
    #axes.set_ylabel("%s index" %outname, size='large')
    
    #HACK
    axes.set_title("Sequencing Saturation", size='xx-large', weight='bold')
    axes.set_xlabel("Sampling size (number of sequences, in millions) ", size='x-large', weight='bold' )
    #axes.set_xlabel("Sampling size (number of sequences, in thousands) ", size='x-large', weight='bold' )
    axes.set_ylabel("Number of clones (in thousands)", size='x-large', weight='bold')
    
    fontP = FontProperties()
    fontP.set_size('medium')
    #rownames = ["A", "B"]
    rownames = [iseqlib.properName(n) for n in rownames]
    legend = axes.legend( lines, rownames, numpoints = 1, loc='best', ncol = 1, prop=fontP)
    legend._drawFrame = False

    if len(xindices) > 0:
        xticklabels = []
        for i in sorted(xindices):
            xticklabels.append( index2colname[i] )
        axes.xaxis.set_ticks( xrange(len(xticklabels)) )
        axes.xaxis.set_ticklabels( xticklabels )

    #HACK:
    #xticks = [ 10000, 50000, 100000, 200000, 300000, 400000, 500000]
    #xticklabels = [ str(x/1000) for x in xticks]
    xticks = [0, 1000000,2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000]
    xticklabels = [ str(x) for x in xrange(0, 10) ]
    axes.xaxis.set_ticks(xticks)
    axes.xaxis.set_ticklabels( xticklabels )
    
    #yticks = xrange(0, 121000, 20000)
    #yticklabels = [ str(y) for y in xrange(0, 121, 20) ]
    yticks = xrange(0, 250000, 50000)
    yticklabels = [ str(y) for y in xrange(0, 250, 50) ]
    axes.yaxis.set_ticks(yticks)
    axes.yaxis.set_ticklabels(yticklabels)

    
    for label in axes.get_xticklabels():
        label.set_fontsize('large')
        label.set_fontweight ('bold')
        #label.set_rotation(45)
    for label in axes.get_yticklabels():
        label.set_fontsize('large')
        label.set_fontweight ('bold')
        
    axes.xaxis.grid(b=True, color='#3F3F3F', linestyle='-', linewidth=0.05)
    axes.yaxis.grid(b=True, color='#3F3F3F', linestyle='-', linewidth=0.05)
    #HACK
    #axes.set_ylim(0.996, 1)
    xspan = xmax - xmin
    #axes.set_xlim(xmin - xspan*0.01, xmax + xspan*0.01)
    #axes.set_ylim(10000, 121000)
    
    #axes.set_xlim(-10, 110000)
    #axes.set_ylim(-10, 60000)
    
    axes.set_xlim(-10, 5100000)
    axes.set_ylim(-10, 205000)
    
    iseqlib.writeImage(fig, pdf, options)

def drawPlots(options, outdir, rowname2cells, index2colname):
    drawAll(options, outdir, rowname2cells, index2colname)

#================ LATEX TABLE ====================#
def tabHeader(f, colnames):
    f.write("\\begin{table}\n")
    f.write("\\centering\n")
    f.write("\\scalebox{1.0}{%\n")
    f.write("\\begin{tabular}{c%s}\n" %("|c"*len(colnames)))
    f.write("\\hline\n")
    f.write("Sample & %s \\\\\n" %(" & ".join(colnames)))
    f.write("\\hline\n")

def tab(f, rowname2cells, index2colname):
    for rowname in sorted( rowname2cells.keys() ):
        row = rowname2cells[rowname]
        f.write("%s" %rowname)
        for i, m in enumerate( row ):
            if i % 2 == 0 and m != 'NA' and m != '' and m != '-':
                colname = index2colname[i]
                if colname in ['uniqClones', 'fisherAlpha']:
                    f.write( " & %s " % iseqlib.prettyInt( int(float(m)) ) )
                else:
                    f.write( " & %.3f " % float(m) )
            elif i%2 == 1 and m != 'NA' and m != '' and m != '-':
                colname = index2colname[i -1]
                if colname in ['uniqClones', 'fisherAlpha']:
                    f.write( " $\\pm$ %s " % iseqlib.prettyInt(int(float(m))) )
                else:
                    f.write( " $\\pm$ %.3f " % float(m) )
        f.write("\\\\\n")
        f.write("\\hline\n")

def getLatexTab(outfile, rowname2cells, index2colname):
    f = open(outfile, 'w')
    iseqlib.writeDocumentStart(f)
    colnames = []
    for i in sorted( index2colname.keys() ):
        colnames.append( index2colname[i] )
    tabHeader(f, colnames)
    tab(f, rowname2cells, index2colname)
    label = ''
    captionStr = ''
    iseqlib.tableCloser(f, captionStr, label)
    iseqlib.writeDocumentEnd(f)
    f.close()

#============== END LATEX TABLE ==================#


def readfile(infile):
    index2colname = {} #key = column index, val = column name
    name2row = {} #key = rowname, vals = list of row's cells
    infh = sys.stdin
    if infile != '-':
        infh = open(infile, 'r')
    
    firstRow = True
    for line in infh:
        items = line.strip().split('\t')
        if firstRow:
            for i in xrange(1, len(items), 2):
                index2colname[ i -1 ] = items[i]
            firstRow = False
        else:
            rowname = items[0]
            name2row[rowname] = items[1:]

    if infile != '-':
        infh.close()
    return name2row, index2colname

def main():
    parser = iseqlib.initOptions()
    parser.add_option('-i', '--infile', dest='infile', default='-', help='Input file. Default is stdin')
    parser.add_option('-o', '--outdir', dest='outdir', default='.', help='Output directory. Default is current directory')
    parser.add_option('--noplot', dest='noplot', action='store_true', default=False, help='If specified, do not draw plot')
    parser.add_option('-t', '--table', dest='table', action='store_true', default=False, help='If specified, make latex table')
    #parser.add_option('-o', '--outfile', dest='outfile', default='-', help='Output file. Default is stdout')

    iseqlib.initPlotOptions(parser)
    options, args = parser.parse_args()
    iseqlib.checkPlotOptions( options, parser )

    rowname2row, index2colname = readfile(options.infile)
    if not options.noplot:
        drawPlots(options, options.outdir, rowname2row, index2colname)

    if options.table:
        outfile = os.path.join(options.outdir, "%s.tex" %( os.path.basename(options.infile).split('.')[0] ))
        getLatexTab(outfile, rowname2row, index2colname)

if __name__ == '__main__':
    main()


