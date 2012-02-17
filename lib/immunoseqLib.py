#nknguyen soe ucsc edu
#Library of common functions used in the immunoseq pipeline

from numpy import *
from matplotlib.ticker import LogLocator
from optparse import OptionParser

#============== UTILITIES =========
def sortDictByValue( dictionary ):
    items = [ (v, k) for k,v in dictionary.items() ]
    items.sort()
    items.reverse()
    return [ (k, v) for v, k in items ]

#============== LATEX =========
def prettyInt( number ):
    numstr = str( number )
    prettyStr = ''
    for i in range( len(numstr) ):
        if i > 0 and (len(numstr) - i) %3 == 0:
            prettyStr = prettyStr + ','
        prettyStr = prettyStr + numstr[i]

    return prettyStr

def prettyFloat( num ):
    if num == 0:
        return "0"
    else:
        return "%.2e" %num

def writeDocumentStart( f ):
    f.write( "\\documentclass[11pt]{article}\n" )

    f.write( "\\usepackage{epsfig}\n" )
    f.write( "\\usepackage{multirow}\n" )
    f.write( "\\usepackage{graphicx}\n" )
    f.write( "\\usepackage{array}\n" )
    f.write( "\\usepackage{color, colortbl}\n" )
    f.write( "\\usepackage{rotating}\n" )
    f.write( "\\usepackage[table]{xcolor}\n" )
    f.write( "\\definecolor{LightCyan}{rgb}{0.88,1,1}")
    f.write( "\n" )

    f.write( "\\newcommand{\\figref}[1]{Figure~\\ref{fig:#1}}\n" )
    f.write( "\\newcommand{\\tabref}[1]{Table~\\ref{tab:#1}}\n" )
    f.write( "\n" )

    f.write("\\textwidth=6.5in\n")
    f.write("\\textheight=9in\n")
    f.write("\\oddsidemargin=0in\n")
    f.write("\\evensidemargin=0in\n")
    f.write("\\topmargin=0in\n")
    f.write("\\topskip=0in\n")
    f.write("\\headheight=0in\n")
    f.write("\\headsep=0in\n")
    f.write("\n")

    f.write("\\begin{document}\n")
    return

def writeDocumentEnd( f ):
    f.write( "\\end{document}\n" )

def tableCloser(f, captionStr, label):
    f.write("\\end{tabular}\n")
    f.write("}\n")
    f.write("\\caption{%s}\n" %captionStr)
    #f.write("\\end{center}\n")
    f.write("\\label{%s}" %label)
    f.write("\\end{table}\n\n")

#======= ERROR classes ===========
class FileFormatError( Exception ):
    pass

class InputOptionError( Exception ):
    pass

#======= initialize options =======
def initOptions():
    usage = "usage: %prog [options]\n"
    parser = OptionParser( usage=usage )
    return parser

#======= matplotlib functions ============
def initPlotOptions( parser ):
    parser.add_option( '--dpi', dest='dpi', default=300, type='int', help='Dots per inch. Default=%default')
    parser.add_option( '--outFormat', dest='outFormat', default='pdf', type='string',
                             help='Output plot format [pdf|png|all|eps]. default=%default' )

def checkPlotOptions( options, parser ):
    if options.dpi < 72:
        parser.error('dpi must >= than screen res, 72. Got (%d)' % options.dpi )
    if options.outFormat not in ('pdf', 'png', 'eps', 'all'):
        parser.error('Unrecognized output format: %s. Choose one from: pdf png eps all.' % options.outFormat )
    options.out = 'plot'

def initImage( width, height, options ):
    """
    initImage takes a width and height and returns
    both a fig and pdf object. options must contain outFormat,
    and dpi
    """
    import matplotlib.backends.backend_pdf as pltBack
    import matplotlib.pyplot as plt
    pdf = None
    if options.outFormat == 'pdf' or options.outFormat == 'all':
        pdf = pltBack.PdfPages( options.out + '.pdf' )
    fig = plt.figure( figsize=( width, height ), dpi=options.dpi, facecolor='w' )
    return ( fig, pdf )

def writeImage( fig, pdf, options ):
    """
    writeImage assumes options contains outFormat and dpi.
    """
    if options.outFormat == 'pdf':
        fig.savefig( pdf, format='pdf' )
        pdf.close()
    elif options.outFormat == 'png':
        fig.savefig( options.out + '.png', format='png', dpi=options.dpi )
    elif options.outFormat == 'all':
        fig.savefig( pdf, format='pdf' )
        pdf.close()
        fig.savefig( options.out + '.png', format='png', dpi=options.dpi )
        fig.savefig( options.out + '.eps', format='eps' )
    elif options.outFormat == 'eps':
        fig.savefig( options.out + '.eps', format='eps' )

def getColors0():
    colors = [ "#1f77b4", "#aec7e8", # blues 
               "#9467bd", "#c5b0d5", # lavenders
               "#ff7f0e", "#ffbb78", # oranges
               "#2ca02c", "#98df8a", # greens
               "#d62728", "#ff9896", # reds
               "#4A766E", "#48D1CC", #
               "#DAA520", "#DBDB70", #gold
               "#302B54", "#5959AB", #purple - blueish
               "#F6A8B6", "#F6CCDA"] #pink

    return colors

#AS color:
def getColors3():
    colors = [ "#FF3333", 
               "#00008B", 
               "#D02090", 
               "#FF6EB4", 
               "#EEAD0E", 
               "#CD2626", 
               "#FF7D40", 
               "#0000FF" ]
    return colors

def setAxes( fig ):                                                                     
    return fig.add_axes( [0.12, 0.1, 0.83, 0.85] )                                      
                                                                                        
def editSpine( axes ):                                                                 
    for loc, spine in axes.spines.iteritems():                                          
        if loc in [ 'left', 'bottom' ]:                                                 
            spine.set_position( ('outward', 10) )                                       
        elif loc in ['right', 'top']:                                                   
            spine.set_color('none')                                                     
        else:                                                                           
            raise ValueError( 'Unknown spine location %s\n' % loc )                     
                                                                                        
def setTicks( axes ):                                                                   
    axes.xaxis.set_ticks_position( 'bottom' )                                           
    axes.yaxis.set_ticks_position( 'left' )                                             
    minorLocator = LogLocator( base=10, subs = range(1, 10) )                           
    axes.xaxis.set_minor_locator( minorLocator )                                        
                                                
def bihist( y1, y2, axes, bins, orientation, color=None ):
    n1, bins1, patch1 = axes.hist( y1, bins=bins, orientation=orientation, color=color ) #Top hist
    n2, bins2, patch2 = axes.hist( y2, bins=bins, orientation=orientation, color=color ) #bottom hist
    #set ymax:
    if orientation == 'vertical':
        ymax = max( [ i.get_height() for i in patch1 ] )
        for i in patch2:
            i.set_height( -i.get_height() )
        ymin = min( [ i.get_height() for i in patch2 ] )
    elif orientation == 'horizontal':
        ymax = max( [ i.get_width() for i in patch1 ] )
        for i in patch2:
            i.set_width( -i.get_width() )
        ymin = min( [ i.get_width() for i in patch2 ] )
    
    #axes.set_ylim( ymin*1.1, ymax*1.1 )
    return ymin, ymax

