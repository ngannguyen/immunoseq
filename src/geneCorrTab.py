#!/usr/bin/env python

'''
'''
import os, sys, re
import immunoseq.lib.immunoseqLib as iseqlib

def readfile(file):
    f = open(file, 'r')
    colnames = f.readline().split('\t')
    sample2row = {}
    for line in f:
        items = line.split('\t')
        sample2row[items[0]] = items[1:]
    f.close()
    return colnames, sample2row

def readfiles(files):
    cat2colnames = {}
    cat2sample2row = {}
    for file in files:
        category = os.path.basename(file).split('.')[0]
        colnames, sample2row = readfile(file)
        cat2colnames[category] = colnames
        cat2sample2row[category] = sample2row
    return cat2colnames, cat2sample2row

def getRownames(cat2sample2row):
    return cat2sample2row.values()[0].keys()

def myTabHeader(f, categories):
    cols = ['P-value', 'AS, B27+', 'Healthy, B27+']
    numcats = len(categories)
    numcols = len(cols)*numcats

    f.write("\\begin{table}\n")
    f.write("\\centering\n")
    f.write("\\scalebox{1}{%\n")
    f.write("\\begin{tabular}{c%s}\n" %( ("|%s" %("|c"*len(cols)) )*numcats ))
    f.write("\\hline\n")
    f.write(" \\multirow{2}{*}{Gene} & %s \\\\\n" %(" & ".join( ["\\multicolumn{%d}{c||}{%s}" % (len(cols), c) for c in categories] )))
    #f.write(" \\cellcolor[gray]{0.9} \\multirow{2}{*}{Gene} & %s \\\\\n" %(" & ".join( ["\\cellcolor[gray]{0.9} \\multicolumn{%d}{c||}{%s}" % (len(cols), c) for c in categories] )))
    f.write("\\cline{2-%d}\n" %( numcols + 1 ))
    for i in xrange(numcats):
        f.write(" & %s" % " & ".join(cols))
    f.write("\\\\\n")
    f.write("\\hline\n")

def tab(f, cat2colnames, cat2sample2row):
    #rownames = sorted( getRownames(cat2sample2row) )
    rownames = ['v', 'vj', 'vdj', 'd', 'j', 'dj']
    cats = sorted( cat2sample2row.keys() )
    for rowname in rownames:
        #grayrows = ['v', 'vj', 'vdj']
        grayrows = []
        if rowname in grayrows:
            f.write("\\cellcolor[gray]{0.9} %s" % rowname.upper())
        else:
            f.write("%s" % rowname.upper())
        
        for c in cats:
            row = cat2sample2row[c][rowname]
            if rowname in grayrows:
                f.write(" & \\cellcolor[gray]{0.9} %.3f " % float(row[0]))
            else:
                f.write(" & %.3f " % float(row[0]))

            for i in xrange(1, len(row) - 1, 2):
                avr = float(row[i])
                std = float(row[i+1])
                if rowname in grayrows:
                    f.write( " & \\cellcolor[gray]{0.9} %.3f $\\pm$ %.3f " %(avr, std) )
                else:
                    f.write( " & %.3f $\\pm$ %.3f " %(avr, std) )
        f.write("\\\\\n")
        f.write("\\hline\n")

def getTab(outfile, cat2colnames, cat2sample2row):
    f = open(outfile, 'w')
    iseqlib.writeDocumentStart(f)
    myTabHeader(f, cat2colnames.keys())
    tab(f, cat2colnames, cat2sample2row)
    label = ''
    captionStr = ''
    iseqlib.tableCloser(f, captionStr, label)
    iseqlib.writeDocumentEnd(f)
    f.close()

def main():
    parser = iseqlib.initOptions()
    parser.add_option('-i', '--infiles', dest = 'infiles', help='Input files')
    parser.add_option('-o', '--outfile', dest= 'outfile', help='Output file')

    options, args = parser.parse_args()
    cat2colnames, cat2sample2row = readfiles(options.infiles.split(','))
    getTab(options.outfile, cat2colnames, cat2sample2row)

if __name__ == '__main__':
    main()

