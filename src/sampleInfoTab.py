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

def myTabHeader(f, colnames):
    f.write("\\begin{sidewaystable}\n")
    f.write("\\centering\n")
    f.write("\\scalebox{0.9}{%\n")
    f.write("\\begin{tabular}{%s}\n" %("|".join(["c" for c in colnames]) ))
    f.write("\\hline\n")
    f.write("\\textbf{%s}\\\\\n" %("} & \\textbf{".join( [c.replace("_", " ") for c in colnames] )))
    f.write("\\hline\n")

def tab(f, sample2row):
    #Print total:
    for sample in sorted(sample2row.keys()):
        f.write("%s" % iseqlib.properName(sample))
        row = sample2row[sample]
        for i, cell in enumerate(row):
            if cell in ['+', '-']:
                cell = "$%s$" %cell
            f.write( " & %s " %(cell) )
        f.write("\\\\\n")
        f.write("\\hline\n")

def getTab(outfile, colnames, sample2row):
    f = open(outfile, 'w')
    iseqlib.writeDocumentStart(f)
    myTabHeader(f, colnames)
    tab(f, sample2row)
    label = ''
    captionStr = ''
    iseqlib.sidewaystableCloser(f, captionStr, label)
    iseqlib.writeDocumentEnd(f)
    f.close()

def main():
    parser = iseqlib.initOptions()
    parser.add_option('-i', '--infile', dest = 'infile', help='Input file')
    parser.add_option('-o', '--outfile', dest= 'outfile', help='Output file')

    options, args = parser.parse_args()
    colnames, sample2row = readfile(options.infile)
    getTab(options.outfile, colnames, sample2row)

if __name__ == '__main__':
    main()

