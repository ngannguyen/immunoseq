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
    f.write("\\begin{table}\n")
    f.write("\\centering\n")
    f.write("\\scalebox{0.9}{%\n")
    f.write("\\begin{tabular}{c%s}\n" %( "|c"*(len(colnames)*2) ))
    f.write("\\hline\n")
    f.write(" \\multirow{2}{*}{Sample} & %s \\\\\n" %(" & ".join( ["\\multicolumn{2}{c|}{%s}" % col.upper() for col in colnames] )))
    f.write("\\cline{2-%d}\n" %( len(colnames)*2 + 1 ))
    #f.write("\\hline\n")
    f.write(" %s \\\\\n" %("& Count & \% "*(len(colnames)) ))
    f.write("\\hline\n")

def tab(f, colnames, sample2row):
    #Print total:
    totalRow = sample2row['Total']
    f.write("%s & %s \\\\\n" % ("Total", " & ".join( ["%s & 100.00" %cell for cell in totalRow] )) )
    f.write("\\hline\n")
    for sample in sorted(sample2row.keys()):
        if sample == 'Total' or sample == 'controls' or sample == 'patients' :
            continue
        f.write("%s" % iseqlib.properName(sample))
        row = sample2row[sample]
        for i, cell in enumerate(row):
            total = totalRow[i]
            pc = iseqlib.getPc( int(cell), int(total) )
            f.write( " & %s & %.2f " %(cell, pc) )
        f.write("\\\\\n")
        f.write("\\hline\n")

def getTab(outfile, colnames, sample2row):
    f = open(outfile, 'w')
    iseqlib.writeDocumentStart(f)
    myTabHeader(f, colnames[1:])
    tab(f, colnames, sample2row)
    label = ''
    captionStr = ''
    iseqlib.tableCloser(f, captionStr, label)
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

