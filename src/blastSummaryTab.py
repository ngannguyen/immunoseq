#!/usr/bin/env python

'''
'''
import os, sys, re
import immunoseq.lib.immunoseqLib as iseqlib

def readfile(file):
    f = open(file, 'r')
    f.readline()
    rows = []
    for line in f:
        items = line.split()
        
        group = items[0]
        numsam = items[1]
        total = items[2]
        autohit = items[5]
        pchit = "%.2f" % float( items[6] )
        #if group == 'AS':
        #    pchit = "\\textbf{%s}" %pchit
        row = [group, numsam, total, autohit, pchit]

        rows.append(row)
    f.close()
    return rows

def myTabHeader(f):
    f.write("\\begin{table}\n")
    f.write("\\centering\n")
    f.write("\\scalebox{0.9}{%\n")
    f.write("\\begin{tabular}{c|c|c|c|c}\n")
    f.write("\\hline\n")
    f.write(" \\textbf{Group} & \\textbf{Samples}  &  \\textbf{Total}  & \\textbf{Hits} & \\textbf{ \% Hits/Total } \\\\\n" )
    f.write("\\hline\n")

def tab(f, rows):
    for row in rows:
        f.write("%s\\\\\n" % " & ".join(row))
        f.write("\\hline\n")

def getTab(outfile, rows):
    f = open(outfile, 'w')
    iseqlib.writeDocumentStart(f)
    myTabHeader(f)
    tab(f, rows)
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
    rows = readfile(options.infile)
    getTab(options.outfile, rows)

if __name__ == '__main__':
    main()

