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

class Sample:
    def __init__(self, name):
        self.name = name
        self.totalCount = 0
        self.totalClones = 0
        self.size2clones = {}
        self.simpson = 0

    def setTotal(self):
        self.totalCount = sum( [ s*c for s,c in self.size2clones.iteritems() ] )
        self.totalClones = sum( [ c for s,c in self.size2clones.iteritems() ] )

    def avrTotal(self, num):
        self.totalCount /= num
        self.totalClones /= num

    def calcSimpsonIndex(self):
        denominator = self.totalCount*(self.totalCount - 1)
        numerator = 0
        for size, clones in self.size2clones.iteritems():
            numerator += size*(size - 1)*clones
        
        self.simpson = 1 - float(numerator)/denominator

def readfiles(indir, cutoff):
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

                #Filter out sequences with size < cutoff:
                if count < cutoff:
                    continue
                
                if count not in s.size2clones:
                    s.size2clones[ count ] = 1
                else:
                    s.size2clones[ count ] += 1
        f.close()
        s.setTotal()
        s.calcSimpsonIndex()
        samples.append(s)
    return samples

def getAvrSamples(samples, group2samples):
    avrSamples = []
    for group, names in group2samples.iteritems():
        groupAvrSample = Sample(group)
        for name in names:
            for s in samples:
                if s.name == name:
                    groupAvrSample.simpson += s.simpson 
                    for size, clones in s.size2clones.iteritems():
                        if size not in groupAvrSample.size2clones:
                            groupAvrSample.size2clones[size] = clones
                        else:
                            groupAvrSample.size2clones[size] += clones
        #for size in groupAvrSample.size2clones:
        #    groupAvrSample.size2clones[size] /= len(names)
        groupAvrSample.setTotal()
        groupAvrSample.simpson /= len(names)
        #average
        #print group, names
        #print groupAvrSample.totalCount, groupAvrSample.totalClones
        groupAvrSample.avrTotal(len(names))
        #print groupAvrSample.totalCount, groupAvrSample.totalClones
        avrSamples.append(groupAvrSample)
    return avrSamples

def tabHeader(f, simpson):
    f.write("\\begin{table}\n") 
    f.write("\\centering\n")
    f.write("\\scalebox{1}{%\n")
    if simpson:
        f.write("\\begin{tabular}{c|r|r|r}\n")
        f.write("\\hline\n")
        f.write("\\cellcolor[gray]{0.9} Sample &\\cellcolor[gray]{0.9} Clones &\\cellcolor[gray]{0.9} Total Reads &\\cellcolor[gray]{0.9} Simpson Index\\\\\n")
    else:
        f.write("\\begin{tabular}{c|r|r}\n")
        f.write("\\hline\n")
        f.write("\\cellcolor[gray]{0.9} Sample &\\cellcolor[gray]{0.9} Clones &\\cellcolor[gray]{0.9} Total Reads \\\\\n")
    f.write("\\hline\n")

def sortSamplesByGroup(samples, groupsamples, group2samples):
    sortedSamples = []
    for g in sorted( group2samples.keys() ):
        gsamples = group2samples[g]
        gsamples.sort()
        for s in gsamples:
            for sample in samples:
                if sample.name == s:
                    sortedSamples.append(sample)
        for gs in groupsamples:
            if gs.name == g:
                sortedSamples.append(gs)
    return sortedSamples
    
def tab(f, samples, groupsamples, group2samples, simpson):
    if len(groupsamples) > 0:
        samples = sortSamplesByGroup(samples, groupsamples, group2samples)

    for s in samples:
        if s.name in group2samples:
            f.write("\\cellcolor[gray]{0.9} Avr %s & \\cellcolor[gray]{0.9} %s & \\cellcolor[gray]{0.9} %s" %(s.name.replace('_', "-"), iseqlib.prettyInt(s.totalClones), iseqlib.prettyInt(s.totalCount)))
            if simpson:
                f.write("& \\cellcolor[gray]{0.9} %.3f" % s.simpson)
            f.write("\\\\\n")
            f.write("\\hline\n")
        else:
            f.write("%s & %s & %s" %(s.name.replace('_', "-"), iseqlib.prettyInt(s.totalClones), iseqlib.prettyInt(s.totalCount)))
            if simpson:
                f.write("& %.3f" %s.simpson)
            f.write("\\\\\n")

    if len(groupsamples) == 0:
        f.write("\\hline\n")

def getClonesizeLatexTab(samples, groupsamples, group2samples, outfile, simpson): 
    f = open(outfile, 'w')
    iseqlib.writeDocumentStart(f)

    tabHeader(f, simpson)
    tab(f, samples, groupsamples, group2samples, simpson)

    label = ''

    captionStr = "TCRB sequence statistics. Columns: `Sample': sample name, `Clones': number of unique clones, `Total Reads': number of total reads. Rows: different samples, where the shaded rows show the average of statistics of each group."
    if simpson:
        captionStr = "TCRB sequence statistics. Columns: `Sample': sample name, `Clones': number of unique clones, `Total Reads': number of total reads, `Simpson Index': Simpson diversity index, which is equivalent with the probability that any two sequences (reads) belong to two different clones. Rows: different samples, where the shaded rows show the average of statistics of each group."
    iseqlib.tableCloser(f, captionStr, label)
    iseqlib.writeDocumentEnd(f)
    f.close()

def getClonesizeTextTab(samples, groupsamples, group2samples, outfile, simpson):
    if len(groupsamples) >0:
        samples = sortSamplesByGroup(samples, groupsamples, group2samples)
    else:
        samples = sorted(samples, key=lambda sample: sample.totalCount) 
    f = open(outfile, 'w')
    if simpson:
        f.write("Sample\tClones\tTotal Reads\tSimpson Index\n")
    else:
        f.write("Sample\tClones\tTotal Reads\n")
    for s in samples:
        if s.name in group2samples:
            f.write("Avr %s\t%s\t%s" %(s.name.replace('_', "-"), iseqlib.prettyInt(s.totalClones), iseqlib.prettyInt(s.totalCount)))
            if simpson:
                f.write("\t%.3f" %s.simpson)
            f.write("\n")
            f.write("#\n")
        else:
            f.write("%s\t%s\t%s" %(s.name.replace('_', "-"), iseqlib.prettyInt(s.totalClones), iseqlib.prettyInt(s.totalCount)))
            if simpson:
                f.write("\t%.3f" %s.simpson)
            f.write("\n")
    f.close()

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

def initOptions( parser ):
    #parser.add_option('-o', '--outdir', dest='outdir', default='.', help="Output directory. Default = current direcotory")
    parser.add_option('-o', '--outfile', dest='outfile', default='clonesize.tex', help="Output directory. Default = %default")
    parser.add_option('-c', '--cutoff', dest='cutoff', type='int', default=1, help="Minimum count. Default = %default")
    parser.add_option('-g', '--group', dest='group', help='Group2samples file (Format: group<space>commaSeparatedListOfSamples. If specified, only draw average of each group.')
    parser.add_option('-t', '--text', dest='text', action='store_true', default=False, help='If specified, output is a tab limited textfile instead of a latex file. Default = %default')
    parser.add_option('-s', '--simpsonIndex', dest='simpson', action='store_true', default=False, help='If specified, reports the Simpson Index of each sample')

def main():
    usage = ('Usage: %prog [options] inputDirectory\n')

    parser = OptionParser(usage = usage)
    initOptions( parser )
    options, args = parser.parse_args()
    checkOptions( args, options, parser )

    samples = readfiles( args[0], options.cutoff )
    groupsamples = []
    group2samples = {}
    if options.group:
        group2samples = readGroup2Samples(options.group) 
        groupsamples = getAvrSamples(samples, group2samples)
    if options.text:
        getClonesizeTextTab(samples, groupsamples, group2samples, options.outfile, options.simpson)
    else:
        getClonesizeLatexTab(samples, groupsamples, group2samples, options.outfile, options.simpson)
     
if __name__ == "__main__":
    main()
