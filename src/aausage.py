#!/usr/bin/env python

#
#nknguyen soe ucsc edu
#Apr 26 2011
#Profile the aa usage of input fasta files
#
#Input: directory containing input fasta files
#Output: For each sequence length, returns a matrix whose columns are the 
#samples, and rows are usage of different amino acid at different positions.
#
#

import sys, os, re, time
from optparse import OptionParser
from sonLib.bioio import system
from immunoseq.lib.immunoseqLib import *
import numpy as np

def setLen2seqs(samples):
    for sample in samples:
        len2seqs = {}
        for header, seq in sample.seqs.iteritems():
            l = len(seq.seq)
            if l not in len2seqs:
                len2seqs[l] = {header:seq}
            else:
                len2seqs[l][header] = seq
        sample.len2seqs = len2seqs
    return

def aaSampleProfile(seqs, length, sizes):
    size2pos2aa2usage = {}
    totalCount = len(seqs)
    
    #sizes = [1, 2, 3, 4] #single-aa, di-aa, tri-aa, 4-aa
    
    for size in sizes:
        size2pos2aa2usage[size] = {}
        for i in xrange(0, length - size + 1):
            aa2usage = {}
	    for seq in seqs.values():
                aa = seq.seq[i:i+size]
                if aa not in aa2usage:
                    aa2usage[aa] = 1.0
                else:
                    aa2usage[aa] += 1.0
            #Normalize:
            for aa in aa2usage:
                aa2usage[aa] /= float(totalCount)
            size2pos2aa2usage[size][ '-'.join( [str(p) for p in xrange(i, i+size)] ) ] = aa2usage 
    
    return size2pos2aa2usage

def aaProfilePerLen(samples, length, verbose):
    size2pos2aa2sample2usage = {}
    sizes = [1, 2, 3, 4] #single-aa, di-aa, tri-aa, 4-aa
    for sample in samples:
        stime = time.clock()
        if verbose:
            sys.stderr.write('Get aaUsage for length %d, sample %s\n' %(length, sample.name))

        if length not in sample.len2seqs:
            continue

        seqs = sample.len2seqs[length]
        size2pos2aa2usage =  aaSampleProfile(seqs, length, sizes)

        if verbose:
            sys.stderr.write('Done in %d seconds.\n' %(time.clock() - stime) )
        stime = time.clock()

        #add to pos2aa2sample2usage
        for size, pos2aa2usage in size2pos2aa2usage.iteritems():
            pos2aa2sample2usage = {}
            for pos, aa2usage in pos2aa2usage.iteritems():
                if pos not in pos2aa2sample2usage:
                    pos2aa2sample2usage[pos] = {}
                    for aa, usage in aa2usage.iteritems():
                        pos2aa2sample2usage[pos][aa] = {sample.name: usage}
                else:
                    for aa, usage in aa2usage.iteritems():
                        if aa not in pos2aa2sample2usage[pos]:
                            pos2aa2sample2usage[pos][aa] = {sample.name: usage}
                        else:
                            pos2aa2sample2usage[pos][aa][sample.name] = usage
            size2pos2aa2sample2usage[size] = pos2aa2sample2usage
        if verbose:
            sys.stderr.write('add to hash in %d seconds.\n' %(time.clock() - stime))

    return size2pos2aa2sample2usage

def printProfile(size2pos2aa2sample2usage, sampleNames, outdir):
    for size, pos2aa2sample2usage in size2pos2aa2sample2usage.iteritems():
        outfile = os.path.join(outdir, "%d.txt" %size)
        f = open(outfile, 'w')
        f.write("#Pos-AA\t%s\n" %( "\t".join(sampleNames) ))
        for pos, aa2sample2usage in pos2aa2sample2usage.iteritems():
            for aa, sample2usage in aa2sample2usage.iteritems():
                rowname = "%s_%s" %(pos, aa)
                usage = []
                for sample in sampleNames:
                    if sample not in sample2usage:
                        usage.append('0')
                    else:
                        usage.append( str( sample2usage[sample] ) )
                f.write("%s\t%s\n" %(rowname, "\t".join(usage) ) )
        f.close()

def getAaUsage(samples, outdir, verbose):
    lens = xrange(11, 17)
    sampleNames = sorted( [s.name for s in samples] )

    for l in lens:
        stime = time.clock()
        if verbose:
            sys.stderr.write("Getting aaUsage for length %d.\n" %l)

        #outfile = os.path.join(outdir, "%d.txt" %l)
        lendir = os.path.join(outdir, "%d" %l)
        system("mkdir -p %s" % lendir)
        size2pos2aa2sample2usage = aaProfilePerLen(samples, l, verbose)
        
        if verbose:
            sys.stderr.write("Got usage in %d seconds.\n" %(time.clock() - stime) )
        stime = time.clock()

        printProfile(size2pos2aa2sample2usage, sampleNames, lendir)
    
        if verbose:
            sys.stderr.write("Printed usage in %d seconds.\n" %(time.clock() -stime) )
    return

#################### DRAW USAGE ##########
#def drawSampleUsage(sample, length, outdir):
#    seqs = sample.len2seqs[length]
#    size2pos2aa2usage = aaSampleProfile(seqs, length, [1])
#    pos2aa2usage = size2pos2aa2usage[1]

    

def checkOptions(parser, args, options):
    if not options.indir:
        raise InputOptionError("Input directory is required. None was given.\n")
    if not os.path.exists(options.indir):
        raise InputOptionError("Input directory %s does not exist\n" %options.indir)
    if not os.path.exists(options.outdir):
        system("mkdir -p %s" %options.outdir)

def initOptions(parser):
    parser.add_option('-i', '--indir', dest='indir', help='Required argument. Input directory containing the fasta files')
    parser.add_option('-o', '--outdir', dest='outdir', default='.', help='Output directory. Default=%default')
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true', default=False, help='If specified, print detailed messages of what steps the program are at')

def main():
    usage = "Usage: %prog [options]\n" 
    parser = OptionParser(usage=usage)
    initOptions(parser)
    options, args = parser.parse_args()
    checkOptions(parser, args, options)

    starttime = time.clock()

    samples = readfiles( options.indir, 1, 1 )
    if options.verbose:
        sys.stderr.write('Read files in %d seconds.\n' %(time.clock() - starttime) )
    stime = time.clock()

    setLen2seqs(samples)
    if options.verbose:
        sys.stderr.write("setLen2seqs in %d seconds.\n" %(time.clock() - stime) )

    stime = time.clock()
    getAaUsage(samples, options.outdir, options.verbose)
    if options.verbose:
        sys.stderr.write('getAaUsage in %d seconds.\n' %(time.clock() -stime))
        sys.stderr.write('Total run time: %d seconds.\n' %(time.clock() - starttime) )

if __name__ == '__main__':
    main()

