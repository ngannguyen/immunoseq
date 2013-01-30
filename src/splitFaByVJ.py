#!/usr/bin/env python

#nknguyen at soe ucsc edu
#Aug 24 2011
#Input: Directory contains input fasta files (one per sample) with the following format:
#>sampleName;id|V[,V2,V3...]|J[,J2,...]|D[,D2];size=###
#
#Output:
#Split the sequences per the V-J combination that they mapped to, and 
#concatenate sequences from all samples to that v.j.fa file
#

import os, sys, re
from optparse import OptionParser

def sortDictByValue( dictionary ):
    items = [ (v, k) for k,v in dictionary.items() ]
    items.sort()
    items.reverse()
    return [ (k, v) for v, k in items ]

def getSeqInfo( header ):
    items = header.split(";")
    if len(items) < 3:
        sys.stderr.write("Input file with wrong header format: %s\n. The format expected is:\n>sampleName;id|v[,v1,v2..]|j[,j1,j2,..];size=##\n" %header)
        sys.exit(1)
    size = int( items[-1].lstrip('size=') )
    list1 = items[1].split("|")
    if len(list1) < 3:
        sys.stderr.write("Input file with wrong header format: %s\n. The format expected is:\n>sampleName;id|v[,v1,v2..]|j[,j1,j2,..];size=##\n" %header)
        sys.exit(1)
    id = list1[0]
    vlist = list1[1].split(',')
    jlist = list1[2].split(',')
    
    return id, vlist, jlist, size

def writeSeq( outdir, header, seq, sampleName, v2seq ):
    id, vlist, jlist, size = getSeqInfo( header )
    for v in vlist:
        v = v.replace("/", "-")
        currseq = seq
        if v2seq and v in v2seq:
            currseq = v2seq[v] + seq

        for j in jlist:
            j = j.replace("/", "-")
            outfile = os.path.join(outdir, "%s.%s.fa" %(v, j))
            outfh = open( outfile, "a")
            if not os.path.exists( outfile ):
                outfh = open( outfile, "w")
            outfh.write( ">%s;%s;size=%d\n" %(sampleName, id, size) )
            outfh.write( "%s\n" % currseq )
            outfh.close()
    return

def splitSeqs(indir, outdir, v2seq):
    for infile in os.listdir( indir ):
        infh = open( os.path.join(indir, infile), "r" )
        sampleName = infile.split('.')[0]
        seq = ""
        header = ""
        for line in infh.readlines():
            line = line.strip()
            if len(line) <= 0:
                continue
            if line[0] == ">":
                if seq != "":
                    writeSeq( outdir, header, seq, sampleName, v2seq )
                header = line.lstrip('>')
                seq = ""
            else:
                seq += line
        if seq != "":
            writeSeq( outdir, header, seq, sampleName, v2seq )
        infh.close()
    return

def readVfile(file):
    v2seq = {}
    f = open(file, 'r')
    header = ''
    seq = ''
    for line in f:
        line = line.strip()
        if len(line) == 0:
            continue
        if line[0] == '>':
            if header != '' and seq != '' and header not in v2seq:
                v2seq[header] = seq
            header = line.lstrip('>')
            seq = ''
        else:
            seq += line
    if header != '' and seq != '' and header not in v2seq:
        v2seq[header] = seq
    f.close()
    return v2seq

def main():
    usage = "Usage: %prog inputDirectory outputDirectory\n"
    parser = OptionParser(usage = usage)
    parser.add_option('-v', '--addV', dest='vfile', help='If specified, add the rest of V gene to each sequence in the output fasta files. Default = None')
    options, args = parser.parse_args()

    v2seq = {}
    if options.vfile:
        v2seq = readVfile(options.vfile)

    splitSeqs( args[0], args[1], v2seq )

if __name__ == "__main__":
    main()
