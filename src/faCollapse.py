#!/usr/bin/env python

#nknguyen at soe ucsc edu
#Aug 11 2011
#Input: fasta file with the header formatted as:
#>sample;id|v1[,v2,...]|j1[,j2,...]
#Example: >SBC2;Cluster5|TRBV6-3,TRBV6-2|TRBJ2-3;size=365304
#
#Output: collapse identical sequences (regardless of V, J genes), and update the count
#Output header format: "sampleName;Clusters<CommaSeparatedListOfClusterHeaders>;size=###"
#Example: >SBC2;Clusters76255|TRBV6-1|TRBJ1-2,49098|TRBV6-5|TRBJ1-2;size=6

import os, sys, re
from optparse import OptionParser

def sortDictByValue( dictionary ):
    items = [ (v, k) for k,v in dictionary.items() ]
    items.sort()
    items.reverse()
    return [ (k, v) for v, k in items ]

def addSeq(seqDict, header, seq):
    items = header.split(';')
    #count = int( items[2].lstrip('size=') )
    count = int( items[ len(items) -1 ].lstrip('size=') )
    sample = items[0]
    cluster = items[1].lstrip('Cluster')
    if sample not in seqDict:
        seqDict[ sample ] = { seq: (count, [cluster]) }
    else:
        if seq not in seqDict[ sample ]:
            seqDict[ sample ][ seq ] = ( count, [cluster] )
        else:
            (currCount, currList) = seqDict[ sample ][ seq ]
            currList.append( cluster )
            seqDict[ sample ][ seq ] = ( currCount + count, currList )
    return

def collapse(infh, outfh):
    seq = ""
    header = ""
    seqDict = {} #key = the sequence. Value = (count, [(V, J), (V, J), ...])

    for line in infh.readlines():
        if not re.search('\w', line):
            continue
        line = line.strip()
        if line[0] == '>':
            if header != "" and seq != "":
                addSeq( seqDict, header, seq )
            header = line.lstrip('>')
            seq = ""
        else:
            seq += line

    if header != "" and seq != "":
        addSeq( seqDict, header, seq )

    #Print out collapsed sequences:
    for sample in seqDict:
        for seq in seqDict[ sample ]:
            (count, clusters) = seqDict[ sample ][ seq ]
            header = "%s;Clusters%s;size=%d" %( sample, ','.join(clusters), count )
            outfh.write( ">%s\n" %header)
            outfh.write( "%s\n" %seq )

def main():
    usage = "Usage: %prog infile outfile\n"
    parser = OptionParser(usage = usage)
    options, args = parser.parse_args()
    
    infh = open(args[0], 'r')
    outfh = open(args[1], 'w')
    collapse(infh, outfh)
    infh.close()
    outfh.close()

if __name__ == "__main__":
    main()
