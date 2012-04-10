#!/usr/bin/env python

#nknguyen at soe ucsc edu
#Aug 11 2011
#Input: fasta file with the header formatted as:
#>sample;id|v1[,v2,...]|j1[,j2,...]
#Example: >SBC2;Cluster5|TRBV6-3,TRBV6-2|TRBJ2-3;size=365304
#
#Output: 
#Mode 1:
#collapse identical sequences regardless of V, J genes, and update the count
#Output header format: "sampleName;Clusters<DoubleCommaSeparatedListOfClusterHeaders>;size=###"
#Example: >SBC2;Clusters76255|TRBV6-1|TRBJ1-2,,49098|TRBV6-5|TRBJ1-2;size=6
#
#Mode 2:
#Collapse identical sequences who must have at least one V in common and at least one J in common


import os, sys, re
from optparse import OptionParser

class Genes:
    def __init__(self):
        self.vs = []
        self.js = []

    def setGenes(self, vs, js):
        self.vs = vs
        self.js = js

    def updateGenes(vs, js):
        for v in vs:
            if v not in self.vs:
                self.vs.append(v)
        for j in js:
            if j not in self.js:
                self.js.append(j)

def sortDictByValue( dictionary ):
    items = [ (v, k) for k,v in dictionary.items() ]
    items.sort()
    items.reverse()
    return [ (k, v) for v, k in items ]

def getGenes(name):
    items = name.split('|')
    if len(items) < 3:
        raise ValueError("Wrong header format. Expected: id|v1<,v2,...>|j1<,j2,...>, see: %s\n" %name)
    vs = items[1]
    js = items[2]
    #vs = items[1].split(',')
    #js = items[2].split(',')
    return vs, js

def addSeq(seqDict, header, seq, sepGene):
    items = header.split(';')
    #count = int( items[2].lstrip('size=') )
    count = int( items[ len(items) -1 ].lstrip('size=') )
    sample = items[0]
    cluster = items[1].lstrip('Cluster')
    vs, js = getGenes(cluster)

    if sample not in seqDict:
        seqDict[ sample ] = { seq:{(vs, js) : (count, [cluster])} }
    else:
        if seq not in seqDict[ sample ]:
            seqDict[ sample ][ seq ] = {(vs, js): (count, [cluster])}
        else:
            vjs2count = seqDict[sample][seq]
            genes = checkGenes(vs, js, vjs2count) 
            if sepGene and not genes:
                vjs2count[ (vs, js) ] = (count, [cluster])
            else:
                if not sepGene:
                    genes = vjs2count.keys()[0]
                (currCount, currList) = vjs2count[genes]
                currList.append( cluster )
                newvs = mergeLists(vs, genes[0])
                newjs = mergeLists(js, genes[1])
                del seqDict[sample][seq][genes]
                seqDict[sample][seq][ (newvs, newjs) ] = (currCount + count, currList)
                
                #seqDict[ sample ][ seq ][ genes ] = ( currCount + count, currList )
    return

def checkGenes(vs, js, vjs2count):
    vsList = vs.split(',')
    jsList = js.split(',')

    for (currVs, currJs) in vjs2count:
        currVsList = currVs.split(',')
        currJsList = currJs.split(',')
        hasV = hasOverlap(vsList, currVsList)
        hasJ = hasOverlap(jsList, currJsList)
        if hasV and hasJ:#has at least one V and one J overlapped, merge the sequence to this crowd.
            return (currVs, currJs)
    return None

def hasOverlap(list1, list2):
    for i1 in list1:
        if i1 in list2:
            return True
    return False

def mergeLists(genes1, genes2):
    """Merge list1 to list2"""
    list1 = genes1.split(',')
    list2 = genes2.split(',')
    for i1 in list1:
        if i1 not in list2:
            list2.append(i1)
    return ",".join(sorted(list2))
    
def collapse(infh, outfh, sepGene):
    seq = ""
    header = ""
    seqDict = {} #key = sample, val = dictionary where: key = the sequence Value = {([vs], [js]): (count, [(V, J), (V, J), ...])}

    for line in infh.readlines():
        if not re.search('\w', line):
            continue
        line = line.strip()
        if line[0] == '>':
            if header != "" and seq != "":
                addSeq( seqDict, header, seq, sepGene )
            header = line.lstrip('>')
            seq = ""
        else:
            seq += line

    if header != "" and seq != "":
        addSeq( seqDict, header, seq, sepGene )

    #Print out collapsed sequences:
    for sample in seqDict:
        for seq in seqDict[ sample ]:
            for genes in seqDict[sample][seq]:
                (count, clusters) = seqDict[ sample ][ seq ][genes]
                header = "%s;Clusters%s;size=%d" %( sample, ',,'.join(clusters), count )
                outfh.write( ">%s\n" %header)
                outfh.write( "%s\n" %seq )

def main():
    usage = "Collapsing redundant (amino acid) sequences.\nUsage: %prog infile outfile\n"
    parser = OptionParser(usage = usage)
    parser.add_option("-g", "--genes", dest='gene', action = 'store_true', default=False, help='If specified, only collapse sequences with identical CDR3 and have at least 1 V and 1 J in common. Default=%defaul, collapse all identical sequences, regardless of which V and J genes they map to.')
    options, args = parser.parse_args()
    
    infh = open(args[0], 'r')
    outfh = open(args[1], 'w')

    collapse(infh, outfh, options.gene)
    infh.close()
    outfh.close()

if __name__ == "__main__":
    main()
