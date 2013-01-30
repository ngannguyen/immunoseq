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
#Example: >SBC2;Clusters76255|TRBV6-1|TRBJ1-2|TRBD1-1,,49098|TRBV6-5|TRBJ1-2|TRBD1-1;size=6
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

def getGeneFamilies(geneString, sep):
    items = geneString.split(',')
    fams = []
    for gene in items:
        fam = gene.split(sep)[0]
        if fam not in fams:
            fams.append(fam)
    return ','.join(fams)

def getGenes(name, family, dgene):
    items = name.split('|')
    if len(items) < 3:
        raise ValueError("Wrong header format. Expected: id|v1<,v2,...>|j1<,j2,...>, see: %s\n" %name)
    vs = items[1]
    js = items[2]
    ds = ''
    if dgene and len(items) == 4:
        ds = items[3]
    
    if family:
        vs = getGeneFamilies(vs, '-')
        js = getGeneFamilies(js, '-')
    else:
        vs = getGeneFamilies(vs, '*')
        js = getGeneFamilies(js, '*')
    return vs, js, ds

def removeSubtype(cluster):
    items = cluster.split('|')
    if len(items) < 3:
        raise ValueError("Wrong header format. Expected: id|v1<,v2,...>|j1<,j2,...>, see: %s\n" %name)
    vs = getGeneFamilies( items[1], "*" )
    js = getGeneFamilies( items[2], "*" )
    items[1] = vs
    items[2] = js
    return '|'.join(items)

def addSeq(seqDict, header, seq, sepGene, family, dgene):
    items = header.split(';')
    count = int( items[ len(items) -1 ].lstrip('size=') )
    sample = items[0]
    cluster = items[1].lstrip('Cluster')
    cluster = removeSubtype(cluster) #remove subtype of genes
    vs, js, ds = getGenes(cluster, family, dgene)

    if sample not in seqDict:
        seqDict[ sample ] = { seq:{(vs, js, ds) : (count, [cluster])} }
    else:
        if seq not in seqDict[ sample ]:
            seqDict[ sample ][ seq ] = {(vs, js, ds): (count, [cluster])}
        else:
            vjs2count = seqDict[sample][seq]
            genes = checkGenes(vs, js, ds, vjs2count) 
            if sepGene and not genes:
                vjs2count[ (vs, js, ds) ] = (count, [cluster])
            else:
                if not sepGene:
                    genes = vjs2count.keys()[0]
                (currCount, currList) = vjs2count[genes]
                currList.append( cluster )
                newvs = mergeLists(vs, genes[0])
                newjs = mergeLists(js, genes[1])
                newds = mergeLists(ds, genes[2])
                del seqDict[sample][seq][genes]
                seqDict[sample][seq][ (newvs, newjs, newds) ] = (currCount + count, currList)
                
                #seqDict[ sample ][ seq ][ genes ] = ( currCount + count, currList )
    return

#def postCollapse(seqDict):
#Due to greedy merging/collapsing algorithm, sometimes miss-merge clones with common genes. 
#For example: clone1 with V1,V2 and clone2 with V3 and clone 3 with V2,V3. 
#Start at clone1. NExt clone 2, no common V with 1, don't merge. Then get to 3: ah hah, similar with 1, so merge 3 to 1. And 2 left alone.
#postCollapse to take care of these cases...
#    for sample in seqDict:
#        seq2vj2count = seqDict[sample]
#        for seq, vj2count in seq2vj2count.iteritems():#each aa sequence
#            vjSets = vj2count.keys()
#            if len(vjSets) == 1:
#                continue
#
#            #newVjSets = [ vjSets[0] ]
#            newVj2count = { vjSets[0]: vj2count[ vjSets[0] ] }
#
#            for i in xrange(1, len(vjSets)):
#                genes = vjSets[i]
#                vs = genes[0]
#                js = genes[1]
#                vsList = vs.split(',')
#                jsList = js.split(',')
#                flag = False #False = cant merge with anybody in newVjSets, True otherwise
#
#                for (prevVs, prevJs) in newVjSets:
#                    prevVsList = prevVs.split(',')
#                    prevJsList = prevJs.split(',')
#                    hasV = hasOverlap(vsList, prevVsList)
#                    hasJ = hasOverlap(jsList, prevJsList)
#                    if hasV and hasJ: #overlap with (prevVs, prevJs)
#                        (currCount, currList) = newVj2count[ (prevVs, prevJs) ]
#                        currList.append( cluster )
#                        
#                        newvs = mergeLists(prevVs, vs)
#                        newjs = mergeLists(prevJs, js)
#                        del 


def addSeq2(seqDict, header, seq, sepGene, family, dgene):
    #seqDict --> seq --> VJgenes --> (count, [clusters], {sample:currsamplecount}) 
    items = header.split(';')
    count = int( items[ len(items) -1 ].lstrip('size=') )
    sample = items[0]
    cluster = items[1].lstrip('Cluster')
    cluster = removeSubtype(cluster) #remove subtype of genes
    vs, js, ds = getGenes(cluster, family, dgene)

    if seq not in seqDict:
        seqDict[seq] = { (vs, js, ds): (count, [cluster], {sample:count}) }
    else:
        vjs2count = seqDict[seq]
        genes = checkGenes(vs, js, ds, vjs2count) 
        if sepGene and not genes:
            vjs2count[ (vs, js, ds) ] = (count, [cluster], {sample:count})
        else:
            if not sepGene:
                genes = vjs2count.keys()[0]
            (currCount, currList, currSample2count) = vjs2count[genes]
            currList.append( cluster )
            if sample not in currSample2count:
                currSample2count[sample] = count
            else:
                currSample2count[sample] += count

            newvs = mergeLists(vs, genes[0])
            newjs = mergeLists(js, genes[1])
            newds = mergeLists(ds, genes[2])
            del seqDict[seq][genes]
            seqDict[seq][ (newvs, newjs, newds) ] = (currCount + count, currList, currSample2count)
            
    return

def checkGenes(vs, js, ds, vjs2count):
    vsList = vs.split(',')
    jsList = js.split(',')
    dsList = ds.split(',')

    for (currVs, currJs, currDs) in vjs2count:
        currVsList = currVs.split(',')
        currJsList = currJs.split(',')
        currDsList = currDs.split(',')
        hasV = hasOverlap(vsList, currVsList)
        hasJ = hasOverlap(jsList, currJsList)
        hasD = hasOverlap(dsList, currDsList)
        if hasV and hasJ and hasD:#has at least one V and one J overlapped, merge the sequence to this crowd.
            return (currVs, currJs, currDs)
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
    
def collapse(infh, outfh, sepGene, family, name, dgene):
    seq = ""
    header = ""
    seqDict = {} #key = sample, val = dictionary where: key = the sequence Value = {([vs], [js], [ds]): (count, [(V, J, D), (V, J, D), ...])}

    for line in infh.readlines():
        if not re.search('\w', line):
            continue
        line = line.strip()
        if line[0] == '>':
            if header != "" and seq != "":
                if not name:
                    addSeq( seqDict, header, seq, sepGene, family, dgene )
                else:
                    addSeq2( seqDict, header, seq, sepGene, family, dgene )
            header = line.lstrip('>')
            seq = ""
        else:
            seq += line

    if header != "" and seq != "":
        if not name:
            addSeq( seqDict, header, seq, sepGene, family, dgene )
        else:
            addSeq2( seqDict, header, seq, sepGene, family, dgene )

    #Print out collapsed sequences:
    if not name:
        for sample in seqDict:
            for seq in seqDict[ sample ]:
                for genes in seqDict[sample][seq]:
                    (count, clusters) = seqDict[ sample ][ seq ][genes]
                    header = "%s;%s;size=%d" %( sample, ',,'.join(clusters), count )
                    outfh.write( ">%s\n" %header)
                    outfh.write( "%s\n" %seq )
    else:
        #seqDict --> seq --> VJgenes --> (count, [clusters], {sample:count}) 
        for seq in seqDict:
            for genes in seqDict[seq]:
                (aggcount, clusters, sample2count) = seqDict[seq][genes]
                count = max( sample2count.values() ) #June 29 2012, if merge samples, count = max(samplecounts)
                samples = sorted(sample2count.keys())
                countStr = ','.join( [str(sample2count[s]) for s in samples] )
                header = "%s;%s;%s;size=%d" %( ','.join(samples), ',,'.join(clusters), countStr, count )
                #header = "%s;%s;size=%d" %( ','.join(samples), ',,'.join(clusters),  count )
                outfh.write( ">%s\n" %header)
                outfh.write( "%s\n" %seq )

def main():
    usage = "Collapsing redundant (amino acid) sequences.\nUsage: %prog infile outfile\n"
    parser = OptionParser(usage = usage)
    parser.add_option("-g", "--genes", dest='gene', action = 'store_true', default=False, help='If specified, only collapse sequences with identical CDR3 and have at least 1 V and 1 J in common. Default=%defaul, collapse all identical sequences, regardless of which V and J genes they map to.')
    parser.add_option("-n", "--sampleName", dest='name', action='store_true', default=False, help='If specified, collapse same sequences of different samples')
    parser.add_option("-f", "--family", dest='fam', action = 'store_true', default=False, help='If specified, and if option "genes" is true, collapse sequences with identical CDR3 sequence and have at least 1 V family (not gene) and 1J family (not gene) in common. Default=%default')
    parser.add_option("-d", dest="d", action="store_true", default=False, help='Only valid when the option --genes is specified. If specified, the program will merge only clones with identical sequence, at least 1 common V, at least 1 common J, AND at least 1 common D (not just common V and J)')

    options, args = parser.parse_args()
    
    infh = open(args[0], 'r')
    outfh = open(args[1], 'w')

    collapse(infh, outfh, options.gene, options.fam, options.name, options.d)
    infh.close()
    outfh.close()

if __name__ == "__main__":
    main()
