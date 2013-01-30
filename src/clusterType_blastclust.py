#!/usr/bin/env python2.6

"""
Tue Nov 27 14:50:38 PST 2012
Parse blastclust output files and call the different cluster types.
This is applied in searching for similar (but not identical) clones to the expanded clones.

The input blast was done on the top clones of each patients:
For each V-J recombination:
    Blast the top clones with that V-J against all the clones carrying that V-J
    
In a way, we use each top clone as a seed and look for similar sequences using that seed. The hits from the blast that pass a certain similarity cutoff form the cluster.
We distinguish 5 types of cluster:
    #Types: 1: "Multiple expanded clones from multiple samples"
    #       2: "Multiple expanded clones from only one sample"
    #       3: "Expanded clones that do not cluster with any other sequence"
    #       4: "Non-expanded clones that do not cluster with any other sequence"
    #       5: "Others"

Return a summary of:

nknguyen at soe ucsc edu
"""

import os, sys, re
from Bio.Blast import NCBIXML
import immunoseq.lib.immunoseqLib as iseqlib

def getCloneInfo(clonestr):
    #as11D;183042;size=8925
    items = clonestr.lstrip('>').split(';')
    sample = items[0]
    size = int(items[-1].lstrip("size="))
    id = items[1]
    return sample, id, size

class Clone():
    def __init__(self, clonestr):
        sample, id, size = getCloneInfo(clonestr)
        self.desc = clonestr.lstrip('>')
        self.sample = sample
        self.id = id
        self.size = size

    def setFreq(self, total):
        if total == 0:
            raise ValueError("Error: Total sequences of sample %s is 0." %(self.sample))
        else:
            self.freq = 100.0*self.size/total

class Cluster():
    def __init__(self, line, vj, sample2total):
        self.vj = vj
        items = line.strip().split()
        clones = [Clone(item) for item in items]
        #set clone frequency:
        for clone in clones:
            if sample2total:
                if clone.sample in sample2total:
                    clone.setFreq(sample2total[clone.sample])
                else:
                    raise ValueError("The total number of sequences is required but not found for sample %s\n" %clone.sample)

        self.clones = sorted(clones, key=lambda c: c.size, reverse=True)
        self.totalReads = sum([c.size for c in clones])
        self.numClones = len(self.clones)

    def getDesc(self):
        return "%s;%d;%d" %(self.vj, self.numClones, self.totalReads)

def readBlastclustOutfiles(infile, vjname, options): 
    #minPos, minLen, minExpSize, minExpClones, minMotifClones, minExpFreq, sample2total
    #Read in blastclust output file:
    f = open(infile, 'r')
    clusters = []

    for line in f:
        cluster = Cluster(line, vjname, options.sample2total)
        clusters.append(cluster)
    f.close()
    return clusters

def getClusterType(cluster, minExpFreq, minExpSize, minExpClones):
    #Types: 1: "Multiple expanded clones from multiple samples"
    #       2: "Multiple expanded clones from only one sample"
    #       3: "Expanded clones that do not cluster with any other sequence"
    #       4: "Non-expanded clones that do not cluster with any other sequence"
    #       5: "Others"
    numExp = 0
    expSamples = [] #samples with expanded clones besides seed sample
    clones = {} #non-redundant clone of the current cluster, key=cloneStr, value = clone sequence
    type = 5
    if cluster.numClones == 1: #type 3 or 4
        clone = cluster.clones[0]
        if clone.size >= minExpSize and clone.freq >= minExpFreq: #expanded clone, type 3
            numExp += 1
            type = 3
        else: #non-expanded clone, type 4
            type = 4
    else:
        for clone in cluster.clones:
            #expanded
            if clone.size >= minExpSize and clone.freq >= minExpFreq: #expanded clone
                numExp += 1
                if clone.sample not in expSamples:
                    expSamples.append(clone.sample)

        if numExp >= minExpClones: #type 1 or 2
            if len(expSamples) == 1: #no other sample has expanded clones
                type = 2
            else:
                type = 1
    return type, numExp

def getFh(type, fh1, fh2, fh3, fh4, fh5):
    if type == 1:
        return fh1
    elif type == 2:
        return fh2
    elif type == 3:
        return fh3
    elif type == 4:
        return fh4
    else:
        return fh5

def printCluster(cluster, fh):
    fh.write("\n>%s\n" %cluster.getDesc())
    fh.write("\t%s\n" %("\n\t".join([clone.desc for clone in cluster.clones])) )

def clusterId2info():
    id2info = { 1: "Multiple expanded clones from multiple samples",
                2: "Multiple expanded clones from only one sample",
                3: "Expanded clones that do not cluster with any other sequence",
                4: "Non-expanded clones that do not cluster with any other sequence",
                5: "Others",
              }
    return id2info

def separateClusters(clusters, minExpFreq, minExpSize, minExpClones, fh1, fh2, fh3, fh4, fh5, type2count, type2clusters):
    #Types: 1: "Multiple expanded clones from multiple samples"
    #       2: "Multiple expanded clones from only one sample"
    #       3: "Expanded clones that do not cluster with any other sequence"
    #       4: "Non-expanded clones that do not cluster with any other sequence"
    #       5: "Others"
    totalseq = 0
    expseq = 0
    for cluster in clusters:
        type, numExp = getClusterType(cluster, minExpFreq, minExpSize, minExpClones)
        expseq += numExp
        totalseq += cluster.numClones
        fh = getFh(type, fh1, fh2, fh3, fh4, fh5)
        printCluster(cluster, fh) 
        
        #Update type count
        if type not in type2count:
            type2count[type] = 1
        else:
            type2count[type] += 1
        #Update type2clusters:
        if type not in type2clusters:
            type2clusters[type] = [cluster]
        else:
            type2clusters[type].append(cluster)
    return totalseq, expseq

def getPc(count, total):
    if total == 0:
        return 0
    return 100.0*count/total

def readFaFile(file):
    header2seq = {}
    f = open(file, 'r')
    header = ''
    seq = ''
    for line in f:
        line = line.strip()
        if len(line) == 0:
            continue
        if line[0] == '>':
            if header != '' and header not in header2seq:
                header2seq[header] = seq
            header = line.lstrip('>')
            seq = ''
        else:
            seq += line
    if header != '' and seq != '' and header not in header2seq:
        header2seq[header] = seq
    f.close()
    return header2seq

def readSample2total(file):
    sample2total = {}
    f = open(file, 'r')
    for line in f:
        items = line.strip().split()
        sample2total[items[0]] = int(items[1])
    f.close()
    return sample2total

def addOptions(parser):
    parser.add_option('-i', '--indir', dest='indir', help='Input directory, containing blastclust output files')
    parser.add_option('-f', '--fasta', dest='fadir', help='Fasta directory, containing fasta files')
    parser.add_option('-o', '--outdir', dest='outdir', help='Output directory')
    parser.add_option('-b', '--basename', dest='basename', default='out', help='Output files basename. Default=%default')
    parser.add_option('-S', '--minExpandedSize', dest='minExpSize', type='int', default=1000, help='Minimum number of reads for a clone to be called "expanded". Default=%default')
    parser.add_option('-C', '--minExpandedClones', dest='minExpClones', type='int', default=2, help='Minimum number of expanded clones for a cluster to be classified as type 1. Default=%default')
    parser.add_option('-F', '--minExpandedFreq', dest='minExpFreq', type='float', default=0.0, help='Minimun frequency (ranging from 0 - 100%) for a clone to be called "expanded". Default=%default')
    parser.add_option('--sample2total', dest='sample2total', help='Required if --minExpandedFreq is larger than 0. Format: <sample> <totalCount>')
    parser.add_option('-v', '--removeV', dest='vfile', help='If specified, remove the rest of V gene and only print out the CDR3 region of each sequence in the output fasta files. Default = None')
    #parser.add_option('-p', '--positive', dest='minPos', type='float', default=0.9, help='Minimum portion of positive matches. Default=%default')
    #parser.add_option('-l', '--len', dest='minLen', type='int', default=10, help='Minimum sequence length to be included in the output. Default=%default')
    #parser.add_option('-c', '--minMotifClones', dest='minMotifClones', type='int', default=10, help='Minimum number of clones carrying the same motif for a seed and its cluster to be classified as type 2. Default=%default ')

def main():
    parser = iseqlib.initOptions()
    addOptions(parser)
    options, args = parser.parse_args()

    if options.minExpFreq > 0 and not options.sample2total:
        parser.error("--sample2total is required as --minExpandedFreq > 0\n")
    if options.sample2total:
        options.sample2total = readSample2total(options.sample2total)

    ext = 'txt'
    infiles = iseqlib.getfiles(options.indir, ext)
    outbasename = os.path.join(options.outdir, options.basename)
    
    fh1 = open("%s-1.txt" %outbasename, 'w')
    fh2 = open("%s-2.txt" %outbasename, 'w')
    fh3 = open("%s-3.txt" %outbasename, 'w')
    fh4 = open("%s-4.txt" %outbasename, 'w')
    fh5 = open("%s-5.txt" %outbasename, 'w')
    
    type2count = {}
    type2clusters = {}
    totalseq = 0
    expseq = 0
    for file in infiles:
        infile = os.path.join(options.indir, file)
        vjname = file.rstrip(ext).rstrip('.')
        clusters = readBlastclustOutfiles(infile, vjname, options)
        total, exp = separateClusters(clusters, options.minExpFreq, options.minExpSize, options.minExpClones, fh1, fh2, fh3, fh4, fh5, type2count, type2clusters)
        totalseq += total
        expseq += exp
    fh1.close()
    fh2.close()
    fh3.close()
    fh4.close()
    fh5.close()

    #summary stats:
    summaryFile = "%s-summary.txt" %outbasename
    fh = open(summaryFile, 'w')
    fh.write("#min size to be called expanded clones: %d\n#min number of expanded clones required for type 1: %d\n" %(options.minExpSize, options.minExpClones))
    id2info = clusterId2info()
    fh.write("Total number of clones: %d\n" %totalseq)
    fh.write("Total number of expanded clones: %d\n" %expseq)
    for type in sorted(type2count.keys()):
        fh.write("%d\t%s\t%d\n" %(type, id2info[type], type2count[type]))
    fh.close()

    #Print type2clones:
    
    #Read v files:
    v2seq = {}
    if options.vfile:
        v2seq = readFaFile(options.vfile)
    
    #Read fasta files:
    ext = 'fa'
    infiles = iseqlib.getfiles(options.fadir, ext)
    vj2id2seq = {}
    for file in infiles:
        infile = os.path.join(options.fadir, file)
        vjname = file.rstrip(ext).rstrip('.')
        vj2id2seq[vjname] = readFaFile(infile)

    for type, clusters in type2clusters.iteritems():
        outfile = "%s-%d.fa" %(outbasename, type)
        f = open(outfile, 'w')
        f.write("#Type: %s\n" %type)
        clusters = sorted(clusters, key=lambda c: c.totalReads, reverse=True)
        for cluster in clusters:
            f.write("#Cluster %s\n" %cluster.getDesc())
            id2seq = vj2id2seq[cluster.vj]
            v = cluster.vj.split(".")[0]
            vseq = ""
            if v in v2seq:
                vseq = v2seq[v]

            for clone in cluster.clones:
                seq = id2seq[ clone.desc ]
                cdr3seq = seq[ len(vseq): ]
                header = "%s;%s;%d;%s" %(clone.sample, cdr3seq, clone.size, cluster.vj)
                f.write(">%s\n" %header)
                f.write("%s\n" %seq)
        f.close()

if __name__ == '__main__':
    main()


