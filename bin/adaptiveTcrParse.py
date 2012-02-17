#!/usr/bin/env python

#nknguyen at soe ucsc edu
#Jan 10 2011
#Parse adaptive TCR's output tsv file
#
#Input: Directory that contains Adaptive Biotechnologies tsv files:
#sequenceID container nucleotide aminoAcid normalizedFrequency normalizedCopy rawFrequency copy cdr3Length VFamilyName VGeneName VTies DGeneName JGeneName JTies VDeletion d5Deletion d3Deletion JDeletion n2Insertion n1Insertion sequenceStatus VIndex n1Index n2Index DIndex JIndex
#
#Output: Output fasta files of productive sequences to output directory
#   The header of each sequence has the following format: 
#   >sampleName;userDefinedCloneID|Vgenes|Jgenes;size=#

import os, re, sys, copy
from scipy.stats.stats import pearsonr, spearmanr, kendalltau

class Sequence:
    def __init__(self, line):
        items = line.strip().split('\t')
        if len(items) < 27:
            sys.stderr.write('Wrong tsv format. Expected 27 fields, only have %d\n%s\n' %(len(items), line))
            sys.exit(1)
        self.id = items[0]
        self.nuc = items[2]
        self.aa = items[3]
        self.normFreq = -1.0
        self.normCount = -1
        if items[4] != '':
            self.normFreq = float(items[4])
        if items[5] != '':
            self.normCount = int(items[5])
        self.freq = float(items[6])
        self.count = int(items[7])
        self.cdr3len = int(items[8])
        self.vfam = items[9]
        self.vgene = items[10]
        self.vties = items[11]
        self.dgene = items[12]
        self.jgene = items[13]
        self.jties = items[14]
        self.vdel = int(items[15])
        self.d5del = int(items[16])
        self.d3del = int(items[17])
        self.jdel = int(items[18])
        self.n2ins = int(items[19])
        self.n1ins = int(items[20])
        self.status = items[21]
        self.vindex = int(items[22])
        self.n1index = int(items[23])
        self.n2index = int(items[24])
        self.dindex = int(items[25])
        self.jindex = int(items[26])

def readFile(file):
    f = sys.stdin
    hasNorm = True
    if file != '-' and os.path.exists(file):
        f = open(file, 'r')
    seqs = []
    for line in f:
        items = line.strip().split('\t')
        if len(items) < 27 or items[0] =='sequenceID':
            continue
        seq = Sequence(line)
        if seq.normFreq < 0 or seq.normCount < 0:
            hasNorm = False
        if seq.status == 'Productive':
            seqs.append(seq)
    f.close()
    return seqs, hasNorm

def printSeqs(seqs, outfile, useNorm):
    samplename = os.path.basename(outfile).rstrip(".fa")
    f = sys.stdout
    if outfile != '-':
        f = open(outfile, 'w')
    
    i = 0
    for seq in seqs:
        if seq.status != 'Productive':
            continue

        vstr = seq.vgene
        if seq.vties != '':
            vstr = ','.join( seq.vties.split(', ') )
        jstr = seq.jgene
        if seq.jties != '':
            jstr = ','.join( seq.jties.split(', ') )
        
        header = ">%s;%s|%s|%s;size=%d" %(samplename, str(i), vstr, jstr, seq.normCount)
        if not useNorm:
            header = ">%s;%s|%s|%s;size=%d" %(samplename, str(i), vstr, jstr, seq.count)
        f.write("%s\n" %header)
        f.write("%s\n" %seq.aa)
        i+=1

def main():
    indir = sys.argv[1]
    outdir = sys.argv[2]
    useNorm = True

    for file in os.listdir(indir):
        sample = file.split('.')[0]
        if len(file.split('.')) < 2 or file.split('.')[1] != 'tsv':
            continue
        seqs, hasNorm = readFile( os.path.join(indir,file) )
        if not hasNorm:
            useNorm = False
            sys.stderr.write("Some sample/sequence does not have normalized frequency or normalized count. Using raw count instead\n")
        outfile = os.path.join(outdir,"%s.fa" % sample)
        printSeqs(seqs, outfile, useNorm)
        #printSeqs(seqs, outfile, False)


if __name__ == '__main__':
    main()














        

