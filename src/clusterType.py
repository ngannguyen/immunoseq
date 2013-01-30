#!/usr/bin/env python2.6

"""
Tue Nov 27 14:50:38 PST 2012
Parse Blast XML output files and call the different cluster types.
This is applied in searching for similar (but not identical) clones to the expanded clones.

The input blast was done on the top clones of each patients:
For each V-J recombination:
    Blast the top clones with that V-J against all the clones carrying that V-J
    
In a way, we use each top clone as a seed and look for similar sequences using that seed. The hits from the blast that pass a certain similarity cutoff form the cluster.
We distinguish 4 types of cluster:
    1/ Seed has multiple expanded clones as matches.
        1.1: The expanded clones are from the same sample with seed
        1.2: The expanded clones are from at least one sample different from seed sample
    2/ Seed has many small clones with similar motifs
        Cutoffs includes: a/ minimum number of clones contribute to one motif, b/ number of samples
        2.1: the clones carrying the motif are from the same sample with seed
        2.2: the clones carrying the motif are from at least one different sample than seed sample
    3/ Everything else (similar to type 2 but did not pass the cutoffs, i.e seed has a small number of low-frequency hits, or too many motifs but not enough clones to support a single motif)
    4/ Seed has no similar clones

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

def getMatchCloneInfo(matchStr):
    #gnl|BL_ORD_ID|3115 as1D;46376;size=2
    cloneStr = matchStr.split()[-1]
    return getCloneInfo(cloneStr)

def getHitSize(cloneStr):
    sample, id, size = getCloneInfo(cloneStr)
    return size

def isSameClone(clone, hitclone, hit):
    sample, id, size = getCloneInfo(clone)
    hitsample, hitid, hitsize = getCloneInfo(hitclone)
    if sample == hitsample and size == hitsize and hit.query == hit.sbjct:
        return True
    else:
        return False

def readNcbiXml(infile, options): 
    #minPos, minLen, minExpSize, minExpClones, minMotifClones, minExpFreq, sample2total
    #Read in blast-output file:
    rh = open(infile)
    records = NCBIXML.parse( rh)

    clone2hits = {} #key = cloneName, val = [ (list of hit papers, identity) ]

    for record in records: #each seed
        if record.query_length < options.minLen: #ignore seeds that are shorter than minimum required length
            continue
        clone = record.query
        clone2hits[clone] = []
        for aln in record.alignments:
            for hit in aln.hsps: #each hit
                #if float(hit.identities)/len(hit.query) < minPositives:
                if float(hit.positives)/len(hit.query) < options.minPos: #ignore matches with similarity lower than required
                    continue
                if aln.title.split()[-1] == clone or isSameClone(clone, aln.title.split()[-1], hit): #hit is the seed itself, ignore
                    continue

                #if len(hit.match) < options.minLen:
                if len(hit.match) < record.query_length:
                    continue
                clone2hits[ clone ].append( (aln.title.split()[-1], hit.positives, hit.query, hit.match, hit.sbjct) )
            #Sort the hits by size
            hits = clone2hits[clone]
            clone2hits[clone] = sorted(hits, key=lambda h:getHitSize(h[0]), reverse=True)
    return clone2hits

def getClusterType(seedSample, hits, sample2total, minExpFreq, minExpSize, minExpClones, minMotifClones):
    numExp = 0
    expSamples = [] #samples with expanded clones besides seed sample
    motif2count = {}
    motif2samples = {}
    clones = {} #non-redundant clone of the current cluster, key=cloneStr, value = clone sequence
    
    if len(hits) == 0:
        return "No match", motif2count, clones
        #return "No match", numExp, expSamples, motif2count, motif2samples #type 3, noHit

    for hit in hits:
        sample, id, size = getCloneInfo(hit[0])
        #expanded
        total = sample2total[sample]
        if size >= minExpSize and float(size)/total >= minExpFreq: #expanded clone
            numExp += 1
            if sample not in expSamples and sample != seedSample:
                expSamples.append(sample)
        #motif 
        match = hit[3]
        motif = match.replace(" ", ".").replace("-", ".")
        if motif not in motif2count:
            motif2count[motif] = 1
            if sample != seedSample:
                motif2samples[motif] = [sample]
            else:
                motif2samples[motif] = []
        else:
            motif2count[motif] += 1
            if sample != seedSample and sample not in motif2samples[motif]:
                motif2samples[motif].append(sample)
        
        #update clones
        if hit[0] not in clones:
            clones[ hit[0] ] = hit[4].replace("-", "")


    if numExp >= minExpClones: #type 1
        if len(expSamples) == 0: #no other sample has expanded clones
            type = "Matches of expanded clones within the seed sample"
        else:
            type = "Matches of expanded clones in multiple samples"
    else:
        type = "Matches of small clones, different motifs"
        for motif, count in motif2count.iteritems():
            if count >= minMotifClones:
                type = "Multiple small clones with same motif"
                break
    return type, motif2count, clones
    #return type, numExp, expSamples, motif2count, motif2samples

def getFh(type, fh1, fh2, fh3, fh4):
    if type == "Matches of expanded clones within the seed sample" or type == "Matches of expanded clones in multiple samples":
        return fh1
    elif type == "Multiple small clones with same motif":
        return fh2
    elif type == "Matches of small clones, different motifs":
        return fh3
    else:
        return fh4

def getFaFilename(type, outbasename):
    if type == "Matches of expanded clones within the seed sample":
        return "%s-1.2.fa" %(outbasename)
    elif type == "Matches of expanded clones in multiple samples":
        return "%s-1.1.fa" %(outbasename)
    elif type == "Multiple small clones with same motif":
        return "%s-2.fa" %(outbasename)
    elif type == "Matches of small clones, different motifs":
        return "%s-3.fa" %(outbasename)
    else:
        return "%s-4.fa" %(outbasename)

def printCluster(clone, vjnames, hits, motif2count, minMotifClones, fh):
    motifs = []
    for m, c in motif2count.iteritems():
        if c >= minMotifClones:
            motifs.append("%s-%d" %(m, c))

    fh.write("\n>%s;%s;%s\n" %(clone, vjnames, ",".join(motifs)))
    for i, hit in enumerate(hits):
        fh.write("\t%s\n" %hit[0])
        fh.write("\t\t%s\n" %hit[2])
        fh.write("\t\t%s\n" %hit[3])
        fh.write("\t\t%s\n" %hit[4])

def separateClusters(clone2hits, sample2total, minExpFreq, minExpSize, minExpClones, minMotifClones, vjnames, fh1, fh2, fh3, fh4, type2count, type2clones):
    #Types: 1: "Matches of expanded clones within the seed sample" or "Matches of expanded clones in multiple samples"
    #       2: "Multiple small clones with same motif"
    #       3: "Matches of small clones, different motifs"
    #       4: "No match"

    for clone, hits in clone2hits.iteritems():
        sample, id, size = getCloneInfo(clone)
        type, motif2count, clones = getClusterType(sample, hits, sample2total, minExpFreq, minExpSize, minExpClones, minMotifClones)
        fh = getFh(type, fh1, fh2, fh3, fh4)
        printCluster(clone, vjnames, hits, motif2count, minMotifClones, fh) 
        
        #Update type count
        if type not in type2count:
            type2count[type] = 1
        else:
            type2count[type] += 1

        #Update type2clones
        if len(hits) > 0:
            seedSeq = hits[0][2].replace("-", "")
            clones[clone] = seedSeq
            if type not in type2clones:
                type2clones[type] = {}
                for header, seq in clones.iteritems():
                    type2clones[type][header] = {vjnames: seq}
            else:
                header2vj2seq = type2clones[type]
                for header, seq in clones.iteritems():
                    if header in header2vj2seq:
                        if vjnames not in header2vj2seq[header]:
                            header2vj2seq[header][vjnames] = seq
                    else:
                        header2vj2seq[header] = {vjnames: seq}
            
    return

def getPc(count, total):
    if total == 0:
        return 0
    return 100.0*count/total

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
            if header != '' and header not in v2seq:
                v2seq[header] = seq
            header = line.lstrip('>')
            seq = ''
        else:
            seq += line
    if header != '' and seq != '' and header not in v2seq:
        v2seq[header] = seq
    f.close()
    return v2seq

def readSample2total(file):
    sample2total = {}
    f = open(file, 'r')
    for line in f:
        items = line.strip().split()
        sample2total[items[0]] = int(items[1])
    f.close()
    return sample2total

def addOptions(parser):
    parser.add_option('-i', '--indir', dest='indir', help='Input directory')
    parser.add_option('-o', '--outdir', dest='outdir', help='Output directory')
    parser.add_option('-b', '--basename', dest='basename', default='out', help='Output files basename. Default=%default')
    parser.add_option('-p', '--positive', dest='minPos', type='float', default=0.9, help='Minimum portion of positive matches. Default=%default')
    parser.add_option('-l', '--len', dest='minLen', type='int', default=10, help='Minimum sequence length to be included in the output. Default=%default')
    parser.add_option('-S', '--minExpandedSize', dest='minExpSize', type='int', default=1000, help='Minimum number of reads for a clone to be called "expanded". Default=%default')
    parser.add_option('-F', '--minExpandedFreq', dest='minExpFreq', type='float', default=0.0, help='Minimun frequency for a clone to be called "expanded". Default=%default')
    parser.add_option('-C', '--minExpandedClones', dest='minExpClones', type='int', default=1, help='Minimum number of similar expanded clones for a seed and its cluster to be classified as type 1. Default=%default')
    parser.add_option('-c', '--minMotifClones', dest='minMotifClones', type='int', default=10, help='Minimum number of clones carrying the same motif for a seed and its cluster to be classified as type 2. Default=%default ')
    parser.add_option('--sample2total', dest='sample2total', help='Required if --minExpandedFreq is larger than 0. Format: <sample> <totalCount>')
    parser.add_option('-v', '--addV', dest='vfile', help='If specified, add the rest of V gene to each sequence in the output fasta files. Default = None')

def main():
    parser = iseqlib.initOptions()
    addOptions(parser)
    options, args = parser.parse_args()

    if options.minExpFreq > 0 and not options.sample2total:
        parser.error("--sample2total is required as --minExpandedFreq > 0\n")
    if options.sample2total:
        options.sample2total = readSample2total(options.sample2total)

    ext = 'xml'
    infiles = iseqlib.getfiles(options.indir, ext)
    outbasename = os.path.join(options.outdir, options.basename)
    
    fh1 = open("%s-1.txt" %outbasename, 'w')
    fh2 = open("%s-2.txt" %outbasename, 'w')
    fh3 = open("%s-3.txt" %outbasename, 'w')
    fh4 = open("%s-4.txt" %outbasename, 'w')
    
    type2count = {}
    type2clones = {}
    for file in infiles:
        infile = os.path.join(options.indir, file)
        vjname = file.rstrip(ext).rstrip('.')
        clone2hits = readNcbiXml(infile, options)
        separateClusters(clone2hits, options.sample2total, options.minExpFreq, options.minExpSize, options.minExpClones, options.minMotifClones, vjname, fh1, fh2, fh3, fh4, type2count, type2clones)
        
    fh1.close()
    fh2.close()
    fh3.close()
    fh4.close()

    #summary stats:
    summaryFile = "%s-summary.txt" %outbasename
    fh = open(summaryFile, 'w')
    fh.write("#min similarity: %f\n#min length: %d\n#min size to be called expanded clones: %d\n#min number of expanded clones required for type 1: %d\n#min number of clones supported a single motif for type 2: %d\n" %(options.minPos, options.minLen, options.minExpSize, options.minExpClones, options.minMotifClones))
    for type, count in type2count.iteritems():
        fh.write("%s\t%d\n" %(type, count))
    fh.close()

    #Print type2clones:
    v2seq = {}
    if options.vfile:
        v2seq = readVfile(options.vfile)
    type2index = {}
    for type, header2vj2seq in type2clones.iteritems():
        outfile = getFaFilename(type, outbasename)
        f = open(outfile, 'w')
        f.write("#%s\n" %type)
        for header, vj2seq in header2vj2seq.iteritems():
            for vj, seq in vj2seq.iteritems():
                #f.write(">%s;%s\n" %(header, vj))
                #reformat header:
                sample, id, size = getCloneInfo( ";".join(header.split(";")[:3]) ) 
                header = "%s;%s;%d;%s" % (sample.lstrip('as'), seq, size, vj)
                f.write(">%s\n" %header)

                v = vj.split(".")[0]
                if v in v2seq:
                    seq = v2seq[v] + seq
                f.write("%s\n" %seq)
        f.close()

if __name__ == '__main__':
    main()


