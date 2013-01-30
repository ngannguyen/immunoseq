#!/usr/bin/env python2.6

"""
Tue Dec  4 11:54:18 PST 2012
Parse Blast XML output file and cluster sequences using greedy approach.
Input: Blast xml file
Output: Text file, each line = 1 cluster, each element of a cluster is space-separated

Algorithm summary:
Sorted sequences by descending in size
Start with the largest sequence and use that as seed
For each sequence:
    Search for the closest seed that have >= %X similarity cutoff
    If found such seed: add the sequence to the seed's cluster
    else: the sequence becomes a seed of a new cluster

Cluster types:
    1/ Seed has multiple expanded clones as matches.
        1.1: The expanded clones are from the same sample with seed
        1.2: The expanded clones are from at least one sample different from seed sample
    2/ Seed has many small clones with similar motifs
        Cutoffs includes: a/ minimum number of clones contribute to one motif, b/ number of samples
        2.1: the clones carrying the motif are from the same sample with seed
        2.2: the clones carrying the motif are from at least one different sample than seed sample
    3/ Seed has no similar clones:
        3.1: seed is expanded
        3.2: seed is not expanded
    4/ Everything else (similar to type 2 but did not pass the cutoffs, i.e seed has a small number of low-frequency hits, or too many motifs but not enough clones to support a single motif)

"""
import os, re, sys
from Bio.Blast import NCBIXML
from optparse import OptionParser

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
        self.seq = ''
        self.hits = {} #key = hitCloneId, val = Hit

    def setSeq(self, seq):
        self.seq = seq

    def addHit(self, hitid, hit):
        self.hits[hitid] = hit

    def setFreq(self, total):
        if total == 0:
            raise ValueError("Error: Total sequences of sample %s is 0." %(self.sample))
        else:
            self.freq = 100.0*self.size/total

    def __cmp__(self, other):
        return cmp(self.size, other.size)

class Cluster():
    def __init__(self, seed):
        self.clones = [seed]
        self.totalReads = seed.size
        self.numClones = 1
        self.seed = seed
        self.motif2count = {}

    def addClone(self, clone):
        if clone not in self.clones:
            self.totalReads += clone.size
            self.numClones += 1
            self.clones.append(clone)

    def setType(self, type):
        self.type = type
            
    def setMotifs(self, motif2count):
        self.motif2count = motif2count

    def __cmp__(self, other):
        return cmp(self.totalReads, other.totalReads)

def typeid2desc(id):
    id2desc = { 1.1: "Multiple expanded clones from 1 sample", 
                1.2: "Multiple expanded clones from at least 2 samples",
                2.1: "Multiple non-expanded clones carrying the same motif, from 1 sample",
                2.2: "Multiple non-expanded clones carrying the same motif, from >=2 samples",
                3.1: "Clone with no hit, expanded",
                3.2: "Clone with no hit, non-expanded",
                4:   "Others"}
    return id2desc[id]

def isExpanded(clone, minSize, minFreq):
    if clone.size >= minSize and clone.freq >= minFreq: #expanded
        return True
    return False

def isSuper(motif1, motif2):
    #Return True if motif1 is a superset of motif2, otherwise return False
    if motif1 == motif2 or len(motif1) != len(motif2):
        return False

    for i, m1 in enumerate(motif1):
        if m1 != '.' and m1 != motif2[i]:
            return False
    return True

def getClusterType(seed2cluster, options):
    for seed, cluster in seed2cluster.iteritems():
        seedclone = cluster.seed

        if cluster.numClones == 1: #single element cluster
            if isExpanded( seedclone, options.minExpSize, options.minExpFreq ): 
                cluster.setType(3.1)
            else:
                cluster.setType(3.2)
        else: 
            numExp = 0
            expSamples = []
            motif2count = {}
            motif2samples = {}

            if isExpanded(seedclone, options.minExpSize, options.minExpFreq):
                numExp += 1
                expSamples.append(seedclone.sample)

            for hitclone in cluster.clones:
                if hitclone.desc != seed:
                    if isExpanded(hitclone, options.minExpSize, options.minExpFreq):
                        numExp += 1
                        if hitclone.sample not in expSamples:
                            expSamples.append(hitclone.sample)

                    if hitclone.desc not in seedclone.hits:
                        hit = hitclone.hits[seed]
                        #sys.stderr.write("Seed: %s. Hitclone: %s is not in the hits list: %s.\n" %(seed, hitclone.desc, "  ".join(seedclone.hits.keys())))
                    else:
                        hit = seedclone.hits[hitclone.desc]
                    motif = ""
                    for i, q in enumerate(hit.query):
                        s = hit.sbjct[i]
                        if q == s:
                            motif += q
                        else:
                            motif += "."
                    
                    #Search to see if any existing motif is a superset of current motif or if current motif is a super set of existing motif:
                    added = False
                    for prevmotif in motif2count.keys():
                        if motif == prevmotif or isSuper(prevmotif, motif):#prevmotif is a superset of current motif, update its count and don't add curr motif
                            motif2count[prevmotif] += 1
                            if hitclone.sample not in motif2samples[prevmotif]:
                                motif2samples[prevmotif].append(hitclone.sample)
                            added = True
                            break

                    if not added: #no prev motif is a super set of current motif
                        #check if current motif is a superset of prevmotif, add curr motif, remove previous motif
                        for prevmotif in motif2count.keys():
                            if isSuper(motif, prevmotif):
                                if motif not in motif2count:
                                    motif2count[motif] = motif2count[prevmotif]
                                    motif2samples[motif] = [ seedclone.sample ]
                                else:
                                    motif2count[motif] += motif2count[prevmotif]
                                for sample in motif2samples[prevmotif]:
                                    if sample not in motif2samples[motif]:
                                        motif2samples[motif].append(sample)
                                del motif2count[prevmotif]
                                del motif2samples[prevmotif]

                    if motif not in motif2count:
                        motif2count[motif] = 1
                        motif2samples[motif] = [ seedclone.sample ]
                    else:
                        motif2count[motif] += 1
                    if hitclone.sample not in motif2samples[motif]:
                        motif2samples[motif].append(hitclone.sample)

            if numExp >= options.minExpClones: #type 1
                if len(expSamples) == 1: #only the seed clone
                    type = 1.1
                else:
                    type = 1.2
            else:
                type = 4
                for motif, count in motif2count.iteritems():
                    if count >= options.minMotifClones:#type 2
                        if len( motif2samples[motif] ) == 1:
                            type = 2.1
                        else:
                            type = 2.2
                        break
            cluster.setType(type)
            cluster.setMotifs(motif2count)

def getFh(type, fh11, fh12, fh21, fh22, fh31, fh32, fh4):
    if type == 1.1:
        return fh11
    elif type == 1.2:
        return fh12
    elif type == 2.1:
        return fh21
    elif type == 2.2:
        return fh22
    elif type == 3.1:
        return fh31
    elif type == 3.2:
        return fh32
    else:
        return fh4

def printClusters( outfile, seed2cluster ):
    clusters = sorted( seed2cluster.values(), reverse = True )
    fh = open(outfile, 'w')
    outbasename = outfile.rstrip('txt').rstrip('.')
    vjname = os.path.basename(outbasename)

    fh11 = open( "%s_1.1" % outbasename, 'w' )
    fh12 = open( "%s_1.2" % outbasename, 'w' )
    fh21 = open( "%s_2.1" % outbasename, 'w' )
    fh22 = open( "%s_2.2" % outbasename, 'w' )
    fh31 = open( "%s_3.1" % outbasename, 'w' )
    fh32 = open( "%s_3.2" % outbasename, 'w' )
    fh4 = open( "%s_4" % outbasename, 'w' )

    totalClones = 0
    for i, cluster in enumerate(clusters):
        clones = sorted( cluster.clones, key= lambda c:c.size, reverse=True)
        fh.write( "%s\n" %(" ".join([c.desc for c in clones])) )
        
        totalClones += cluster.numClones
        fh_long = getFh( cluster.type, fh11, fh12, fh21, fh22, fh31, fh32, fh4 )
        fh_long.write(">Cluster %d, type %.1f, %s, %d clones, %d totalReads, motifs: %s\n" %(i, cluster.type, vjname, cluster.numClones, cluster.totalReads, ";".join(["%s_%d" %(m,c)  for m,c in cluster.motif2count.iteritems()]) ))
        for c in clones:
            fh_long.write("\t%s\t%s\t%f\n" %(c.seq, c.desc, c.freq))
        fh_long.write("\n")
    #fh_long.write("\nTotal clones: %d\nTotal clusters: %d\n" %(totalClones, len(clusters) ))

    fh.close()
    fh11.close()
    fh12.close()
    fh21.close()
    fh22.close()
    fh31.close()
    fh32.close()
    fh4.close()

def getClusters(clones):
    seed2cluster = {} #key = seed, val = list of clones that cluster to the seed
    if len(clones) == 0:
        return seed2cluster
    
    #First seed:
    firstCluster = Cluster(clones[0])
    seed2cluster[ clones[0].desc ] = firstCluster

    if len(clones) == 1:
        return seed2cluster

    for clone in clones[1:]:
        maxPos = 0
        maxSize = 0
        bestSeed = ''

        for seed in seed2cluster:
            if seed in clone.hits:
                hit = clone.hits[seed]
                #if hit.positives > maxPos:
                if hit.identities > maxPos:
                    #maxPos = hit.positives
                    maxPos = hit.identities
                    maxSize = Clone(seed).size
                    bestSeed = seed
                #elif hit.positives == maxPos:
                elif hit.identities == maxPos:
                    currSize = Clone(seed).size
                    if currSize > maxSize:
                        maxSize = currSize
                        bestSeed = seed

        if bestSeed != '':
            seed2cluster[ bestSeed ].addClone( clone ) #add to the cluster of the closest seed
        else:
            cluster = Cluster(clone)
            seed2cluster[ clone.desc ] = cluster #create new cluster
    return seed2cluster

def readNcbiXml(infile, minLen, minPos, sample2total, sameLen):
    rh = open(infile)
    records = NCBIXML.parse(rh)
    clones = []

    for record in records:
        if record.query_length < minLen: #too short, pass
            continue
        clone = Clone(record.query)

        if sample2total:
            clone.setFreq(sample2total[clone.sample])
        for aln in record.alignments:
            for hit in aln.hsps: # each hit
                
                if len(hit.match) < record.query_length: #ignore local alignment
                    continue
                if clone.seq == '':
                    clone.setSeq(hit.query)

                if sameLen and (re.search('-', hit.query) or re.search('-', hit.sbjct)): #don't allow for gap alignment if sameLen is specifie
                    continue
                #if float(hit.positives)/len(hit.query) < minPos: #low similarity hit, ignore
                if float(hit.identities)/len(hit.query) < minPos: #low identity hit, ignore
                    continue
                hitid = aln.title.split()[-1]
                if hitid == clone.desc: #self alignment, ignore
                    continue
                clone.addHit(hitid, hit)
        clones.append(clone)
    
    #if sample2total:
    #    clones = sorted( clones, key=lambda c:c.freq, reverse=True )
    #else:
    #    clones = sorted( clones, key=lambda c:c.size, reverse=True )
    return clones

def readSample2total(file):
    sample2total = {}
    f = open(file, 'r')
    for line in f:
        items = line.strip().split()
        sample2total[items[0]] = int(items[1])
    f.close()
    return sample2total

def getfiles(indir, ext):
    files = []
    for file in os.listdir(indir):
        items = file.split('.')
        if items[-1] == ext:
            files.append(file)
    return files

def getInfiles(input):
    ext = 'xml'
    if os.path.isdir(input): #input is a directory
        infiles = getfiles(input, ext)
        infiles = [ os.path.join(input, f) for f in infiles ]
    else:
        infiles = [input]
    return infiles

def addOptions(parser):
    parser.add_option('-i', '--input', dest='input', help='Input file or directory')
    parser.add_option('-o', '--outfile', dest='outfile', help='Output file')
    parser.add_option('-p', '--positive', dest='minPos', type='float', default=0.9, help='Minimum portion of positive matches. Default=%default')
    parser.add_option('-l', '--len', dest='minLen', type='int', default=10, help='Minimum sequence length to be included in the output. Default=%default')
    parser.add_option('-L', '--lenRestriction', dest='sameLen', action='store_true', default=False, help='If specified, only sequences of same length can be clustered together. Default=%default')
    parser.add_option('-S', '--minExpandedSize', dest='minExpSize', type='int', default=1000, help='Minimum number of reads for a clone to be called "expanded". Default=%default')
    parser.add_option('-F', '--minExpandedFreq', dest='minExpFreq', type='float', default=0.0, help='Minimun frequency for a clone to be called "expanded".Range from 0 - 100. Default=%default')
    parser.add_option('-C', '--minExpandedClones', dest='minExpClones', type='int', default=1, help='Minimum number of similar expanded clones for a seed and its cluster to be classified as type 1. Default=%default')
    parser.add_option('-c', '--minMotifClones', dest='minMotifClones', type='int', default=10, help='Minimum number of clones carrying the same motif for a seed and its cluster to be classified as type 2. Default=%default ')
    parser.add_option('--sample2total', dest='sample2total', help='Required if --minExpandedFreq is larger than 0. Format: <sample> <totalCount>')
    #parser.add_option('-v', '--addV', dest='vfile', help='If specified, add the rest of V gene to each sequence in the output fasta files. Default = None')

def main():
    usage = "usage: %prog [options]\n"
    parser = OptionParser( usage=usage )
    addOptions(parser)
    options, args = parser.parse_args()
    
    if options.minExpFreq > 0 and not options.sample2total:
        parser.error("--sample2total is required as --minExpandedFreq > 0\n")
    if options.sample2total:
        options.sample2total = readSample2total(options.sample2total)
    
    #Read input XML file(s):
    infiles = getInfiles(options.input)
    clones = []
    for infile in infiles:
        currclones = readNcbiXml(infile, options.minLen, options.minPos, options.sample2total, options.sameLen)
        clones.extend( currclones )
    if options.sample2total:
        clones = sorted( clones, key=lambda c:c.freq, reverse=True )
    else:
        clones = sorted( clones, key=lambda c:c.size, reverse=True )
    sys.stderr.write("Done reading input file. Total %d clones passed minLen\n" %len(clones))
    #Done reading input XML file(s)

    seed2cluster = getClusters(clones)
    getClusterType(seed2cluster, options)
    printClusters( options.outfile, seed2cluster )

if __name__ == '__main__':
    main()

