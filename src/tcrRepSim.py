#!/usr/bin/env python

#nknguyen at soe ucsc edu
#Feb 07 2012
#TCR repertoire simulation

'''
Simulate CDR3 TCR repertoires to assess the overlapping of clones/sequences between two independent samplings of certain size.
(Purpose: use to determine what is the sufficient sampling size to get all the important, clonally expanded CDR3)

There are 4 main steps involved:
1/ Generate a set of m unique CDR3 nucleotide sequences, call this set R (it is estimated that each individual has about 3x10^6 uniq seqs)
   This steps take into account the CDR3-length distribution and the Amino-acid usage at each position of the CDR3.
   Input: a/ len2count.txt with format: <aa-cdr3Length>\t<count>
          b/ aaUsage.txt with format: <aa-cdr3Length>\t<Position>\t<aa-Letter>\t<count>
          c/ m = size of the repertoire to be created
   Generate: m nucleotide sequences such that when translated to amino acid sequences, follows the input length distribution and aaUsage.

2/ Generate the repertoire (say P) of CDR3 of an individual from set R in (1):
   Let M be the number of cells (number of total sequences) an individual has. It is estimated that there are about 6x10^7 CD8+ T cells per individual.
   Populate repertoire P with the n sequences in R until P has M sequences. This step can assume a uniform distribution of the uniq sequences or model clonal expansion.
   The program model clonal expansion by in input pseudo repertoire with format <topFreq>,<secondTopFreq>,...,<ith-topFreq>. A random sequence will be drawn from n sequences in (a), and will populate P with topFreq*M counts. Then a second random seq will be drawn, etc until there are M sequences or the frequency list is done. If freq list is done and there are < M sequences in P, just randomly fill in P with the rest of the m sequences
   Generate: M CDR3 nucleotide sequences

3/ Sampling:
   Repeat step (1) and (2) to generate another TCR repertoire for the second individual.
   For each individual, randomly choose S sequences from the corresponding repertoire.
   Calculate how much overlap the two samplings have.
   Repeat this step with difference sampling size

4/ a. Repeat step (3) -numSamplingPerSim times. Calculate mean & std
(Or/And)
   b. Repeat steps (1,2,3) N times (N simulations). Calculate mean statistics

(5/ Plots)
'''

import os, sys, re, time, random, copy, gzip
import cPickle as pickle
from optparse import OptionParser
import xml.etree.ElementTree as ET

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import getTempDirectory
from sonLib.bioio import setLogLevel
import numpy as np

#from immunoseq.lib import *

################## SEQUENCE OBJECTS ######################
class Seq:
    def __init__(self, nuc, aa):
        self.nuc = nuc
        self.aa = aa
        self.count = 0

    def __cmp__(self, other):
        return cmp(self.nuc, other.nuc)

    def setCount(self, count):
        self.count = count

    def updateCount(self, add):
        self.count += add

class Aaseq:
    def __init__(self, aa):
        self.aa = aa
        self.count = 0

    def __cmp__(self, other):
        return cmp(self.aa, other.aa)

    def setCount(self, count):
        self.count = count

    def updateCount(self, add):
        self.count += add

    def setFreq(self, total):
        if total == 0:
            self.freq = 0
        else:
            self.freq = 100.0*self.count/total

class Seqs(list):
    def add(self, seq):
        leftIndex = 0
        rightIndex = len(self) -1
        
        if len(self) == 0 or seq > self[rightIndex]:
            self.append(seq)
        elif seq < self[leftIndex]:
            self.insert(0, seq)
        elif seq == self[leftIndex] or seq == self[rightIndex]:
            return
        else:
            while True:
                middleIndex = int( (leftIndex + rightIndex)/2 )
                compare = cmp( self[middleIndex], seq )
                if compare == 0: #sequence is already in the list, do not add, just return
                    break
                elif leftIndex == middleIndex: #add seq to the right of leftIndex
                    self.insert(rightIndex, seq)
                    break
                elif compare > 0: #seq lies somewhere btw [middleIndex, rightIndex]
                    rightIndex = middleIndex
                elif compare < 0: 
                    leftIndex = middleIndex
            
    def search(self, seq):
        leftIndex = 0
        rightIndex = len(self) -1
        
        if len(self) == 0 or seq > self[rightIndex] or seq < self[leftIndex]: #not in list
            return -1
        elif seq == self[leftIndex]:
            return leftIndex
        elif seq == self[rightIndex]:
            return rightIndex
        else:
            while True:
                middleIndex = int( (leftIndex + rightIndex)/2 )
                compare = cmp( self[middleIndex], seq )
                if compare == 0: #sequence is already in the list, do not add, just return
                    return middleIndex
                elif leftIndex == middleIndex: #Not found
                    return -1
                elif compare > 0: #seq lies somewhere btw [middleIndex, rightIndex]
                    rightIndex = middleIndex
                elif compare < 0: 
                    leftIndex = middleIndex

######################### MAIN PIPELINE ################################

#============= SETTING UP SIMULATIONS ==============
class Setup( Target ):
    """Setting up simulations
    """
    def __init__(self, options, doneSims):
        Target.__init__(self, time = 0.00025)
        self.options = options
        self.doneSims = doneSims

    def run(self):
        setLogLevel("DEBUG")
        numSim = self.options.numSim
        maxSimPerRun = self.options.maxSimPerRun
        sims = min( [maxSimPerRun, numSim - self.doneSims] )
        
        outdir = os.path.join( self.options.outdir, "sims", str(self.doneSims) ) #outdir/sims/batchId
        system("mkdir -p %s" %outdir)
        self.addChildTarget( SimulationBatch(self.options, outdir, sims, self.doneSims) )
        
        doneSims = self.doneSims + sims
        if doneSims < numSim:
            self.setFollowOnTarget( Setup(self.options, doneSims) ) #recusion call for next batch of simulations
        else:
            summarydir = os.path.join(self.options.outdir, "summary") #outdir/summary
            statsdir = os.path.join(self.options.outdir, "stats") #outdir/stats
            system("mkdir -p %s" % statsdir)
            readPickle = True
            writePickle = False
            writeSummary = True
            self.setFollowOnTarget( Summary(summarydir, self.options.samSize, self.options.numSamples, statsdir, writePickle, writeSummary, readPickle) )

class SimulationBatch( Target ):
    def __init__(self, options, outdir, numSim, batchid):
        Target.__init__(self, time=0.00025)
        self.options = options
        self.outdir = outdir #outdir/sims/batchId
        self.numSim = numSim
        self.batchid = batchid

    def run(self):
        for i in xrange(self.numSim):
            outdir = os.path.join(self.outdir, "sim-%d" %i) #outdir/sims/batchId/sim-Id
            system("mkdir -p %s" %outdir)
            self.addChildTarget( Simulation(self.options, outdir, self.batchid, i) )
         
        #Combine the results
        samplingsDir = os.path.join(self.options.outdir, "samplings", str(self.batchid)) #outdir/samplings/batchid
        sumdir = os.path.join(self.options.outdir, "summary", str(self.batchid)) #outdir/summary/batchid
        system("mkdir -p %s" %sumdir) 
        readPickle = False
        writePickle = True
        writeSummary = False
        #self.setFollowOnTarget( Summary(self.outdir, self.options.samSize, self.options.numSamples, sumdir, writePickle, writeSummary, readPickle) )
        self.setFollowOnTarget( Summary(samplingsDir, self.options.samSize, self.options.numSamples, sumdir, writePickle, writeSummary, readPickle) )

#------------- SIMULATION ---------------------
class Simulation( Target ):
    """Generate s simulations (corresponding to s samples (s individuals))
       Creating TCR repertoire for each individuals
       After got repertoires for All individuals, call following job to sampling from these repertoires
    """
    def __init__(self, options, outdir, batchid, simid):
        Target.__init__(self, time=0.00025)
        self.options = options
        self.outdir = outdir #outdir/sims/batchId/sim-Id
        self.batchid = batchid
        self.simid = simid
    
    def run(self):
        globalTempDir = self.getGlobalTempDir()
        for i in xrange(self.options.numSamples):
            outdir = os.path.join(globalTempDir, str(i))
            system("mkdir -p %s" %outdir) #simTempDir/sampleID/
            self.addChildTarget( SimulationSingle(self.options, outdir) )
        
        #Sampling from these repertoire:
        samplingsOutdir = os.path.join(self.options.outdir, "samplings", str(self.batchid), "sim-%d" %self.simid)
        system("mkdir -p %s" %samplingsOutdir) #outdir/samplings/batchId/sim-id
        self.setFollowOnTarget( Samplings(globalTempDir, self.outdir, self.options.numSamples, self.options.samSize, self.options.cutoffs, self.options.numSamplingPerSim, samplingsOutdir) )
        #Calculate pair-wise overlap between the samples 
        #self.setFollowOnTarget( Overlap(globalTempDir, self.outdir, self.options.numSamples, self.options.samSize, self.options.cutoffs) )

class SimulationSingle( Target ):
    """Generate repertoire for one sample/individual
       Starts by generating the uniq sequences of certain CDR3 length
       Then calls follow on jobs to use these sequences to create the repertoire
    """
    def __init__(self, options, outdir):
        Target.__init__(self, time = 0.00025)
        self.outdir = outdir #simTempDir/sampleID
        self.options = options

    def run(self):
        globalTempDir = self.getGlobalTempDir()
        #Generate repertoire of uniq sequences for each length:
        for l in self.options.len2freq:#each length
            size = self.options.totalClones*self.options.len2freq[l]
            if l not in self.options.len2aaUsage or size == 0:
                continue
            aaUsage = self.options.len2aaUsage[l]
            aa2codons = self.options.aa2codons
            self.addChildTarget( UniqSeqs(globalTempDir, size, aaUsage, aa2codons) )
        
        self.setFollowOnTarget( Repertoire(globalTempDir, self.outdir, self.options.totalSeqs, self.options.topFreqs, self.options.samSize) )

#+++++++++++ REPERTOIRE OF UNIQ SEQS/ CLONES ++++++++++++++
class UniqSeqs( Target ):
    """Generate repertoire of uniq sequences of a certain length
    """
    def __init__(self, outdir, size, aaUsage, aa2codons):
        Target.__init__(self)
        self.outdir = outdir
        self.size = size
        self.aaUsage = aaUsage
        self.aa2codons = aa2codons
    
    def run(self):
        #Generate the sequences
        seqs = getUniqSeqRep(self.size, self.aaUsage, self.aa2codons)
        #Pickling the sequences to outdir
        pickleFile = os.path.join( self.outdir, "uSeqs-l%d.pickle" %len(self.aaUsage) )
        pickle.dump(seqs, gzip.open(pickleFile, "wb"))

#+++++++++++++ PSEUDO REPERTOIRE ++++++++++++++++++
class Repertoire( Target ):
    """Combine the len-specific uniqSeqs repertoires
       From uniqSeq repertoire, generate a pseudo repertoire, using topFreqs if available
       Add jobs to run samplings with different sizes
    """
    def __init__(self, indir, outdir, size, topFreqs, samSizes):
        Target.__init__(self)
        self.indir = indir #simSampleTempDir
        self.outdir = outdir #simTempDir/sampleID
        self.size = size
        self.topFreqs = topFreqs
        self.samSizes = samSizes

    def run(self):
        #Load pickles to get uniqSeqs
        pickleFiles = os.listdir( self.indir )
        seqs = []
        for file in pickleFiles:
            if re.search("pickle", file):
                seqs.extend( pickle.load( gzip.open( os.path.join(self.indir, file) , "rb") ) )
                system( "rm -f %s" %(os.path.join(self.indir, file)) ) #Remove the uniq-sequences pickle file of the specific length after done loading

        getRep( seqs, self.size, self.topFreqs )
        
        #Print repertoire to server so that it can be read by following jobs:
        pickleFile = os.path.join(self.outdir, "rep.pickle") #simTempDir/sampleID/rep.pickle
        pickle.dump( seqs, gzip.open(pickleFile, "wb") )

#++++++++++++++ SAMPLING from pseudoRepertoire ++++++++++++++++
class Samplings( Target ):
    """(For each repertoier simulation,) Set up samplings for all pairs, all sampling sizes, all number of samplings.
    """
    def __init__(self, indir, outdir, numSamples, samSizes, cutoffs, numSamplings, samplingsOutdir):
        Target.__init__(self, time=0.00025)
        self.indir = indir #simTempDir
        self.outdir = outdir #outdir/sims/batchId/sim-Id
        self.numSamples = numSamples
        self.samSizes = samSizes
        self.cutoffs = cutoffs
        self.numSamplings = numSamplings
        self.samplingsOutdir = samplingsOutdir #outdir/samplings/batchId/sim-id

    def run(self):
        for i in xrange(self.numSamplings):#Each sampling
            for samsize in self.samSizes: #Each sampling size
                outdir = os.path.join(self.outdir, "sampling-%d" %i) #outdir/sims/batchId/sim-Id/sampling-Id
                for j in xrange(self.numSamples - 1):
                    for k in xrange(self.numSamples): #Each pair
                        pair = [j, k]
                        self.addChildTarget( SamplingPair(self.indir, outdir, samsize, pair, self.cutoffs) )
        #Summary the overlapping stats:
        readPickle = False
        writePickle = False
        writeSummary = True
        self.setFollowOnTarget( Summary(self.outdir, self.samSizes, self.numSamples, self.samplingsOutdir, writePickle, writeSummary, readPickle) )

class SamplingPair( Target ):
    """
    """
    def __init__(self, indir, outdir, size, pair, cutoffs):
        Target.__init__(self, time=0.00025)
        self.indir = indir #simTempDir (which has simTempDir/sampleId/rep.pickle)
        self.outdir = outdir #outdir/sims/batchId/sim-Id/sampling-Id
        self.size = size
        self.pair = pair
        self.cutoffs = cutoffs

    def run(self):
        globalTempDir = self.getGlobalTempDir()
        seqfiles = []
        for i, p in enumerate(self.pair):
            repertoireFile = os.path.join(self.indir, str(p), "rep.pickle") #simTempDir/sampleId/rep.pickle
            outfile = os.path.join(globalTempDir, "%d.pickle" %i) #samplingPairTempDir/sampleId.pickle
            seqfiles.append(outfile)
            self.addChildTarget( Sampling(outfile, repertoireFile, self.size) )
        
        overlapOutdir = os.path.join( self.outdir, "%d-%d" %(self.pair[0], self.pair[1]) ) #outdir/sims/batchId/sim-Id/sampling-Id/sam1-sam2
        system("mkdir -p %s" %overlapOutdir)
        overlapOutfile = os.path.join( overlapOutdir, "%d.txt" %self.size )
        self.setFollowOnTarget( OverlapPairwise(seqfiles, overlapOutfile, self.cutoffs) )

class Sampling( Target ):
    """Sample from the input repertoire for "size" sequences
    """
    def __init__(self, outfile, repertoireFile, samsize):
        Target.__init__(self)
        self.outfile = outfile
        self.repfile = repertoireFile
        self.size = samsize

    def run(self):
        seqs = pickle.load( gzip.open(self.repfile, "rb") ) 
        samseqs = sampling(seqs, self.size)

        #Pickling seqs:
        pickle.dump( samseqs, gzip.open(self.outfile, "wb") )

#----------- OVERLAP -----------
class OverlapPairwise( Target ):
    """Output file: %cutoff\tclones1\tclones2\toverlap\t%1overlap2\t%2overlap1\t%reads1overlap2\t%reads2overlap1
    """
    def __init__(self, seqsFiles, outfile, cutoffs):
        Target.__init__(self, time=0.001)
        self.outfile = outfile
        self.seqsFile1 = seqsFiles[0]
        self.seqsFile2 = seqsFiles[1]
        self.cutoffs = cutoffs

    def run(self):
        seqs1 = pickle.load( gzip.open(self.seqsFile1, "rb") )
        seqs2 = pickle.load( gzip.open(self.seqsFile2, "rb") )
        reads1, reads2, clones1, clones2, stats1, stats2 = getOverlap(seqs1, seqs2, self.cutoffs)
        
        #Print stats
        f = open(self.outfile, "w")
        f.write("#%Cutoff\tClones1\tClones2\tOverlap1\tOverlap2\t%1overlap2\t%2overlap1\t%reads1overlap2\t%reads2overlap1\n")
        for i,c in enumerate(self.cutoffs):
            oc1 = stats1["oclones"][i]
            or1 = stats1["oreads"][i]
            t1 = clones1[i]
            r1 = reads1[i]
            
            oc2 = stats2["oclones"][i]
            or2 = stats2["oreads"][i]
            t2 = clones2[i]
            r2 = reads2[i]
            
            f.write("%.3f\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n" %( c, t1, t2, oc1, oc2, getPc(oc1, t1), getPc(oc2, t2), getPc(or1, r1), getPc(or2, r2) ))
        f.close()

#================ AVERAGE ACROSS THE SIMULATIONS ====================
class Summary( Target ):
    def __init__(self, indir, samSizes, numSamples, outdir, writePickle, writeSummary, readPickle):
        Target.__init__(self, time=0.001)
        self.indir = indir #outdir/sims/batchId/sim-Id; #outdir/samplings/batchId; #outdir/summary
        self.samSizes = samSizes
        self.numSamples = numSamples
        self.outdir = outdir #outdir/samplings/batchId/sim-Id; #outdir/summary/batchid; #outdir/stats
        self.writePickle = writePickle
        self.writeSummary = writeSummary
        self.readPickle = readPickle

    def run(self):
        #Get the pairs:
        for i in xrange(self.numSamples - 1):
            for j in xrange(self.numSamples):
                pair = "%d-%d" %( i, j )
                outdir = os.path.join(self.outdir, pair) #outdir/samplings/batchId/sim-Id/pair; #outdir/summary/batchid/pair 
                system("mkdir -p %s" %outdir)
                for samsize in self.samSizes:
                    self.addChildTarget( SummarySamsize(self.indir, samsize, pair, outdir, self.writePickle, self.writeSummary, self.readPickle) )
        if not self.writePickle:
            self.setFollowOnTarget( CleanupSummary(self.indir) )

class SummarySamsize( Target ):
    def __init__(self, indir, samsize, pair, outdir, writePickle, writeSummary, readPickle):
        Target.__init__(self, time=0.001)
        self.indir = indir #outdir/sims/batchId/sim-Id; #outdir/samplings/batchId; #outdir/summary
        self.samsize = samsize
        self.pair = pair
        self.outdir = outdir #outdir/samplings/batchId/sim-Id/pair; #outdir/summary/batchid/pair; #outdir/stats/pair
        self.writePickle = writePickle
        self.writeSummary = writeSummary
        self.readPickle = readPickle

    def run(self):
        #Get the stats:
        cutoff2stats = {} #key = cutoff, val = [clones1, clones2, overlap1, overlap2, %1oc2, %2oc1, %1or2, %2or1]
        indirs = os.listdir(self.indir)
        for indir in indirs:
            if not os.path.isdir( os.path.join(self.indir,indir) ):
                continue
            file = os.path.join(os.path.join(self.indir,indir), self.pair, "%d.txt" %self.samsize)
            if self.readPickle:
                file = os.path.join(os.path.join(self.indir,indir), self.pair, "%d.pickle" %self.samsize)
            
            if not os.path.exists(file):
                raise NonExistFileError("File %s does not exist\n" %file)
            
            if not self.readPickle:
                currCutoff2stats = readStats(file)
            else:
                currCutoff2stats = pickle.load( gzip.open(file, "rb") )
            for c, stats in currCutoff2stats.iteritems():
                if c not in cutoff2stats:
                    cutoff2stats[c] = stats
                else:
                    for i, s in enumerate(stats):
                        cutoff2stats[c][i].extend(s)
        
        #Calculate average and standard deviation:
        if self.writePickle:
            outfile = os.path.join(self.outdir, "%d.pickle" %self.samsize)
            pickle.dump(cutoff2stats, gzip.open(outfile, "wb"))
        
        if self.writeSummary:
            for c, stats in cutoff2stats.iteritems():
                for i, s in enumerate(stats):
                    cutoff2stats[c][i] = [np.mean( s ), np.std( s )]
                    #cutoff2stats[c][i] = float(s)/self.denom

            #Print to output file:
            outfile = os.path.join(self.outdir, "%d.txt" %self.samsize)
            f = open(outfile, 'w')
            f.write("#%Cutoff\tClones1\tClones2\toverlap1\toverlap2\t%1overlap2\t%2overlap1\t%reads1overlap2\t%reads2overlap1\t")
            f.write("StdClones1\tStdClones2\tStdOverlap1\tStdOverlap2\tStd%1overlap2\tStd%2overlap1\tStd%reads1overlap2\tStd%reads2overlap1\n")
            for c in sorted( cutoff2stats.keys() ):
                stats = cutoff2stats[c]
                meanStats = [ "%.4f" % s[0] for s in stats ]
                stdStats = [ "%.4f" % s[1] for s in stats ]
                f.write( "%.3f\t%s\t%s\n" %(c, "\t".join(meanStats), "\t".join(stdStats)) )
            f.close()

#=============== CLEANUP =================
class CleanupSummary ( Target ):
    def __init__(self, indir):
        Target.__init__(self, time = 0.00025)
        self.indir = indir

    def run(self):
        system("rm -fR %s" %self.indir)

####################### UTILITIES FUNCTIONS ##########################

#============== calculate overlap =================
def getOverlap( seqs1, seqs2, cutoffs ):
    #Initialize stats:
    stats1 = {"oclones":[], "oreads":[]}
    stats2 = {"oclones":[], "oreads":[]}
    for i in xrange( len(cutoffs) ):
        for k in stats1:
            stats1[k].append(0)
            stats2[k].append(0)

    #get number of clones that pass the cutoff:
    reads1, clones1, total1 = getNumClones(seqs1, cutoffs)
    reads2, clones2, total2 = getNumClones(seqs2, cutoffs)

    if total1 == 0 or total2 == 0:
        return
    #get overlap:
    for s1 in seqs1:
        i2 = seqs2.search(s1)
        if i2 == -1:#not found in repertoire 2
            continue
        s2 = seqs2[i2]
        for i, cutoff in enumerate(cutoffs):
            #if s1.freq >= cutoff and s2.freq >= cutoff:
            if s1.freq >= cutoff:
                stats1["oclones"][i] += 1
                stats1["oreads"][i] += s1.count
            if s2.freq >= cutoff:
                stats2["oclones"][i] += 1
                stats2["oreads"][i] += s2.count 
    return reads1, reads2, clones1, clones2, stats1, stats2

def getNumClones(seqs, cutoffs):
    reads = [ 0 for c in cutoffs ] 
    clones = [ 0 for c in cutoffs ]
    total = sum([ s.count for s in seqs ])
    for s in seqs:
        s.setFreq(total)
        for i,c in enumerate(cutoffs):
            if s.freq >= c:
                clones[i] += 1
                reads[i] += s.count
    return reads, clones, total

#================ GENERATE PSEUDO REPERTOIRE =========
def sampling( seqs, size ):
    samseqs = Seqs()
    indexList = []
    for i, s in enumerate(seqs):
        indexList.extend( [ i for j in xrange( int(s.count)) ] )
    
    for j in xrange(size):
        i = random.randint(0, len(indexList) -1) #randomly pick one sequence from the repertoire
        s = Aaseq( seqs[indexList[i]].aa )
        #s = copy.copy( seqs[indexList[i]] )
        sindex = samseqs.search(s) #Search for this sequence in the current set of picked sequences
        
        #print [ seq.aa for seq in samseqs]
        #print "%s" %s.aa
        #print sindex
        
        if sindex >=0:
            samseqs[sindex].updateCount(1)
        else:#sequence hasn't picked yet
            s.setCount(1)
            samseqs.add(s)
    return samseqs

def getRep( seqs, size, topFreqs ):
    #Generate repertoire:
    topindices = []
    sizeToFill = size
    if topFreqs:
        for i, f in enumerate(topFreqs):
            index = random.randint(0, len(seqs) -1)
            while (index in topindices):
                index = random.randint(0, len(seqs) -1)
            topindices.append(index)
            count = int(f*size)
            seqs[index].setCount( count )
        sizeToFill = size - int( sum(topFreqs)*size )
    
    for i in xrange(sizeToFill):
        index = random.randint(0, len(seqs) -1)
        while (index in topindices):
            index = random.randint(0, len(seqs) -1)
        seqs[index].updateCount(1)

def getUniqSeqRep(size, aaUsage, aa2codons):
    """Generate "size" number of unique nucleotide sequences with length (len(aaUsage))
    Using the amino-acid Usage from aaUsage, and codon usage from aa2codons
    """
    seqs = Seqs() #list of sequences
    while len(seqs) < size:
        aa = getAaSeq(aaUsage)
        nuc = aa2nuc(aa, aa2codons)
        seq = Seq(nuc, aa)
        seqs.add(seq)
    return seqs

def readStats( file ):
    cutoff2stats = {} #key = cutoff, val = [clones1, clones2, overlap, %1oc2, %2oc1, %1or2, %2or1]
    f = open(file, 'r')
    for line in f:
        if line[0] == "#":
            continue
        items = line.split("\t")
        if len(items) < 9:
            raise FileFormatError("Wrong file format, 9 fields are required, see %d. File %s\n" %(len(items), file))
        cutoff = float(items[0])
        #stats = [int(items[1]), int(items[2]), int(items[3]), int(items[4]), float(items[5]), float(items[6]), float(items[7]), float(items[8])]
        stats = [ [float(items[i])] for i in xrange(1,9) ]
        #stats = [float(items[1]), float(items[2]), float(items[3]), float(items[4]), float(items[5]), float(items[6]), float(items[7]), float(items[8])]
        cutoff2stats[cutoff] = stats
    f.close()
    return cutoff2stats

def getPc(count, total):
    if total == 0:
        return 0
    return 100.0*count/total

def getAaSeq(aaUsage):
    aa = ""
    seqlen = len(aaUsage)
    for i in xrange(seqlen):
        letters = aaUsage[i]
        j = random.randint(0, len(letters) -1)
        aa += letters[j] 
    return aa

def aa2nuc(aaseq, aa2codons):
    nucseq = ''
    for aa in aaseq:
        codons = aa2codons[aa]
        i = random.randint( 0, len(codons) -1 )
        nucseq += codons[i]
    return nucseq

#================ ERROR CLASSES ============
class FileFormatError(Exception):
    pass

class TopFreqError(Exception): 
    pass

class NonExistFileError(Exception):
    pass

#================ READ INPUT FILES ================
def readCodonUsage(file):
    """
    """
    aa2codons = {} #key = aa, val = {codon:freq}
    aa2codon2freq = {}
    f = open(file, 'r')
    for line in f:
        line.strip()
        if line == "" or line[0] == "#":
            continue
        items = line.split()
        if len(items) %4 > 0:
            raise FileFormatError("Wrong codonUsage format\n")
        for i in xrange(0, len(items), 4):
            codon = items[i]
            aa = items[i+1]
            freq = float( items[i+2] )
            if aa not in aa2codon2freq:
                aa2codon2freq[aa] = {codon: freq}
            else:
                aa2codon2freq[aa][codon] = freq
    f.close()
    
    #Normalize
    normlen = 100
    for aa in aa2codon2freq:
        aa2codons[aa] = []
        for codon, freq in aa2codon2freq[aa].iteritems():
            for i in xrange( int(freq*normlen) ):
                aa2codons[aa].append(codon)

    return aa2codons

def readAaUsage(file):
    '''
    Return: len2aaUsage = {} #key = length, val= position-array, each item represent a.a usage of that position.
    Example: 
    9 1 C 10
    9 2 A 7
    9 2 R 1
    9 2 S 3
    9 2 V 1
    len2aaUsage = {9:["CCCCCCCCCC", "AAAAAAARSSSV"]}
    '''
    len2aaUsage = {} #key = length, val= position-array, each item represent a.a usage of that position.
    len2listLetter2freq = {}
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if line == "" or line[0] == "#":
            continue
        items = line.split()
        if len(items) != 4:
            raise FileFormatError("Wrong aaUsage file format. Required 4 fields, have %d. Line: %s\n" %(len(items), line) )
        l = int( items[0] )
        pos = int( items[1] ) - 1
        letter = items[2]
        count = int( items[3] )

        if l < 0 or pos < 0 or count < 0:
            raise ValueError("readAaUsage: length, position, count must >=0")
        if pos >= l:
            raise FileFormatError("readAaUsage: Column 2 value (Position) must <= col 1 value (length)\n")
        if l not in len2listLetter2freq:
            len2listLetter2freq[ l ] = [{} for i in xrange(l)]
        len2listLetter2freq[ l ][ pos ][ letter ] = count
    f.close()

    #Make sure that for each length, usage for all positions are provided from input file:
    for l in len2listLetter2freq:
        for i in xrange(l):
            if len( len2listLetter2freq[l][i].keys() ) == 0:
                raise FileFormatError("readAaUsage requires input file to have AA usage for all positions of each length. Could not find usage for length %d, position %d." %(l, i + 1))

    #Normalize counts:
    normLen = 10000 #equivalent the number of digits of the frequencies got taken into account
    for l in len2listLetter2freq:
        posList = len2listLetter2freq[l]
        len2aaUsage[ l ] = []
        for i, letter2count in enumerate(posList):
            total = sum([ letter2count[letter] for letter in letter2count ])
            if total > 0:
                for letter, count in letter2count.iteritems():
                    letter2count[letter] = normLen*count/total
            posLetters = "".join([ letter*count for letter, count in letter2count.iteritems()])
            len2aaUsage[l].append(posLetters)

    return len2aaUsage

def getUniAaUsage(lens):
    aas = "ARNDCEQGHILKFPSTWYV"
    len2aaUsage = {}
    for l in lens:
        if not isinstance(l, int):
            raise ValueError("getUniAaUsage, input lengths contain non-integer items\n")
        len2aaUsage[l] = []
        for i in xrange(l):#each position
            len2aaUsage[l].append(aas)
    return len2aaUsage

def readCdr3LenDist(file):
    '''Get length distribution of CDR3 amino-acid sequences. Format of input file: <length>\t<count>
    Return a dictionary len2freq, where key=length and val = freq of that length
    '''
    len2freq = {} #key = length, val=frequency
    f = open(file, 'r')
    total = 0
    for line in f:
        line = line.strip()
        if line == "" or line[0] == "#":
            continue
        #items = line.split("\t")
        items = line.split()
        if len(items) < 2:
            raise FileFormatError("Wrong cdr3LenDist file format. Require two fields: <length>\\t<count>\n.%s\n" %line)
        if items[0] == "TOTAL":
            continue
        try:
            l = int( (items[0].split("."))[0] )
            c = int( items[1] )
        except ValueError:
            raise ValueError("First and Second columns in cdr3LenDist file %s must be integer\n" %file)

        total += c
        if l not in len2freq:
            len2freq[l] = c
        else:
            len2freq[l] += c
    #Normalize counts to frequencies:
    if total > 0: 
        for l, c in len2freq.iteritems():
            len2freq[l] = float(c)/total
    f.close()
    return len2freq

def getUniCdr3LenDist():
    len2freq = {}
    minlen = 9
    maxlen = 23
    r = maxlen + 1 - minlen
    for l in xrange(minlen, maxlen + 1):
        len2freq[l] = 1.0/r
    return len2freq

def getTopFreqs( topFreqs, totalClones, totalSeqs ):
    f = open(topFreqs, "r")
    freqs = []
    for line in f:
        line = line.strip()
        freq = float(line)
        if not isinstance(freq, float):
            raise ValueError("Wrong format topFreqs file. Requires floats.\n")
        if freq <= 0:
            raise ValueError("Frequencies in topFreqs file (%s) must be positive. Got %f\n" %(topFreqs, freq))

        freqs.append( freq )
    #Make sure that sum of freqs <= 1 - (totalClones - len(freqs))/totalSeqs
    if sum(freqs) > 1 - float((totalClones - len(freqs)))/totalSeqs:
        raise TopFreqError("Frequencies of top clones add up to be larger than permitted")
    topFreqs = freqs
    f.close()
    return topFreqs
    
################### OPTIONS #########################
def checkOptions( parser, options, args ):
    #Check output directory
    if not os.path.exists( options.outdir ):
        system( "mkdir %s\n" %(options.outdir) )
    
    #Number of samples:
    if not isinstance(options.numSamples, int):
        parser.error("Number of samples must be in integer\n")
    if options.numSamples < 2:
        parser.error("Number of samples must > 2. %d was given\n" %options.numSamples)
    
    #Make sure the number of totalClones < totalSeqs:
    if options.totalClones > options.totalSeqs:
        raise ValueError("Total number of clones must be smaller than total number of sequences. %d clones and %d sequences were given" %(options.totalClones, options.totalSeqs) ) 

    #sampling sizes:
    options.samSize = [long(s) for s in options.samSize.split(",")]

    #Top Freqs:
    if options.topFreqs:
        if not os.path.exists( options.topFreqs ):
            parser.error("TopFreqs file %s does not exists\n" %options.topFreqs)
        options.topFreqs = getTopFreqs( options.topFreqs, options.totalClones, options.totalSeqs )

    #Cutoffs:
    options.cutoffs = [ float(c) for c in options.cutoffs.split(",") ]
     
    #CDR3 length distribution
    if options.cdr3LenDist:
        if not os.path.exists(options.cdr3LenDist):
            parser.error("cdr3LenDist file %s does not exist\n" %options.cdr3LenDist)
        options.len2freq = readCdr3LenDist(options.cdr3LenDist)
    else:
        options.len2freq = getUniCdr3LenDist() 

    if options.aaUsage:
        if not os.path.exists(options.aaUsage):
            parser.error("aaUsage file %s does not exist\n" %options.aaUsage)
        options.len2aaUsage = readAaUsage(options.aaUsage)
    else:
        options.len2aaUsage = getUniAaUsage( range(9, 24) )
    
    if not options.codonUsage:
        parser.error("Required argument codonUsage. None was given\n")
    if not os.path.exists(options.codonUsage):
        parser.error("codonUsage file %s does not exist\n" %options.codonUsage)
    options.aa2codons = readCodonUsage(options.codonUsage)

def addOptions( parser ):
    #parser.add_option("-i", "--indir", dest = 'indir', help="Required argument. Input directory that contains input fastq files (each represents a sample)\n")
    #parser.add_option("-n", "--popCdr3", dest='popCdr3', type="long", default=5*(10**11), help="Estimated number of unique CDR3 nucleotide sequences in the population. Default=5x10^11")
    parser.add_option("-o", "--outdir", dest = 'outdir', default=".", help="Output directory. Default: current directory\n")
    parser.add_option("-N", "--numSim", dest="numSim", type="int", default=100, help="Number of repertoire simulations. Default = 100")
    parser.add_option("--numSamplingPerSim", dest="numSamplingPerSim", type="int", default=100, help="Number of samplings per repertoire simulation. Default = %default")
    parser.add_option("-m", "--totalClones", dest="totalClones", type="long", default=3*(10**6), help="Estimated number of unique CDR3 nucleotide sequences in one individual. Default=3x10^6")
    parser.add_option("-M", "--totalSeqs", dest="totalSeqs", type="long", default=6*(10**7), help="Estimated number of CDR3 nucleotide sequences in one individual. Default=6x10^7 (number of CD8+ T cells per person)")
    parser.add_option("-S", "--samplingSize", dest="samSize", default="50000,100000,150000,200000", help="Comma separated string of sampling sizes. Default=%default")
    parser.add_option("-a", "--aaUsage", dest="aaUsage", help="Amino acid usage of the CDR3 sequences. Format:<aa-cdr3-length>\t<position>\t<letter>\t<count>. If not specified, use uniform distribution")
    parser.add_option("-l", "--cdr3LenDist", dest="cdr3LenDist", help="CDR3 length distribution. Format: <aa-cdr3-length>\t<count>")
    parser.add_option("-c", "--codonUsage", dest="codonUsage", default="~/codonUsage.txt", help="Codon usage table")
    parser.add_option("-f", "--topFreqs", dest="topFreqs", help="Frequencies of x top clones. Sum of these freqs must <= (M - (m-x))/M. If not specified, use uniform distribution")
    parser.add_option("-s", "--numSamples", dest="numSamples", type="int", default=2, help="Number of samples (individuals/repertoires) to compare. Default=2.")
    parser.add_option("-p", "--cutoffs", dest="cutoffs", default="0,0.001,0.01,0.05,0.1,0.5,1,5,10", help="Comma separated clonesize percentage cutoffs. Default='0,0.001,0.01,0.05,0.1,0.5,1,5'")
    parser.add_option("--maxSimPerRun", dest="maxSimPerRun", type='int', default=100, help="Maximum number of simulations per batch. Default=100")

################# MAIN ##############################
def main():
    parser = OptionParser()
    Stack.addJobTreeOptions( parser )

    addOptions( parser )

    options, args = parser.parse_args()
    checkOptions( parser, options, args )

    i = Stack( Setup(options, 0) ).startJobTree( options )
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" %i)

if __name__ == "__main__":
    #from tcrRepSim import *
    from immunoseq.src.tcrRepSim import *
    main()




