#!/usr/bin/env python

#nknguyen at soe ucsc edu
#Tue Jul 17 09:09:24 PDT 2012
#Immuno-seq pipeline

'''
Overlap stats 
'''

import os, sys, re, time, gzip, random
import random as rand
import cPickle as pickle
from optparse import OptionParser

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import logger
from sonLib.bioio import system

import immunoseqLib as iseqlib
from immunoseq.lib.overlapLib import *
from immunoseq.lib.overlap2Lib import *
from immunoseq.lib.lendistLib import *
import numpy as np
from scipy.stats import ttest_ind

###################################
############# OBJECTS #############
###################################
class PairOverlapStats():
    def __init__(self, reads1, reads2, clones1, clones2, oreads1, oreads2, oclones1, oclones2):
        self.reads1 = reads1
        self.reads2 = reads2
        self.oreads1 = oreads1
        self.oreads2 = oreads2
        self.oreadsAvr = float(oreads1 + oreads2)/2.0

        self.oreads1pc = iseqlib.getPc(oreads1, reads1)
        self.oreads2pc = iseqlib.getPc(oreads2, reads2)
        self.oreadsAvrPc = (self.oreads1pc + self.oreads2pc)/2.0

        self.clones1 = clones1
        self.clones2 = clones2
        self.oclones1 = oclones1
        self.oclones2 = oclones2
        self.oclonesAvr = float(oclones1 + oclones2)/2.0
        
        self.oclones1pc = iseqlib.getPc(oclones1, clones1)
        self.oclones2pc = iseqlib.getPc(oclones2, clones2)
        self.oclonesAvrPc = (self.oclones1pc + self.oclones2pc)/2.0



##############################################################################################
#========================= MAIN PIPELINE ====================================================#
##############################################################################################
class Setup(Target):
    '''Start of the main pipeline.
       Firstly, read input fasta files
       When done, call following steps to do sampling if needed and all the analyses
    '''
    def __init__(self, options):
        Target.__init__(self)
        self.options = options
    
    def run(self):
        #read input fasta files:
        ext = 'fa'
        files = iseqlib.getfiles(self.options.indir, ext)
        globalTempDir = self.getGlobalTempDir()
        for file in files:
            filepath = os.path.join(self.options.indir, file)
            self.addChildTarget( ReadFasta(filepath, globalTempDir, self.options.minReadCount) )
        
        self.setFollowOnTarget( SamplingAndAnalyses(globalTempDir, self.options) )

class ReadFasta(Target):
    '''Read input fasta file, return an pickle file of the sample object
       Where sample.seqs = {header: Seq}
    '''
    def __init__(self, file, outdir, mincount):
        Target.__init__(self)
        self.file = file
        self.outdir = outdir
        self.mincount = mincount
    
    def run(self):
        name = os.path.basename(self.file).split('.')[0]
        sample = iseqlib.Sample(name)
        mode = 1 #seqs = {header: seq}
        seqs, total = iseqlib.readFile( self.file, self.mincount, mode)
        for s in seqs.values():
            s.setFreq(total)
        sample.seqs = seqs
        sample.setTotal(total)

        picklefile = os.path.join(self.outdir, "%s.pickle" %name)
        pickle.dump( sample, gzip.open(picklefile, "wb") )

class SamplingAndAnalyses(Target):
    '''Sampling if needed and do the analyses. 
    '''
    def __init__(self, indir, options):
        Target.__init__(self)
        self.indir = indir
        self.options = options
        
    def run(self):
        globalTempDir = self.getGlobalTempDir()
        #Before sampling, chose 'numcombs' of sets of samples:
        group2samplesList = [] #list of group2subsetofsamples
        if self.options.numsamples:
            for i in xrange(self.options.numcombs):
                pickedGroup2samples = pickSamples(self.options.group2samples, self.options.numsamples)
                group2samplesList.append( pickedGroup2samples )
        else:
            group2samplesList.append(self.options.group2samples)
        self.options.g2slist = group2samplesList
        
        #print sampleSetId, group2subsetofsamples of that set:
        outfile = os.path.join(self.options.outdir, 'sampleSets.txt')
        f = open(outfile, 'w')
        for i, g2s in enumerate(group2samplesList):
            f.write("#Set %d\n" %i)
            for g, subsamples in g2s.iteritems():
                f.write('%s\t%s\n' %(g, ','.join(subsamples)))
            f.write("#\n")
        f.close()

        if self.options.sampling: 
            for i in xrange(self.options.numsam): #sampling a number of times
                samplingdir = os.path.join(globalTempDir, "%d" %i) #temp/i/
                self.addChildTarget( SamplingSamples(self.indir, samplingdir, self.options) ) #seqdir, outdir, options
        else:
            tempoutdir = os.path.join(globalTempDir, "0")
            system("mkdir -p %s" %tempoutdir)
            filterSamples(self.indir, self.options.vs, self.options.js)
            self.addChildTarget( Analyses(self.indir, tempoutdir, self.options) )

        #Calculate means & standard deviations of samplings
        self.setFollowOnTarget( AverageResults(globalTempDir, self.options) )

class SamplingSamples(Target):
    '''Sampling from all samples. When done, call followon jobs to calculate the stats needed
    '''
    def __init__(self, indir, statsOutdir, options):
        Target.__init__(self)
        self.indir = indir
        self.outdir = statsOutdir
        self.options = options

    def run(self):
        globalTempDir = self.getGlobalTempDir()
        ext = "pickle"
        files = iseqlib.getfiles(self.indir, ext)
        samples = [ '.'.join(file.split('.')[:-1]) for file in files ]
        
        for sample in samples:
            samplefile = os.path.join(self.indir, "%s.%s" %(sample, ext))
            outfile = os.path.join(globalTempDir, "%s.pickle" %sample) #temp/sample.pickle
            self.addChildTarget( Sampling(samplefile, outfile, self.options) )

        self.setFollowOnTarget( Analyses(globalTempDir, self.outdir, self.options) )

class Sampling(Target):
    '''Sampling and and call analyses
    '''
    def __init__(self, samplefile, outfile, options):
        Target.__init__(self)
        self.samplefile = samplefile
        self.outfile = outfile
        self.options = options
        self.size = options.sampling

    def run(self):
        sample = pickle.load( gzip.open(self.samplefile, "rb") )
        if self.options.uniq:
            subsample = iseqlib.samplingSample_weightedUniq(sample, self.size)
            #subsample = iseqlib.samplingSample_uniq(sample, self.size)
        else:
            subsample = iseqlib.samplingSample(sample, self.size)
        
        #filtering if selected Vs and/or selected Js were specified
        subsample = iseqlib.filterSampleByGenes(subsample, self.options.vs, self.options.js)
        
        pickle.dump(subsample, gzip.open(self.outfile, "wb"))

class  Analyses(Target):
    '''Call jobs to compute different analyses
    '''
    def __init__(self, indir, outdir, options):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir #temp/sampling#
        self.options = options

    def run(self):
        for i, g2s in enumerate(self.options.g2slist):
            outdir = os.path.join(self.outdir, "%d" %i) #outdir/sampling#/set#
            system('mkdir -p %s' %outdir)
            self.addChildTarget( RunSampleSetAnalyses(self.indir, outdir, g2s, self.options) )

####################################
############# ANALYSES #############
####################################
class RunSampleSetAnalyses(Target):
    '''
    '''
    def __init__(self, indir, outdir, group2samples, options):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.g2s = group2samples
        self.options = options

    def run(self):
        #Pairwise overlap:
        if 'pairwiseOverlap' in self.options.analyses:
            pairOutfile = os.path.join(self.outdir, "pairOverlap.pickle")
            self.addChildTarget( PairwiseOverlap(self.indir, pairOutfile, self.g2s, self.options) )
        
        #number of samples versus percentage of clones
        if 'numsamVsClones' in self.options.analyses:
            numSamOutfile = os.path.join(self.outdir, "numSamVsClones.pickle")
            lendistdir = os.path.join(self.outdir, 'lendist')
            system('mkdir -p %s' %lendistdir)
            self.addChildTarget( NumSamplesVsClones(self.indir, numSamOutfile, lendistdir, self.g2s, self.options) )
        
class PairwiseOverlap(Target):
    '''
    '''
    def __init__(self, indir, outfile, group2samples, options):
        Target.__init__(self)
        self.indir = indir
        self.outfile = outfile
        self.g2s = group2samples
        self.options = options

    def run(self):
        #Load the samples:
        g2s = {} #key = group, val = list of Sample objects
        for g, names in self.g2s.iteritems():
            g2s[g] = []
            for s in names:
                picklefile = os.path.join(self.indir, "%s.pickle" %s)
                sample = pickle.load( gzip.open(picklefile, "rb") )
                g2s[g].append(sample)

        g2stats = {} #key = group, val = {pair: stats}
        cutoffs = [0.0]
        mode = 2
        discrete = False
        #Pairs of intra-group samples
        for g, samples in g2s.iteritems():
            g2stats[g] = {}
            sample2aa2v2j = iseqlib.getsample2aa2v2j(samples)
            #all posible pairs
            for i in xrange(0, len(samples) -1):
                s1 = samples[i]
                aa2v2j1 = sample2aa2v2j[ s1.name ] 
                for j in xrange(i+1, len(samples)):
                    s2 = samples[j]
                    pair = '-'.join( sorted([s1.name, s2.name]) ) #pair = name1,name2
                    aa2v2j2 = sample2aa2v2j[ s2.name ]
                    reads1, reads2, clones1, clones2, stats1, stats2 = getPairwiseOverlap(s1.seqs, s2.seqs, aa2v2j1, aa2v2j2, cutoffs, mode, discrete)
                    stats = PairOverlapStats( reads1[0], reads2[0], clones1[0], clones2[0], stats1['oreads'][0], stats2['oreads'][0], stats1['oclones'][0], stats2['oclones'][0] )
                    g2stats[g][pair] = stats
        
        #Pair of inter-group samples:
        if self.options.crossGroup:
            groups = g2s.keys()
            for i1 in xrange( len(groups) -1 ):
                g1 = groups[i1]
                samples1 = g2s[g1]
                sample2aa2v2j1 = iseqlib.getsample2aa2v2j(samples1)
                for i2 in xrange( i1+1, len(groups) ):
                    g2 = groups[i2]
                    samples2 = g2s[g2]
                    sample2aa2v2j2 = iseqlib.getsample2aa2v2j(samples2)

                    g = '-'.join( sorted([g1, g2]) ) 
                    g2stats[g] = {}
                    for s1 in samples1:
                        aa2v2j1 = sample2aa2v2j1[s1.name]
                        for s2 in samples2:
                            aa2v2j2 = sample2aa2v2j2[s2.name]
                            pair = '-'.join( sorted([s1.name, s2.name]) ) #pair = name1-name2
                            reads1, reads2, clones1, clones2, stats1, stats2 = getPairwiseOverlap(s1.seqs, s2.seqs, aa2v2j1, aa2v2j2, cutoffs, mode, discrete)
                            stats = PairOverlapStats( reads1[0], reads2[0], clones1[0], clones2[0], stats1['oreads'][0], stats2['oreads'][0], stats1['oclones'][0], stats2['oclones'][0] )
                            g2stats[g][pair] = stats
        
        #pickle group2stats
        pickle.dump(g2stats, gzip.open(self.outfile, "wb"))

class NumSamplesVsClones(Target):
    '''
    '''
    def __init__(self, indir, outfile, lendistdir, group2samples, options):
        Target.__init__(self)
        self.indir = indir
        self.outfile = outfile
        self.lendistdir = lendistdir
        self.g2s = group2samples
        self.options = options

    def run(self):
        #Load the samples:
        g2s = {} #key = group, val = list of Sample objects
        for g, names in self.g2s.iteritems():
            g2s[g] = []
            for s in names:
                picklefile = os.path.join(self.indir, "%s.pickle" %s)
                sample = pickle.load( gzip.open(picklefile, "rb") )
                g2s[g].append(sample)
        
        #For each group, calculate the clone vs number of sample shared distribution
        g2stats = {}
        groupsamples = []
        for g, samples in g2s.iteritems():
            groupsample = combineSamples(samples, g)
            groupsamples.append(groupsample)
        uniq = True
        group2numsam2uniq = getSharedSeqDist(groupsamples, uniq)
        group2numsam2count = getSharedSeqDist(groupsamples, not uniq)
        stats = {'uniq': group2numsam2uniq, 'reads': group2numsam2count}
        pickle.dump(stats, gzip.open(self.outfile, 'wb'))

        #Length distribution for sequences that are shared by 1, 2, 3, ... samples
        g2numsam2lenfreq = getGroup2numsam2len2freq(groupsamples)
        for g, numsam2lenfreq in g2numsam2lenfreq.iteritems():
            for numsam, lenfreq in numsam2lenfreq.iteritems():
                lendistdir = os.path.join(self.lendistdir, '%d' %numsam)
                system("mkdir -p %s" % lendistdir)
                picklefile = os.path.join(lendistdir, "%s.pickle" %g)
                pickle.dump(lenfreq, gzip.open(picklefile, 'wb'))

##########################################################
####### COMPUTE MEANS & STDS OF THE SAMPLINGS ############
##########################################################
class AverageResults(Target):
    '''
    Indir/
        Sampling#/
            SampleSet#/
                pairOverlap.pickle
                numSamVsClones.pickle
                lendist/
                    numsam1/
                        group1.pickle
                        group2.pickle
    '''
    def __init__(self, indir, options):
        Target.__init__(self)
        self.indir = indir
        self.options = options

    def run(self):
        if 'pairwiseOverlap' in self.options.analyses:
            self.addChildTarget( PairwiseOverlapSummary(self.indir, self.options) )
        if 'numsamVsClones' in self.options.analyses:
            self.addChildTarget( NumSamplesVsClonesSummary(self.indir, self.options) )
            self.addChildTarget( NumSamplesLendistSummary(self.indir, self.options) )

########################################
########### SUMMARY ####################
########################################
class PairwiseOverlapSummary(Target):
    '''
    Indir/
        Sampling#/
            SampleSet#/
                pairOverlap.pickle
    '''
    def __init__(self, indir, options):
        Target.__init__(self)
        self.indir = indir
        self.options = options

    def run(self):
        #uniq = True
        outdir = os.path.join(self.options.outdir, "pairwiseOverlap")
        system("mkdir -p %s" % outdir)
        
        aggSet2group2pair2stats = {}
        #samplings = os.listdir(self.indir)
        samplings = finddirs(self.indir)
        #Accumulate samplings
        for i, sampling in enumerate(samplings):
            samplingdir = os.path.join(self.indir, sampling)
            #samplesets = os.listdir(samplingdir)
            samplesets = finddirs(samplingdir)
            for samset in samplesets:
                samsetdir = os.path.join(samplingdir, samset)
                picklefile = os.path.join(samsetdir, 'pairOverlap.pickle')
                group2pair2stats = pickle.load( gzip.open(picklefile, 'rb') )
                
                if samset not in aggSet2group2pair2stats:
                    aggSet2group2pair2stats[samset] = {}
                
                for group, pair2stats in group2pair2stats.iteritems():
                    if group not in aggSet2group2pair2stats[samset]:
                        aggSet2group2pair2stats[samset][group] = {}
                    
                    for pair, stats in pair2stats.iteritems():
                        count = stats.oclonesAvr
                        if pair not in aggSet2group2pair2stats[samset][group]:
                            aggSet2group2pair2stats[samset][group][pair] = [0]*i + [count]
                        else:
                            aggSet2group2pair2stats[samset][group][pair].append( count )

        #Compute means & stds for the samplings
        avrSet2stats = {}
        stdSet2stats = {}
        for samset, g2p2v in aggSet2group2pair2stats.iteritems():
            avrSet2stats[samset] = {}
            stdSet2stats[samset] = {}
            for g, p2v in g2p2v.iteritems():
                avrSet2stats[samset][g] = {}
                stdSet2stats[samset][g] = {}
                for p, vector in p2v.iteritems():
                    avrSet2stats[samset][g][p] = np.mean(vector)
                    stdSet2stats[samset][g][p] = np.std(vector)
        
        #Print to output files:
        for samset, stats in avrSet2stats.iteritems():
            setoutdir = os.path.join(outdir, samset)
            system("mkdir -p %s" %setoutdir)
            avroutfile = os.path.join(setoutdir, 'pairOverlapAvr.txt')
            printPairOverlapTable(avroutfile, stats)
            stdoutfile = os.path.join(setoutdir, 'pairOverlapStd.txt')
            printPairOverlapTable(stdoutfile, stdSet2stats[samset])
            
            #Draw:
            plotfile = os.path.join(setoutdir, 'pairOverlap')
            g2avr = {}
            g2std = {}
            g2p = {}
            stdstats = stdSet2stats[samset]
            for g, p2c in stats.iteritems():
                g2avr[g] = []
                g2std[g] = []
                g2p[g] = []
                for p, c in p2c.iteritems():
                    g2avr[g].append(c)
                    g2p[g].append(p)
                    g2std[g].append(stdstats[g][p])

            #drawPairOverlap(plotfile, g2avr, g2std, g2p)
            #drawPairOverlap_hack(plotfile, g2avr, g2std, g2p)

            #HACK ... Instead of drawing, just pickle the data:
            pickleFile = os.path.join(setoutdir, "data.pickle")
            pickle.dump( (g2avr, g2std, g2p), gzip.open(pickleFile, "wb") )

def printPairOverlapTable(outfile, group2pair2count):
    f = open(outfile, 'w')
    f.write("Group\tPairwiseStats\tMean\tStd\tPairs\n")
    for group, pair2count in group2pair2count.iteritems():
        pairs = pair2count.keys()
        counts = [pair2count[p] for p in pairs]
        f.write("%s\t%s\t%f\t%f\t%s\n" %(group, ",".join([str(c) for c in counts]), np.mean(counts), np.std(counts), ",".join(pairs)))
    
    #ttest:
    f.write("\n#ttest\n")
    f.write("Group1-group2\tPvalue\n")
    groups = group2pair2count.keys()
    for i1 in xrange( len(groups) - 1 ):
        g1 = groups[i1]
        counts1 = group2pair2count[g1].values()
        for i2 in xrange( i1+1, len(groups) ):
            g2 = groups[i2]
            counts2 = group2pair2count[g2].values()
            tval, pval = ttest_ind(counts1, counts2)
            f.write("%s-%s\t%f\n" %(g1, g2, pval))
    f.close()

class NumSamplesVsClonesSummary(Target):
    '''
    Indir/
        Sampling#/
            SampleSet#/
                numSamVsClones.pickle
    '''
    def __init__(self, indir, options):
        Target.__init__(self)
        self.indir = indir
        self.options = options

    def run(self):
        outdir = os.path.join(self.options.outdir, "numsamVsClones")
        system("mkdir -p %s" % outdir)
        
        aggSet2group2numsam2uniq = {}
        aggSet2group2numsam2count = {}
        samplings = finddirs(self.indir)
        #Accumulate samplings
        for i, sampling in enumerate(samplings):
            samplingdir = os.path.join(self.indir, sampling)
            samplesets = finddirs(samplingdir)
            for samset in samplesets:
                samsetdir = os.path.join(samplingdir, samset)
                picklefile = os.path.join(samsetdir, 'numSamVsClones.pickle')
                #stats = {'uniq': group2numsam2uniq, 'reads': group2numsam2count}
                stats = pickle.load( gzip.open(picklefile, 'rb') )
                ustats = stats['uniq']
                rstats = stats['reads']
                
                if samset not in aggSet2group2numsam2uniq:
                    aggSet2group2numsam2uniq[samset] = {}
                    aggSet2group2numsam2count[samset] = {}
                for group, numsam2uniq in ustats.iteritems():
                    if group not in aggSet2group2numsam2uniq[samset]:
                        aggSet2group2numsam2uniq[samset][group] = {}
                        aggSet2group2numsam2count[samset][group] = {}
                    for numsam, uniq in numsam2uniq.iteritems():
                        count = rstats[group][numsam]
                        if numsam not in aggSet2group2numsam2uniq[samset][group]:
                            aggSet2group2numsam2uniq[samset][group][numsam] = [0]*i + [uniq]
                            aggSet2group2numsam2count[samset][group][numsam] = [0]*i + [count]
                        else:
                            aggSet2group2numsam2uniq[samset][group][numsam].append(uniq)
                            aggSet2group2numsam2count[samset][group][numsam].append(count)
        #Means & Stds
        avrSet2statsU = {}
        stdSet2statsU = {}
        avrSet2statsR = {}
        stdSet2statsR = {}
        for samset, g2n2v in aggSet2group2numsam2uniq.iteritems():
            avrSet2statsU[samset] = {}
            stdSet2statsU[samset] = {}
            avrSet2statsR[samset] = {}
            stdSet2statsR[samset] = {}
            for g, n2v in g2n2v.iteritems():
                avrSet2statsU[samset][g] = {}
                stdSet2statsU[samset][g] = {}
                avrSet2statsR[samset][g] = {}
                stdSet2statsR[samset][g] = {}
                for n, vector in n2v.iteritems():
                    avrSet2statsU[samset][g][n] = np.mean(vector)
                    stdSet2statsU[samset][g][n] = np.std(vector)
                    readvector = aggSet2group2numsam2count[samset][g][n]
                    avrSet2statsR[samset][g][n] = np.mean(readvector)
                    stdSet2statsR[samset][g][n] = np.std(readvector)
        
        #Print to output files:
        for samset, g2dist in avrSet2statsU.iteritems():
            setoutdir = os.path.join(outdir, samset)
            system("mkdir -p %s" %setoutdir)
            #uniq
            uniqoutfile = os.path.join(setoutdir, 'numsamVsClones')
            drawNumsamVsClonesDist2(g2dist, uniqoutfile)
            #count
            readoutfile = os.path.join(setoutdir, 'numsamVsReads')
            drawNumsamVsClonesDist2(avrSet2statsR[samset], readoutfile)

            #table:
            uniqtabfile = os.path.join(setoutdir, 'numsamVsClones.txt')
            printNumsamVsClonesDist(uniqtabfile, g2dist, False)
            uniqtabfile = os.path.join(setoutdir, 'numsamVsClones_rel.txt')
            printNumsamVsClonesDist(uniqtabfile, g2dist, True)
            readtabfile = os.path.join(setoutdir, 'numsamVsReads.txt')
            printNumsamVsClonesDist(readtabfile, avrSet2statsR[samset], False)
            readtabfile = os.path.join(setoutdir, 'numsamVsReads_rel.txt')
            printNumsamVsClonesDist(readtabfile, avrSet2statsR[samset], True)


##########################################
########## CDR3 LENGTH DISTRIBUTION ######
##########################################
class NumSamplesLendistSummary(Target):
    '''Length distribution
    Indir/
        Sampling#/
            SampleSet#/
                lendist/
                    numsam1/
                        group1.pickle
                        group2.pickle
    '''
    def __init__(self, indir, options):
        Target.__init__(self)
        self.indir = indir
        self.options = options

    def run(self):
        samplings = finddirs(self.indir)
        
        samset2intersectNumsam = {} 
        for sampling in samplings: #each sampling
            samplingdir = os.path.join(self.indir, sampling)
            samsets = finddirs(samplingdir)
            for samset in samsets:
                lendistdir = os.path.join(samplingdir, samset, "lendist")
                numsams = finddirs(lendistdir)
                numsams = [int(n) for n in numsams]
                maxnumsam = max(numsams)
                if samset not in samset2intersectNumsam or samset2intersectNumsam[samset] > maxnumsam:
                    samset2intersectNumsam[samset] = maxnumsam
        
        for samset, numsam in samset2intersectNumsam.iteritems():
            for i in xrange(1, numsam + 1):
                self.addChildTarget( Lendist(self.indir, samset, i, self.options) )

class Lendist(Target):
    def __init__(self, indir, samset, numsam, options):
        Target.__init__(self)
        self.indir = indir
        self.samset = samset
        self.numsam = numsam
        self.options = options
    
    def run(self):
        group2avr = {}
        group2std = {}

        group2agglenfreq = {}
        samplings = finddirs(self.indir)
        for i, sampling in enumerate(samplings):
            samplingdir = os.path.join(self.indir, sampling)
            indir = os.path.join(samplingdir, self.samset, "lendist", str(self.numsam))
            #get stats for all samplings
            groups = os.listdir(indir)
            groups = [g.split('.')[0] for g in groups]
            for g in groups:
                picklefile = os.path.join(indir, "%s.pickle" %g)
                len2freq = pickle.load( gzip.open(picklefile, 'rb') )
                #update aggType2gene2count
                if g not in group2agglenfreq:
                    group2agglenfreq[g] = {}
                addSamplingLenDist(len2freq, group2agglenfreq[g], i)
            
        #Average stats of the samplings:
        for g, agglenfreq in group2agglenfreq.iteritems():
            l2avr, l2std = avrSamplingLenDist(agglenfreq)
            group2avr[g] = l2avr
            group2std[g] = l2std
        #Get union lendist:
        getUnionLenDist(group2avr)
        getUnionLenDist(group2std)

        #Draw lendist:
        outdir = os.path.join(self.options.outdir, "lendist", self.samset, str(self.numsam))
        system("mkdir -p %s" %outdir)
        uniq = True
        uniqoutfile = os.path.join(outdir, 'lendist_uniq')
        drawLenDist2(group2avr, uniqoutfile, uniq)
        readoutfile = os.path.join(outdir, 'lendist_read')
        drawLenDist2(group2avr, readoutfile, not uniq)

        #Print tables:
        taboutdir = os.path.join(outdir, 'tables')
        system('mkdir -p %s' %taboutdir)
        printLenDistTab(group2avr, group2std, taboutdir)

        #t-test
        #if self.options.group2samples:
        #    ttestoutfile = os.path.join(outdir, 'ttests.txt')
        #    lenDistTtests(sample2avr, self.options.group2samples, ttestoutfile, self.options.ttestTargetGroup)

        #clonesize distribution for each length

###################### UTILITIES FUNCTIONS #########################
def finddirs(indir):
    files = os.listdir(indir)
    dirs = []
    for file in files:
        if os.path.isdir( os.path.join(indir, file) ):
            dirs.append(file)
    return dirs

def pickSamples(group2samples, samPerGroup):
    #randomly pick 'samPerGroup' number of samples for each group
    g2s = {}
    for g, samples in group2samples.iteritems():
        if len(samples) < samPerGroup:
            raise ValueError("Not enough number of sample for group %s. Asked for %d, only has %d\n" %(g, samPerGroup, len(samples)))
        subsamples = rand.sample( samples, samPerGroup )
        g2s[g] = subsamples
    return g2s

def getGroup2numsam2len2freq(groupsamples):
    g2stats = {} # g2stats = {groupname: {numSampleShared: {len2freq}}}, where len2freq={len: [%reads,%clones]}
    for sample in groupsamples:#each group
        numsam2len2freq = {}
        for header, name2seq in sample.seqs.iteritems():
            numsam = len( name2seq.keys() )
            seq0 = name2seq.values()[0]
            seqlen = len(seq0.seq) 
            totalreads = sum([seq.count for seq in name2seq.values()])
            if numsam not in numsam2len2freq:
                numsam2len2freq[numsam] = {seqlen: [totalreads, 1]}
            else:
                if seqlen not in numsam2len2freq[numsam]:
                    numsam2len2freq[numsam][seqlen] = [totalreads, 1]
                else:
                    counts = numsam2len2freq[numsam][seqlen]
                    numsam2len2freq[numsam][seqlen] = [counts[0] + totalreads, counts[1] + 1]
        
        #Calculate frequencies:
        for numsam, len2freq in numsam2len2freq.iteritems():
            totalUniqs = sum([counts[1] for counts in len2freq.values()])
            totalReads = sum([counts[0] for counts in len2freq.values()])
            if totalUniqs > 0 and totalReads > 0:
                for l, counts in len2freq.iteritems():
                    len2freq[l] = [ counts[0]*100.0/totalReads, counts[1]*100.0/totalUniqs ]
        g2stats[sample.name] = numsam2len2freq
    return g2stats

def filterSamples(indir, vs, js):
    if not vs and not js:
        return
    ext = 'pickle'
    files = iseqlib.getfiles(indir, ext)
    for file in files:
        filepath = os.path.join(indir, file)
        sample = pickle.load( gzip.open(filepath, "rb") )
        subsample = iseqlib.filterSampleByGenes(sample, vs, js)
        system("rm %s" %filepath)
        pickle.dump(subsample, gzip.open(filepath, "wb"))
    return

############## ERROR CLASSES ######################
class InputOptionError(Exception):
    pass

def checkOptions( parser, options, args ):
    if options.group2samplesFile:
        options.group2samples, options.sample2group = iseqlib.readGroup2samples(options.group2samplesFile)
    if options.vs:
        options.vs = options.vs.split(',')
    if options.js:
        options.js = options.js.split(',')
    analyses = ["pairwiseOverlap", "numsamVsClones"]
    if options.analyses == 'all':
        options.analyses = analyses
    else:
        options.analyses = options.analyses.split(',')
        for a in options.analyses:
            if a not in analyses:
                raise ValueError("Invalid analysis choice %s. Valid choices are: pairwiseOverlap, numsamVsClones")

def addOptions( parser ):
    parser.add_option("-i", "--indir", dest = 'indir', help="Required argument. Input directory that contains fasta files of productive CDR3 sequences (each file represents a sample)\n")
    parser.add_option("-o", "--outdir", dest = 'outdir', default='immunoseqOut', help="Required argument. Output directory. Default=%default\n")
    parser.add_option('--sampling', dest='sampling', type='int', help='Sampling size if wish to normalize the read counts of all samples (randomly select the specified size of reads from each sample). Default is no sampling.')
    parser.add_option("--numsam", dest='numsam', type='int', default=1, help='If sampling is specified, this is the number of samplings to perform for each sample')
    parser.add_option("--group2samples", dest= 'group2samplesFile', help='File specifying which samples belong to which group. Format: <group> <comma-separted-list-of-samples>')
    parser.add_option("--minReadCount", dest='minReadCount', type='int', default=1, help='Minimum read count: clones with < this number will be filtered out. Default=%default')
    parser.add_option('--crossGroup', dest='crossGroup', action='store_true', default=False, help='If specified, include group "controlVersusPatient", which includes pairs of patient-control')
    parser.add_option('--ttestTargetGroup', dest='ttestTargetGroup', help='If specified, will compare samples belong to this group versus all other samples')
    parser.add_option('--numsamples', dest='numsamples', type='int', help='If specified, randomly select numsamples of samples from each group')
    parser.add_option('--numcombinations', dest='numcombs', type='int', default = 1, help='If specified, randomly select numcomb set of numsamples of samples from each group')
    #parser.add_option('--avrcombinations', dest='avrcombinations', action='store_true', default=False, help='If specified, calculate average stats across different combinations of samples')
    parser.add_option('--vs', dest='vs', help='Comma separated list of selected V genes. If specified, only calculate stats for clones with at least 1 of these genes, filter out all other clones')
    parser.add_option('--js', dest='js', help='Comma separated list of selected J genes. If specified, only calculate stats for clones with at least 1 of these genes, filter out all other clones')
    parser.add_option('--analyses', dest='analyses', default='all', help='Comma separated list of analyses to be performed. Valid values are ["pairwiseOverlap", "numsamVsClones", "lendist"].Default=%default.')
    parser.add_option('-u', '--uniq', dest='uniq', action='store_true', default=False, help='If specified, will perform the weighted-unique (clones) samplings insteads of reads samplings')
    #parser.add_option("--minReadPercentage", dest='minReadPercentage', default='0,0.001,0.01,0.05,0.1,0.5,1,5', help='Comma separated list of read percentage cutoffs. Clone must have at least this percentage of total reads to be included. Default=%default')

def main():
    parser = iseqlib.initOptions()
    Stack.addJobTreeOptions( parser )

    addOptions( parser )

    options, args = parser.parse_args()
    checkOptions( parser, options, args )

    i = Stack( Setup(options) ).startJobTree( options )
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" %i)

if __name__ == "__main__":
    from immunoseq.src.overlapWithSampling import *
    main()

