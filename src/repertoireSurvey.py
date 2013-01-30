#!/usr/bin/env python2.6

#Mon Jul  2 10:17:47 PDT 2012
#nknguyen soe ucsc edu
#
#Repertoire survey

'''
Report different statistics (#uniq clones, diversity index (both for species richness and abundance)) of input repertoires, 
as well as pairwise similarity indices of them, through an increasing range of sampling space.

Including:
    0/ diversity indices summary
    1/ rarefraction analyses
    2/ similarity indices
    
'''

import os, sys, re, time, copy, gzip
import random as rnd
import cPickle as pickle
from optparse import OptionParser

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import logger
from sonLib.bioio import system
import numpy as np
#import rpy2.robjects as robjs

import immunoseqLib as iseqlib

###################### OBJECTS ##############################
class SingleSamplingStats:
    '''Different stats of one sampling from a specific sample with a specific size 
    '''
    def __init__(self):
        self.uniqClones = None
        #R vegan: diversity
        self.simpson = None 
        self.invsimpson = None
        self.shannon = None
        self.fisherAlpha = None
        #self.rarefy = None

        #Standard deviation
        self.uniqClonesStd = None
        self.simpsonStd = None 
        self.invsimpsonStd = None
        self.shannonStd = None
        self.fisherAlphaStd = None

    def __getitem__(self, name):
        #name2val = {'uniqClones': self.uniqClones, 'uniqClonesStd': self.uniqCloneStd,
        #            'simpson': self.simpson, 'simpsonStd': self.simpsonStd,
        #            'invsimpson': self.invsimpson, 'invsimpsonStd': self.invsimpsonStd,
        #            'shannon': self.shannon, 'shannonStd': self.shannonStd,
        #            'fisherAlpha': self.fisherAlpha, 'fisherAlphaStd': self.fisherAlphaStd}
        if name not in self.__dict__:
            raise KeyError("SingleSamplingStats does not have attribute %s\n" %name)
        return self.__dict__[name]

    def __setitem__(self, name, val):
        self.__dict__[name] = val

class PairSamplingStats:
    '''
    '''
    def __init__(self):
        self.horn = None
        self.chao = None
        self.bray = None
        self.morisita = None
        self.manhattan = None
        self.euclidean = None
        self.canberra = None
        self.kulczynski = None
        self.jaccard = None
        self.gower = None
        self.mountford = None
        self.raup = None
        self.binomial = None
        #"manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn", "mountford", "raup" , "binomial" or "chao"
        self.hornStd = None
        self.chaoStd = None
        self.brayStd = None
        self.morisitaStd = None
        self.manhattanStd = None
        self.euclideanStd = None
        self.canberraStd = None
        self.kulczynskiStd = None
        self.jaccardStd = None
        self.gowerStd = None
        self.mountfordStd = None
        self.raupStd = None
        self.binomialStd = None
    
    def __getitem__(self, name):
        if name not in self.__dict__:
            raise KeyError("SingleSamplingStats does not have attribute %s\n" %name)
        return self.__dict__[name]

    def __setitem__(self, name, val):
        self.__dict__[name] = val

###################### MAIN PIPELINE ###############################
class Setup(Target):
    '''Setting up the pipeline
    '''
    def __init__(self, options):
        Target.__init__(self)
        self.options = options
        
    def run(self):
        ext = 'fa'
        samples = iseqlib.getfiles(self.options.indir, ext)
        globalTempDir = self.getGlobalTempDir()

        #Read input fasta files and write pickle files into globalTempDir:
        for sample in samples:
            name = sample.rstrip(ext).rstrip('.')
            infile = os.path.join(self.options.indir, sample)
            self.addChildTarget( ReadFasta(infile, globalTempDir, name) )
        
        #After done reading fastas, move to the analyses
        self.setFollowOnTarget( Analyses(globalTempDir, self.options) )

class Analyses(Target):
    '''Set up jobs to do two main analyses: stats for single samples, and similarity measurement for pair of samples
    '''
    def __init__(self, samdir, options):
        Target.__init__(self)
        self.samdir = samdir #directory that contains pickled Sample obj
        self.options = options

    def run(self):
        if self.options.diversity:
            self.addChildTarget( SingleAnalyses(self.samdir, self.options) ) #Stats for single samples
        if self.options.similarity:
            self.addChildTarget( PairAnalyses(self.samdir, self.options) ) #Compare pair of samples
        
        #self.setFollowOnTarget( )

#===================== DIVERSITY METRICS FOR SINGLE SAMPLES =============================
class SingleAnalyses(Target):
    '''Set up jobs for analyses of single samples, including: rarefraction, diversity indices: shannon, simpson, fisher
    '''
    def __init__(self, samdir, options):
        Target.__init__(self)
        self.samdir = samdir
        self.options = options

    def run(self):
        singleOutdir = os.path.join(self.options.outdir, "diversity")
        system("mkdir -p %s" %singleOutdir)
        globalTempDir = self.getGlobalTempDir()
        ext = 'pickle'
        samples = iseqlib.getfiles(self.samdir, ext)
        for sample in samples: #Each sample
            samplename = sample.rstrip(ext).rstrip('.')
            self.addChildTarget( SampleSingleAnalyses(globalTempDir, samplename, self.samdir, self.options, singleOutdir) )
        # R --no-save --no-restore --args adapt16D-adapt11D.txt < diversityPlot.R 
        self.setFollowOnTarget( SummarySingle(globalTempDir, singleOutdir, self.options.diversityIndices) )

class SummarySingle(Target):
    '''For each statistic, print out a table of all samples (Rows = samples, columns = different sampling sizes
    '''
    def __init__(self, indir, outdir, metrics):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.metrics = metrics

    def run(self):
        ext = 'pickle'
        files = iseqlib.getfiles(self.indir, ext)
        sample2size2stats = {}
        sizes = []
        for file in files:
            name = file.rstrip(ext).rstrip('.')
            size2stats = pickle.load( gzip.open(os.path.join(self.indir, file), "rb") )
            for size in size2stats:
                if size not in sizes:
                    sizes.append(size)
            sample2size2stats[name] = size2stats
        sizes.sort()

        #Print summary of each statistic to output files (1 file/statistic where row = samples, columns = sampling size)
        metrics = self.metrics
        metricsStd = [m + "Std" for m in metrics]

        for i, metric in enumerate(metrics):
            outfile = os.path.join(self.outdir, "%s.txt" %metric)
            f = open(outfile, 'w')
            f.write("Sample")
            for size in sizes:
                f.write("\t%d\tStd" %size)
            for sample, size2stats in sample2size2stats.iteritems():
                f.write("\n%s" %sample)
                for size in sizes:
                    if size not in size2stats:
                        f.write("\tNA\tNA")
                    else:
                        s = size2stats[size]
                        f.write("\t%f\t%f" % (s[metric], s[metricsStd[i]]) )
            f.write("\n")
            f.close()

        #Summary of all statictics for each sampling size (each file per sampling size, row=samples, columns = statistics
        for size in sizes:
            outfile = os.path.join(self.outdir, "%d.txt" %size)
            f = open(outfile, 'w')
            f.write("Sample")
            for m in metrics:
                f.write("\t%s\tStd" %m)
            
            for sample, size2stats in sample2size2stats.iteritems():
                f.write("\n%s" %sample)
                for i, metric in enumerate(metrics):
                    if size not in size2stats:
                        f.write("\tNA\tNA")
                    else:
                        s = size2stats[size]
                        f.write("\t%f\t%f" %(s[metric], s[metricsStd[i]]) )
            f.write("\n")
            f.close()
        
class SampleSingleAnalyses(Target):
    '''Single analyses for a specific sample
    '''
    def __init__(self, tempOutdir, samplename, samdir, options, outdir):
        Target.__init__(self)
        self.tempOutdir = tempOutdir
        self.name = samplename
        self.samdir = samdir
        self.options = options
        self.outdir = outdir

    def run(self):
        globalTempDir = self.getGlobalTempDir() #will contain samsize.pickle with average stats of each sampling size
        picklefile = os.path.join(self.samdir, "%s.pickle" %self.name)
        
        #Get sampling sizes
        samSizes = getSamplingSizes(self.options, picklefile)
        for size in samSizes: #each sampling size
            self.addChildTarget( SamplingSampleSingle(globalTempDir, picklefile, size, self.options) )

        self.setFollowOnTarget( SummarySampleSingle(globalTempDir, self.tempOutdir, self.name, self.outdir, self.options.diversityIndices) )
        
class SummarySampleSingle(Target):
    '''Summary stats for the sample across different sampling sizes
    '''
    def __init__(self, indir, tempOutdir, name, outdir, metrics):
        Target.__init__(self)
        self.indir = indir
        self.tempOutdir = tempOutdir
        self.name = name
        self.outdir = outdir
        self.metrics = metrics
        
    def run(self):
        ext = 'pickle'
        files = iseqlib.getfiles(self.indir, ext)
        sizes = []
        size2stats = {} #key = size, val = averStats
        for file in files:
            size = int(file.rstrip(ext).rstrip('.'))
            sizes.append(size)
            stats = pickle.load( gzip.open(os.path.join(self.indir, file), "rb") )
            size2stats[size] = stats
        sizes.sort()

        #output summary file of the sample:
        outfile = os.path.join( self.outdir, "%s.txt" %self.name )
        f = open(outfile, 'w')
        f.write("Index")
        for s in sizes:
            f.write("\t%d\tStd" %s)
        metrics = self.metrics
        metricsStd = [m + "Std" for m in metrics]
        for i, metric in enumerate(metrics):
            f.write("\n%s" %metric)
            for size in sizes:
                avr = size2stats[size][metric]
                std = size2stats[size][metricsStd[i]]
                f.write("\t%f\t%f" %(avr, std))
        f.write("\n")
        f.close()

        #pickle size2stats to temporary output directory
        picklefile = os.path.join(self.tempOutdir, "%s.pickle" % self.name)
        pickle.dump( size2stats, gzip.open(picklefile, "wb") )

class SamplingSampleSingle(Target):
    '''Sampling the sample with a specific size for a specified number of times, report mean and std
    '''
    def __init__(self, outdir, sampleobjfile, size, options):
        Target.__init__(self)
        self.outdir = outdir
        self.sampleobjfile = sampleobjfile
        self.size = size
        self.options = options

    def run(self):
        globalTempDir = self.getGlobalTempDir()
        for i in xrange(0, self.options.numSamplings):
            self.addChildTarget( SamplingSingle(globalTempDir, self.sampleobjfile, self.size, i, self.options) )
        self.setFollowOnTarget( AverageSampling(globalTempDir, self.outdir, self.size, self.options.diversityIndices) )
        
class SamplingSingle(Target):
    '''One sampling of the sample with a specific size
    '''
    def __init__(self, outdir, sampleobjfile, size, id, options):
        Target.__init__(self)
        self.outdir = outdir
        self.sampleobjfile = sampleobjfile
        self.size = size
        self.id = id
        self.options = options

    def run(self):
        sample = pickle.load( gzip.open(self.sampleobjfile, "rb") )
        #sampling "size" sequences from sample:
        subsample = iseqlib.samplingSample( sample, self.size ) #sampled sample
        
        #Get different statistics:
        stats = SingleSamplingStats() #initialize stats for the current subsample

        #1/ Number of uniq clones:
        stats.uniqClones = len( subsample.seqs )

        #2/ Diversity Indices
        counts = [ seq.count for seq in subsample.seqs.values() ] #list of clone sizes

        ## Initialize the R interface
        import rpy2.rinterface as rinterface
        rinterface.set_initoptions(('rpy2', '--no-save'))
        rinterface.initr()

        import rpy2.robjects as robjs
        from rpy2.robjects.packages import importr
        vegan = importr('vegan')
        rcounts = robjs.IntVector( counts ) #convert python list into R vector
        
        ##Diversity indices (Simpson, InverseSimpson, Shannon)
        indices = self.options.diversityIndices
        #['simpson', 'invsimpson', 'shannon']
        for index in indices:
            if index == 'uniqClones' or index == 'fisherAlpha':
                continue
            rval = vegan.diversity( rcounts, index )
            stats[index] = rval[0]
            
        ## Fisher Alpha:
        rfisher = vegan.fisher_alpha( rcounts )
        stats.fisherAlpha = rfisher[0]

        #Write to temp file:
        picklefile = os.path.join( self.outdir, "%d.pickle" % self.id)
        pickle.dump( stats, gzip.open(picklefile, "wb") )

class AverageSampling(Target):
    '''Mean and standard deviation of current sampling size
    '''
    def __init__(self, indir, outdir, size, metrics):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.samplingsize = size
        self.metrics = metrics

    def run(self):
        ext = 'pickle'
        picklefiles = iseqlib.getfiles( self.indir, ext )
        statsList = []
        for file in picklefiles:
            stats = pickle.load( gzip.open(os.path.join(self.indir, file), "rb") )
            statsList.append(stats)
        
        avrstats = SingleSamplingStats() #initialize avrstats
        #Calculate mean and std using numpy
        #metrics = ['uniqClones', 'simpson', 'invsimpson', 'shannon', 'fisherAlpha']
        #stds = ['uniqClonesStd', 'simpsonStd', 'invsimpsonStd', 'shannonStd', 'fisherAlphaStd']
        metrics = self.metrics
        stds = [m + "Std" for m in metrics]

        for i in xrange( len(metrics) ):
            vals = [s[metrics[i]] for s in statsList]
            avrstats[ metrics[i] ] = np.mean( vals )
            avrstats[ stds[i] ] = np.std( vals )

        #Pickle the average stat of this sampling size 
        picklefile = os.path.join(self.outdir, "%d.pickle" %self.samplingsize)
        pickle.dump( avrstats, gzip.open(picklefile, "wb") )
#===================== END METRICS FOR SINGLE SAMPLES =============================

#===================== PAIRWISE COMPARISONS, SIMILARITY INDICES ===================
class PairAnalyses(Target):
    '''Pairwise comparison
    '''
    def __init__(self, samdir, options):
        Target.__init__(self)
        self.samdir = samdir
        self.options = options

    def run(self):
        pairOutdir = os.path.join(self.options.outdir, 'similarity')
        system("mkdir -p %s" %pairOutdir)

        globalTempDir = self.getGlobalTempDir()
        ext = 'pickle'
        samples = iseqlib.getfiles(self.samdir, ext)
        samplenames = [s.rstrip(ext).rstrip('.') for s in samples]
        for i in xrange( len(samples) - 1 ):
            s1 = samples[i]
            s1name = samplenames[i] 
            for j in xrange( i+1, len(samples) ):
                s2 = samples[j]
                s2name = samplenames[j]
                self.addChildTarget( SamplePairAnalyses(globalTempDir, s1name, s2name, self.samdir, self.options, pairOutdir) )
        
        self.setFollowOnTarget( SummaryPair(globalTempDir, pairOutdir, self.options.similarityIndices) )

class SummaryPair(Target):
    '''
    '''
    def __init__(self, indir, outdir, metrics):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.metrics = metrics

    def run(self):
        ext = 'pickle'
        files = iseqlib.getfiles(self.indir, ext)
        sample2mate2size2stats = {}
        sizes = []

        for file in files:
            name = file.rstrip(ext).rstrip('.')
            samples = name.split('-')
            size2stats = pickle.load( gzip.open(os.path.join(self.indir, file), "rb") )
            for size in size2stats:
                if size not in sizes:
                    sizes.append(size)

            s1 = samples[0]
            s2 = samples[1]
            if s1 not in sample2mate2size2stats:
                sample2mate2size2stats[s1] = { s2: size2stats }
            else:
                sample2mate2size2stats[s1][s2] = size2stats

            if s2 not in sample2mate2size2stats:
                sample2mate2size2stats[s2] = { s1: size2stats }
            else:
                sample2mate2size2stats[s2][s1] = size2stats

            #for i, sample in enumerate(samples):
            #    if sample not in sample2mate2size2stats:
            #        sample2mate2size2stats[sample] = {samples[(i+1) %2]: size2stats}
            #    else:
            #        sample2mate2size2stats[sample][samples[(i+1)%2]] = size2stats
        sizes.sort()

        #Print summary of each statistic to output files (1 file/1 statistic, 1 sampling size where row = samples, cols = samples)
        metrics = self.metrics
        metricsStd = [m + "Std" for m in metrics]

        samples = sorted(sample2mate2size2stats.keys())
        for i, metric in enumerate(metrics):
            for size in sizes:
                outfile = os.path.join(self.outdir, "%s-%d.txt" %(metric, size))
                f = open(outfile, 'w')
                f.write("Sample")
                for s in samples:
                    f.write("\t%s\tStd" %(s) )
                for s in samples:
                    f.write("\n%s" %s)
                    for s2 in samples:
                        if s == s2:
                            f.write("\t-\t-")
                        else:
                            if s in sample2mate2size2stats and s2 in sample2mate2size2stats[s]:
                                size2stats = sample2mate2size2stats[s][s2]
                                if size not in size2stats:
                                    f.write("\t-\t-")
                                else:
                                    stat= size2stats[size]
                                    f.write("\t%f\t%f" %(stat[metric], stat[ metricsStd[i] ] ))
                            else:
                                f.write("\t-\t-")
                f.write("\n")
                f.close()
            
class SamplePairAnalyses(Target):
    '''Comparing a specific pair of samples
    '''
    def __init__(self, tempOutdir, samplename1, samplename2, samdir, options, outdir):
        Target.__init__(self)
        self.tempOutdir = tempOutdir
        self.name1 = samplename1
        self.name2 = samplename2
        self.samdir = samdir
        self.options = options
        self.outdir = outdir

    def run(self):
        globalTempDir = self.getGlobalTempDir() #will contain samsize.pickle with average stats of each sampling size
        picklefile1 = os.path.join(self.samdir, "%s.pickle" %self.name1)
        picklefile2 = os.path.join(self.samdir, "%s.pickle" %self.name2)

        #Get sampling sizes
        samSizes1 = getSamplingSizes(self.options, picklefile1)
        samSizes2 = getSamplingSizes(self.options, picklefile2)
        samSizes = samSizes1
        if max(samSizes1) > max(samSizes2):
            samSizes = samSizes2

        for size in samSizes: #each sampling size
            self.addChildTarget( SamplingSamplePair(globalTempDir, picklefile1, picklefile2, size, self.options) )

        self.setFollowOnTarget( SummarySamplePair(globalTempDir, self.tempOutdir, self.name1, self.name2, self.outdir, self.options.similarityIndices) )

class SummarySamplePair(Target):
    '''Summary stats for the sample across different sampling sizes
    '''
    def __init__(self, indir, tempOutdir, name1, name2, outdir, metrics):
        Target.__init__(self)
        self.indir = indir
        self.tempOutdir = tempOutdir
        self.name1 = name1
        self.name2 = name2
        self.outdir = outdir
        self.metrics = metrics
        
    def run(self):
        ext = 'pickle'
        files = iseqlib.getfiles(self.indir, ext)
        sizes = []
        size2stats = {} #key = size, val = averStats
        print self.name1, self.name2
        for file in files:
            size = int(file.rstrip(ext).rstrip('.'))
            sizes.append(size)
            stats = pickle.load( gzip.open(os.path.join(self.indir, file), "rb") )
            size2stats[size] = stats
            print sizes
            print stats.horn, stats.chao
        sizes.sort()

        #output summary file of the sample:
        outfile = os.path.join( self.outdir, "%s-%s.txt" %(self.name1, self.name2) )
        f = open(outfile, 'w')
        f.write("Index")
        for s in sizes:
            f.write("\t%d\tStd" %s)
        metrics = self.metrics
        metricsStd = [m + 'Std' for m in metrics]

        for i, metric in enumerate(metrics):
            f.write("\n%s" %metric)
            for size in sizes:
                avr = size2stats[size][metric]
                std = size2stats[size][metricsStd[i]]
                f.write("\t%f\t%f" %(avr, std))
        f.write("\n")
        f.close()

        #pickle size2stats to temporary output directory
        picklefile = os.path.join(self.tempOutdir, "%s-%s.pickle" % (self.name1, self.name2))
        pickle.dump( size2stats, gzip.open(picklefile, "wb") )

class SamplingSamplePair(Target):
    '''Sampling the sample pair with a specific size for a specified number of times, report mean and std
    '''
    def __init__(self, outdir, sampleobjfile1, sampleobjfile2, size, options):
        Target.__init__(self)
        self.outdir = outdir
        self.sampleobjfile1 = sampleobjfile1
        self.sampleobjfile2 = sampleobjfile2
        self.size = size
        self.options = options

    def run(self):
        globalTempDir = self.getGlobalTempDir()
        for i in xrange(0, self.options.numSamplings):
            self.addChildTarget( SamplingPair(globalTempDir, self.sampleobjfile1, self.sampleobjfile2, self.size, i, self.options) )
        self.setFollowOnTarget( AverageSamplingPair(globalTempDir, self.outdir, self.size, self.options.similarityIndices) )
        
class SamplingPair(Target):
    '''One sampling of the sample with a specific size
    '''
    def __init__(self, outdir, sampleobjfile1, sampleobjfile2, size, id, options):
        Target.__init__(self)
        self.outdir = outdir
        self.sampleobjfile1 = sampleobjfile1
        self.sampleobjfile2 = sampleobjfile2
        self.size = size
        self.id = id
        self.options = options

    def run(self):
        sample1 = pickle.load( gzip.open(self.sampleobjfile1, "rb") )
        sample2 = pickle.load( gzip.open(self.sampleobjfile2, "rb") )
        
        print self.size, os.path.basename(self.sampleobjfile1).split('.')[0], sample1.total, os.path.basename(self.sampleobjfile2).split('.')[0], sample2.total
        #sampling "size" sequences from sample:
        subsample1 = iseqlib.samplingSample( sample1, self.size ) #sampled sample
        subsample2 = iseqlib.samplingSample( sample2, self.size ) #sampled sample
        
        #Get different statistics:
        stats = PairSamplingStats() #initialize stats for the current subsample

        #Find Sizes of the Union of clones between subsample1 and subsample2:
        counts1 = []
        counts2 = []
        for header, seq in subsample1.seqs.iteritems():
            counts1.append( seq.count )
            if header in subsample2.seqs:
                counts2.append( subsample2.seqs[header].count )
            else:
                counts2.append( 0 )
        for header, seq in subsample2.seqs.iteritems():
            if header not in subsample1.seqs:
                counts1.append( 0 )
                counts2.append( seq.count )
        
        if len(counts1) != len(counts2):
            raise ValueError("The count vectors of two compared samples must have equal length. Found %d and %d\n" %(len(counts1), len(counts2)))

        ## Initialize the R interface
        import rpy2.rinterface as rinterface
        rinterface.set_initoptions(('rpy2', '--no-save'))
        rinterface.initr()

        import rpy2.robjects as robjs
        from rpy2.robjects.packages import importr
        vegan = importr('vegan')
        vegdist = vegan.vegdist

        counts = counts1
        counts.extend(counts2)
        rcountsVec = robjs.IntVector( counts ) #convert python list into R vector
        rcountsMatrix = robjs.r['matrix'](rcountsVec, nrow = 2, byrow=True)
        
        #indices = ['bray', 'horn', 'mountford', 'chao']
        indices = self.options.similarityIndices
        #"manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn", "mountford", "raup" , "binomial" or "chao"
        for index in indices:
            disimilarity = vegdist( rcountsMatrix, method=index )
            stats[index] = disimilarity[0]
            
        #Write to temp file:
        picklefile = os.path.join( self.outdir, "%d.pickle" % self.id)
        pickle.dump( stats, gzip.open(picklefile, "wb") )

class AverageSamplingPair(Target):
    '''Mean and standard deviation of current sampling size
    '''
    def __init__(self, indir, outdir, size, metrics):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.samplingsize = size
        self.metrics = metrics

    def run(self):
        ext = 'pickle'
        picklefiles = iseqlib.getfiles( self.indir, ext )
        statsList = []
        for file in picklefiles:
            stats = pickle.load( gzip.open(os.path.join(self.indir, file), "rb") )
            statsList.append(stats)
        
        avrstats = PairSamplingStats() #initialize avrstats
        #Calculate mean and std using numpy
        #metrics = ['bray', 'horn', 'mountford', 'chao']
        #stds = ['brayStd', 'hornStd', 'mountfordStd', 'chaoStd']
        metrics = self.metrics
        stds = [m + "Std" for m in metrics]

        for i in xrange( len(metrics) ):
            vals = [s[metrics[i]] for s in statsList]
            avrstats[ metrics[i] ] = np.mean( vals )
            avrstats[ stds[i] ] = np.std( vals )

        #Pickle the average stat of this sampling size 
        picklefile = os.path.join(self.outdir, "%d.pickle" %self.samplingsize)
        pickle.dump( avrstats, gzip.open(picklefile, "wb") )

#===================== END PAIRWISE COMPARISONS, SIMILARITY INDICES ===================

class ReadFasta(Target):
    '''Read input fasta file and pickle the sample object to output directory
    '''
    def __init__(self, infile, outdir, name):
        Target.__init__(self)
        self.infile = infile
        self.outdir = outdir
        self.name = name

    def run(self):
        mode = 1 #seqs = { header: Seq() }
        mincount = 1
        seqs, total = iseqlib.readFile(self.infile, mincount, mode)
        total = sum([seq.count for seq in seqs.values()])
        for s in seqs.values():
            s.setFreq(total)

        sample = iseqlib.Sample(self.name)
        sample.seqs = seqs
        sample.setTotal(total)

        pickleFile = os.path.join(self.outdir, "%s.pickle" %self.name)
        pickle.dump( sample, gzip.open(pickleFile, "wb") )

####################### UTILITIES FUNCTIONS #########################
def getSamplingSizes(options, picklefile):
    samSizes = []
    sample = pickle.load( gzip.open(picklefile, "rb") )
    total = sample.total
    
    if options.samSize:
        for s in options.samSize:
            if s <= total:
                samSizes.append(s)
    else:
        if options.bin > total:
            samSizes = [total]
        else:
            samSizes = xrange(options.bin, total, options.bin)
    return samSizes


####################### OPTIONS ##############################
def checkOptions( parser, options, args ):
    #Check input directory
    if not os.path.exists( options.indir ):
        raise ValueError("Input directory %s does not exist\n" %options.indir)
    elif not os.path.isdir( options.indir ):
        raise ValueError("Input directory %s is not a directory\n" %options.indir)

    #Check output directory
    if not os.path.exists( options.outdir ):
        system( "mkdir %s\n" %(options.outdir) )
    
    #sampling sizes:
    if options.samSize:
        options.samSize = [long(s) for s in options.samSize.split(",")]

    options.diversityIndices = options.diversityIndices.split(',')
    options.similarityIndices = options.similarityIndices.split(',')

def addOptions( parser ):
    parser.add_option("-i", "--indir", dest = 'indir', help="Required argument. Input directory that contains input fasta files (each represents a sample)\n")
    parser.add_option("-o", "--outdir", dest = 'outdir', default=".", help="Output directory. Default: current directory\n")
    parser.add_option("-n", "--numSamplings", dest="numSamplings", type="int", default=100, help="Number of samplings. Default = %default")
    #parser.add_option("-s", "--samplingSize", dest="samSize", default="1000,10000,50000,100000,150000,200000,250000,300000,350000,400000,450000,500000,550000,600000,650000,700000,750000,800000,850000,900000,950000,1000000,1050000,1100000,1150000,1200000,1250000,1300000,1350000,1400000,1450000,1500000,1550000,1600000,1650000,1700000,1750000,1800000,1850000,1900000,1950000,2000000", help="Comma separated string of sampling sizes. Default=%default")
    parser.add_option('-b', '--bin', dest='bin', type='int', default=50000, help='Increment between sampling sizes. Default=%default. For example if bin=50k, will sample from the samples with size 50k, 100k, 150k, etc.')
    parser.add_option('-S', '--samplingSize', dest='samSize', help='Optional. Comma separted string of sampling sizes. If specified, will ignore "bin" and sample samples with these sizes. Default=None')
    parser.add_option('-d', '--diversity', dest='diversity', action='store_true', default=False, help='Diversity analyses if specified.')
    parser.add_option('-s', '--similarity', dest='similarity', action='store_true', default=False, help='Similarity analyses if specified.')
    parser.add_option('--diversityIndices', dest='diversityIndices', default='uniqClones,simpson,invsimpson,shannon,fisherAlpha', help='Comma separated list of diversity indices. Valid indices include uniqClones, simpson, invsimpson, shannon, fisherAlpha. Default=%default')
    parser.add_option('--similarityIndices', dest='similarityIndices', default='horn,mountford', help='Comma separated list of similarity indices. Valid indices include <manhattan, euclidean, canberra, bray, kulczynski, jaccard, gower, morisita, horn, mountford, raup, binomial, chao>. Default=%default')
    #parser.add_option("--maxSimPerRun", dest="maxSimPerRun", type='int', default=100, help="Maximum number of simulations per batch. Default=100")
    #parser.add_option("-s", "--numSamples", dest="numSamples", type="int", default=2, help="Number of samples (individuals/repertoires) to compare. Default=2.")

######################## MAIN ################################
def main():
    parser = iseqlib.initOptions()
    Stack.addJobTreeOptions(parser)

    addOptions(parser)

    options, args = parser.parse_args()
    checkOptions(parser, options, args)

    i = Stack( Setup(options) ).startJobTree(options)
    if i :
        raise RuntimeError("The jobtree contains %d failed jobs.\n" %i)

if __name__ == '__main__':
    from immunoseq.src.repertoireSurvey import *
    main()

