#!/usr/bin/env python2.6

#nknguyen at soe ucsc edu
#Feb 14 2012
#Immuno-seq pipeline

'''
Immuno seq pipepline
Inputs: 
    1/ Directory contains fasta files of productive CDR3 sequences (1 file per sample/patient)  in the following format:
        >sampleName;userDefinedCloneID|Vgenes|Jgenes;size=#
        Sequence
    Vgenes and Jgenes are comma separated, size = clone size (number of reads)
    Example: 
        >SBC1;Clusters6438|TRBV7-2|TRBJ2-3;size=1
        CASSFTQGETDTQYF

    2/ File that maps sample to patient, with format: <SampleName> <PatientName>
    Example: SBC2 B

    3/ 

Outputs/Functions:
    1/ V, J, VJ usage
    2/ Clone Sizes
    3/ CDR3 length distribution
    4/ Pairwise-overlapping
    5/ sequence Logo
    6/ (Differential expression of clones: heatmap where columns = samples and rows = clones.) 
    7/ Clustering...
'''

import os, sys, re, time
import xml.etree.ElementTree as ET

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import getTempDirectory
from sonLib.bioio import setLogLevel

#import immunoseq.lib.immunoseqLib as iseqlib
import immunoseqLib as iseqlib

##############################################################################################
#============= ADAPTIVE BIOTECHNOLOGIES FILE PROCESSING: TSVs TO FASTAs =====================#
##############################################################################################
class AdaptiveBioPreProcess( Target ):
    '''Process Adaptive Biotechnologies output files (.tsv) into productive CDR3 fasta files
    '''
    def __init__(self, indir, outdir, collapse, options):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.outdir = outdir
        self.collapse = collapse
        self.options = options

    def run(self):
        files = os.listdir(self.indir)
        #sequenceStatus = ['productive', 'nonProductive']
        status = 'productive'
        uniqv = False
        uniqj = False
        outdir = os.path.join( self.outdir, "%s-%s-%s" %(status, str(uniqv), str(uniqj)) )
        system("mkdir -p %s" %outdir)
        for f in files:
            self.addChildTarget( AdaptiveBioTsvToFa(self.indir, outdir, self.collapse, f, status, uniqv, uniqj) )
        self.setFollowOnTarget( MainPipeline(outdir, self.options) )

class AdaptiveBioTsvToFa( Target ):
    def __init__(self, indir, outdir, collapse, tsvfile, sequenceStatus, uniqV, uniqJ):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.outdir = outdir
        self.collapse = collapse
        self.tsvfile = tsvfile
        self.status = sequenceStatus
        self.uniqv = uniqV
        self.uniqj = uniqJ

    def run(self):
        samplename = self.tsvfile.rstrip(".tsv")
        localTempDir = self.getLocalTempDir()
        system("cp %s %s" %( os.path.join(self.indir, self.tsvfile), localTempDir ))
        command = "adaptiveTcrParse.py -i %s -o %s" %(localTempDir, localTempDir)
        if self.status == 'productive':
            command += " -p"
        elif self.status == 'nonProductive':
            command += " -b"
        else:
            sys.stderr.write("The sequence status is either 'productive' or 'nonProducitve'. Therefore print all sequences.\n")
        if self.uniqv:
            command += " -v"
        if self.uniqj:
            command += " -j"
        system(command)

        localfa = os.path.join( localTempDir, "%s.fa" %samplename )
        if not self.collapse:
            system("mv %s %s" %( localfa, self.outdir ))
        else:
            collapseFa = os.path.join( localTempDir, "%sCol.fa" %samplename )
            system("faCollapse.py -g %s %s" %(localfa, collapseFa))
            outfile = os.path.join(self.outdir, "%s.fa" %samplename)
            system("mv %s %s" %(collapseFa, outfile))

##############################################################################################
#========================= MAIN PIPELINE ====================================================#
##############################################################################################
class Setup( Target ):
    """Setting up jobs
    """
    def __init__(self, options):
        Target.__init__(self)
        self.options = options

    def run(self):
        if self.options.tsvdir:
            collapsed = True
            self.addChildTarget( AdaptiveBioPreProcess(self.options.tsvdir, self.options.outdir, collapsed, self.options) )
        else:
            self.addChildTarget( MainPipeline(self.options.seqdir, self.options) )

class MainPipeline(Target):
    def __init__(self, indir, options):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.options = options
    
    def run(self):
        #indir = self.options.seqdir
        indir = self.indir
        #sampling
        if self.options.sampling:
            samplingdir = os.path.join(self.options.outdir, "sampling", "%d" %self.options.sampling)
            system("mkdir -p %s" %samplingdir)
            self.addChildTarget( Sampling(self.indir, samplingdir, self.options.sampling, self.options.minReadCount) )
            indir = samplingdir #IF the sampling option was specified, use the sampled fasta files for all analyses 

        #Clone sizes
        clonesizeOutdir = os.path.join( self.options.outdir, "clonesize")
        system("mkdir -p %s" %(clonesizeOutdir))
        
        #--Clonesize Table
        self.addChildTarget( ClonesizeTab(indir, clonesizeOutdir, 1, self.options.group2samples) ) #all
        if self.options.minReadCount > 1: #only clones that pass minumum cutoffs
            self.addChildTarget( ClonesizeTab(indir, clonesizeOutdir, self.options.minReadCount, self.options.group2samples) )
        
        #--Clonesize Plots
        discretedir = os.path.join(clonesizeOutdir, "discrete")
        self.addChildTarget( ClonesizeDist(indir, discretedir, False, None) )
        cumulativedir = os.path.join(clonesizeOutdir, "cumulative")
        self.addChildTarget( ClonesizeDist(indir, cumulativedir, True, None) )
        
        if self.options.group2samples: #Average per group
            avrDiscretedir = os.path.join(clonesizeOutdir, "avrDiscrete")
            self.addChildTarget( ClonesizeDist(indir, avrDiscretedir, False, self.options.group2samples) )
            avrCumulativedir = os.path.join(clonesizeOutdir, "avrCumulative")
            self.addChildTarget( ClonesizeDist(indir, avrCumulativedir, True, self.options.group2samples) )

        #V, J, VJ usage:
        usageOutdir = os.path.join(self.options.outdir, "vjusage")
        system("mkdir -p %s" %usageOutdir)
        self.addChildTarget( VJusage(indir, usageOutdir, self.options.group2samples) )
        
        #CDR3 length distribution
        lendistOutdir = os.path.join( self.options.outdir, "lendist")
        system("mkdir -p %s" %(lendistOutdir))
        self.addChildTarget( Lendist(indir, lendistOutdir) )

        #Pairwise overlapping
        overlapOutdir = os.path.join( self.options.outdir, "overlap" ) #outdir/overlap
        system("mkdir -p %s" %overlapOutdir)
        #discrete
        discreteOverlapOutdir = os.path.join(overlapOutdir, 'discrete') #outdir/overlap/discrete
        system("mkdir -p %s" %discreteOverlapOutdir)
        self.addChildTarget( Overlap(indir, discreteOverlapOutdir, self.options.group2samples, self.options.minReadCount, self.options.minReadPercentage, True) )
        #cumulative
        cumulativeOverlapOutdir = os.path.join(overlapOutdir, 'cumulative') #outdir/overlap/cumulative
        system("mkdir -p %s" %cumulativeOverlapOutdir)
        self.addChildTarget( Overlap(indir, cumulativeOverlapOutdir, self.options.group2samples, self.options.minReadCount, self.options.minReadPercentage, False) )

class Sampling( Target ):
    '''Sampling a specific number of reads out of each sample.
    '''
    def __init__(self, indir, outdir, size, minReadCount):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.size = size
        self.minReadCount = minReadCount

    def run(self):
        cmd = "sampling.py -i %s -o %s -c %d -s %d" %(self.indir, self.outdir, selfminReadCount, self.size)
        system(cmd)

class ClonesizeTab( Target ):
    def __init__(self, indir, outdir, minReadCount, group2samples):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.outdir = outdir
        self.minReadCount = minReadCount
        self.group2samples = group2samples

    def run(self):
        #Latex table that summarize number of uniq clones, number of total reads, (simpson index)
        texfile = os.path.join(self.outdir, "clonesize-%d.tex" %self.minReadCount)
        cmd = "cloneSizeTab.py -o %s -c %d " %(texfile, self.minReadCount)
        if self.group2samples:
            cmd += " -g %s" %self.group2samples
        cmd += " %s" %self.indir
        system(cmd)

class ClonesizeDist(Target):
    def __init__(self, indir, outdir, isCumulative, group2samples):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.isCumulative = isCumulative
        self.group2samples = group2samples

    def run(self):
        system("mkdir -p %s" %self.outdir)
        cmd = "cloneSizeDist.py -o %s " %self.outdir
        if self.isCumulative:
            cmd += " -c"
        if self.group2samples:
            cmd += " -g %s" %self.group2samples
        cmd += " %s" %self.indir
        system(cmd)

class VJusage( Target ):
    '''Get V, J, VJ usage
    '''
    def __init__(self, indir, outdir, group2samples):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.group2samples = group2samples

    def run(self):
        #=====ALL samples=====:
        system("vjusage.py -i %s -o %s" %(self.indir, self.outdir))
        vjdir = os.path.join(self.outdir, "vj")
        
        #move average and std files to avr dir
        avrdir = os.path.join(vjdir, 'avrAll')
        system("mkdir -p %s" %avrdir)
        system("mv %s/average_* %s" %(vjdir, avrdir))
        system("mv %s/std_* %s" %(vjdir, avrdir))
        
        #usage plots and ttest:
        plotdir = os.path.join(self.outdir, "plots") #outdir/vjusage/plots
        system("mkdir -p %s" %plotdir)
        
        #uniq sequences, absolute count:
        absUniqDir = os.path.join(plotdir, "absUniq")
        absUniqVjDir = os.path.join(plotdir, "absUniq", 'vj')
        system("mkdir -p %s" %absUniqVjDir)
        system("mv %s/*_abs_uniq-vj.txt %s" %(vjdir, absUniqVjDir))
        self.addChildTarget( VJusagePlot(absUniqVjDir, absUniqDir, None, True, self.group2samples, True) )

        #uniq sequences, relative count:
        relUniqDir = os.path.join(plotdir, "relUniq")
        relUniqVjDir = os.path.join(plotdir, "relUniq", 'vj')
        system("mkdir -p %s" %relUniqVjDir)
        system("mv %s/*_uniq-vj.txt %s" % (vjdir, relUniqVjDir))
        self.addChildTarget( VJusagePlot(relUniqVjDir, relUniqDir, None, True, self.group2samples, True) )

        #reads, absolute count
        absAllDir = os.path.join(plotdir, "absAll")
        absAllVjDir = os.path.join(plotdir, "absAll", 'vj')
        system("mkdir -p %s" %absAllVjDir)
        system("mv %s/*_abs-vj.txt %s" %(vjdir, absAllVjDir))
        self.addChildTarget( VJusagePlot(absAllVjDir, absAllDir, None, True, self.group2samples, True) )

        #reads, relative count
        relAllDir = os.path.join(plotdir, "relAll")
        relAllVjDir = os.path.join(plotdir, "relAll", 'vj')
        system("mkdir -p %s" %relAllVjDir)
        system("mv %s/*-vj.txt %s" %( vjdir, relAllVjDir))
        self.addChildTarget( VJusagePlot(relAllVjDir, relAllDir, None, True, self.group2samples, True) )

        #=========GROUP AVERAGE ONLY=======:
        #read group2samples:
        group2samples = readGroup2samples(self.group2samples)
        #get vj usage:
        avrPlotdir = os.path.join(self.outdir, "averagePlots")
        system("mkdir -p %s" %avrPlotdir)
        
        stdfile = os.path.join(avrPlotdir, 'avr2std.txt')
        f = open(stdfile, 'w')

        absUniqDir = os.path.join(avrPlotdir, "absUniq") #outdir/averagePlots/absUniq
        relUniqDir = os.path.join(avrPlotdir, "relUniq") #outdir/averagePlots/relUniq
        absAllDir = os.path.join(avrPlotdir, "absAll") #outdir/averagePlots/absAll
        relAllDir = os.path.join(avrPlotdir, "relAll") #outdir/averagePlots/relAll
        
        for g, samples in group2samples.iteritems():
            groupdir = os.path.join(self.outdir, g, 'fasta')
            system("mkdir -p %s" %groupdir)
            for sample in samples:
                system("cp %s %s" %( os.path.join(self.indir, "%s.fa" %sample)   , groupdir))
            system("vjusage.py -i %s -o %s" %(groupdir, os.path.join(self.outdir, g)))
            groupvjdir = os.path.join(self.outdir, g, 'vj')
            #uniq sequences, absolute count:
            system("mkdir -p %s/vj" %absUniqDir)
            system("mv %s/average_abs_uniq-vj.txt %s/vj/%sAverage" %(groupvjdir, absUniqDir, g) )
            system("mv %s/std_abs_uniq-vj.txt %s/vj/std_%s" %(groupvjdir, absUniqDir, g))
            f.write("%sAverage std_%s\n" %(g, g))

            #uniq sequences, relative count:
            system("mkdir -p %s/vj" %relUniqDir)
            system("mv %s/average_uniq-vj.txt %s/vj/%sAverage" %(groupvjdir, relUniqDir, g) )
            system("mv %s/std_uniq-vj.txt %s/vj/std_%s" %(groupvjdir, relUniqDir, g))
            f.write("%sAverage std_%s\n" %(g, g))
            
            #all reads, absolute count:
            system("mkdir -p %s/vj" %absAllDir)
            system("mv %s/average_abs-vj.txt %s/vj/%sAverage" %(groupvjdir, absAllDir, g) )
            system("mv %s/std_abs-vj.txt %s/vj/std_%s" %(groupvjdir, absAllDir, g))
            f.write("%sAverage std_%s\n" %(g, g))
            
            #all reads, relative count:
            system("mkdir -p %s/vj" %relAllDir)
            system("mv %s/average-vj.txt %s/vj/%sAverage" %(groupvjdir, relAllDir, g) )
            system("mv %s/std-vj.txt %s/vj/std_%s" %(groupvjdir, relAllDir, g))
            f.write("%sAverage std_%s\n" %(g, g))
            
            system("rm -Rf %s" %groupdir)
        f.close()

        self.addChildTarget( VJusagePlot( os.path.join(absUniqDir, 'vj'), absUniqDir, stdfile, False, self.group2samples, True) )
        self.addChildTarget( VJusagePlot( os.path.join(relUniqDir, 'vj'), relUniqDir, stdfile, False, self.group2samples, True) )
        self.addChildTarget( VJusagePlot( os.path.join(absAllDir, 'vj'), absAllDir, stdfile, False, self.group2samples, True) )
        self.addChildTarget( VJusagePlot( os.path.join(relAllDir, 'vj'), relAllDir, stdfile, False, self.group2samples, True) )

class VJusagePlot( Target ):
    '''Plot V, J, VJ usage
    '''
    def __init__(self, indir, outdir, std, ttest, group2samples, vj):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.std = std
        self.ttest = ttest
        self.group2samples = group2samples
        self.vj = vj

    def run(self):
        cmd = 'usagePlot.py -i %s -o %s' %(self.indir, self.outdir)
        if self.vj:
            cmd += ' --vj'
        if self.std:
            cmd += ' --std %s' %self.std
        if self.ttest:
            cmd += ' --ttest'
        if self.group2samples:
            cmd += ' --group2samples %s' %self.group2samples
        system(cmd)
        #system("usagePlot.py -i %s -o %s --vj" %(self.indir, self.outdir))

class Lendist( Target ):
    def __init__(self, indir, outdir):
        Target.__init__(self, time = 0.00025)
        self.indir = indir
        self.outdir = outdir

    def run(self):
        files = iseqlib.getfiles(self.indir, 'fa')
        
        uniqOutdir = os.path.join(self.outdir, 'uniqSeqs')
        system("mkdir -p %s" %uniqOutdir)
        allOutdir = os.path.join(self.outdir, 'allReads')
        system("mkdir -p %s" %allOutdir)

        for file in files:
            filepath = os.path.join(self.indir, file)
            filename = file.rstrip('fa').rstrip('.')
            #
            uniqfile = os.path.join(uniqOutdir, '%s.txt' %filename)
            system("awk '$0 !~/>/ {print length($1)}' %s > %s" % (filepath, uniqfile))
            drawHist( uniqfile, os.path.join(uniqOutdir, '%s.pdf' %filename), True )
            #drawHist( uniqfile, os.path.join(uniqOutdir, '%s-abs.pdf' %filename), False )
            #
            allfile = os.path.join(allOutdir, '%s.txt' %filename)
            system("awk '{ if($0 ~/>/){ split($1, arr, \";size=\"); count = arr[2]}else{ for(i=0; i< count; i++){ print length($1)}} }' %s > %s" %(filepath, allfile) )
            drawHist( allfile, os.path.join(allOutdir, '%s.pdf' %filename), True )
            #drawHist( allfile, os.path.join(allOutdir, '%s-abs.pdf' %filename), False )

class Overlap( Target ):
    def __init__(self, indir, outdir, group2samples, minReadCount, cutoffs, discrete):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir #outdir/overlap/discrete[cumulative]
        self.group2samples = group2samples
        self.discrete = discrete
        self.minReadCount = minReadCount
        self.cutoffs = cutoffs
        
    def run(self):
        #Get pairwise overlap:
        getOverlapStats = True
        getPairwiseSharedSeqs = False
        getUniqAndSharedSeqs = None
        getOverlapPlot = False
        self.addChildTarget( Overlap2(self.indir, self.outdir, self.discrete, self.minReadCount, self.cutoffs, getOverlapStats, getPairwiseSharedSeqs, getUniqAndSharedSeqs, getOverlapPlot) )
        
        #Get pairwise shared sequences:
        getOverlapStats = False
        getPairwiseSharedSeqs = True
        getUniqAndSharedSeqs = None
        getOverlapPlot = False
        self.addChildTarget( Overlap2(self.indir, self.outdir, self.discrete, self.minReadCount, self.cutoffs, getOverlapStats, getPairwiseSharedSeqs, getUniqAndSharedSeqs, getOverlapPlot) )

        #Get uniq and shared sequences:
        getOverlapStats = False
        getPairwiseSharedSeqs = False
        getUniqAndSharedSeqs = self.group2samples
        getOverlapPlot = False
        self.addChildTarget( Overlap2(self.indir, self.outdir, self.discrete, self.minReadCount, self.cutoffs, getOverlapStats, getPairwiseSharedSeqs, getUniqAndSharedSeqs, getOverlapPlot) )

class Overlap2( Target ):
    def __init__(self, indir, outdir, discrete, minReadCount, cutoffs, getOverlapStats, getPairwiseSharedSeqs, group2samples, getOverlapPlot):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir #outdir/overlap/discrete[cumulative]/overlapStats
        self.discrete = discrete
        self.minReadCount = minReadCount
        self.cutoffs = cutoffs
        self.getOverlapStats = getOverlapStats
        self.getPairwiseSharedSeqs = getPairwiseSharedSeqs
        #self.getUniqAndSharedSeqs = getUniqAndSharedSeqs
        self.group2samples = group2samples
        self.getOverlapPlot = getOverlapPlot
        
    def run(self):
        #raise ValueError("MinPercentageCutoffs: %s\n" %self.cutoffs) #HACK
        cmd = "overlap.py -i %s -o %s -m 2 -c %d -p %s" %(self.indir, self.outdir, self.minReadCount, self.cutoffs) 
        if self.discrete:
            cmd += " -d"
        if not self.getOverlapStats:
            cmd += " --noPairwiseOverlap"
        if self.getPairwiseSharedSeqs:
            cmd += " -f"
        if self.getOverlapPlot:
            cmd += " --overlapPlot"
        if self.group2samples:
            cmd += " --group2samples %s" %self.group2samples
        system(cmd)

        #Set follow on jobs:
        if self.getOverlapStats:
            indir = os.path.join(self.outdir, "overlapStats")
            outdir = os.path.join(indir, "tables")
            system("mkdir -p %s" %outdir)
            self.setFollowOnTarget( OverlapStatTables(indir, outdir, self.group2samples, self.cutoffs) )
        #elif self.getPairwiseSharedSeqs:
        #    indir = os.path.join(self.outdir, "pairwiseSharedSeqs")
        #    outdir = os.path.join(indir, "")

class OverlapStatTables(Target):
    def __init__(self, indir, outdir, group2samples, cutoffs):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.group2samples = group2samples
        self.cutoffs = cutoffs

    def run(self):
        #mode 1: one statistics of interest across different cutoffs
        stats_types = ["clone1", "clone2", "cloneAvr", "read1", "read2", "readAvr"]
        for type in stats_types:
            cmd = "getOverlapTab.py -i %s -o %s -s %s -m 1 -l " %(self.indir, self.outdir, type)
            system(cmd)

        #mode 2: different statistics for each cutoff
        cmd = "getOverlapTab.py -i %s -o %s -m 2 -l -c %s" %(self.indir, self.outdir, self.cutoffs)
        if self.group2samples:
            cmd += " -g %s" %self.group2samples
        system(cmd)

def drawHist(file, outfile, rel):
    tempR = "%s-TEMP.R" %file
    f = open(tempR, 'w')
    f.write("data = read.table(\"%s\", header=FALSE)\n" % file)
    f.write("pdf('%s')\n" %outfile)
    if rel:
        f.write("hist(data[,1], freq=FALSE, lwd=5, lend='square', col='darkblue', xlab='Length (number of amino acids)', ylab='Frequency', main='CDR3 Length Distribution')\n" )
    else:
        f.write("hist(data[,1], freq=TRUE, lwd=5, lend='square', col='darkblue', xlab='Length (number of amino acids)', ylab='Number of sequences', main='CDR3 Length Distribution')\n" )
    f.write("dev.off()\n")
    f.close()
    system("R --no-save < %s" %tempR)
    system("rm %s" %tempR)

def readGroup2samples(file):
    group2samples = {}
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        items = line.split()
        if len(items) != 2:
            raise ValueError("Wrong group2samples file format. Required <group> <comma separated list of samples>. Got: %s\n" %line)
        group2samples[items[0]] = items[1].split(',')
    f.close()
    return group2samples

############## ERROR CLASSES ######################
class InputOptionError(Exception):
    pass

def checkOptions( parser, options, args ):
    if not options.seqdir:
        raise InputOptionError("Input directory is required, none was given.\n")
    if not os.path.isdir( options.seqdir ):
        raise InputOptionError("Input directory (%s) is not a directory\n" %options.seqdir)
    if not os.path.exists( options.outdir ):
        system( "mkdir %s\n" %(options.outdir) )
    #options.cutoffs = [ float(p) for p in options.minReadPercentage.strip(',').split(',') ]

def addOptions( parser ):
    parser.add_option("-s", "--sequences", dest = 'seqdir', help="Required argument. Input directory that contains fasta files of productive CDR3 sequences (each file represents a sample)\n")
    parser.add_option("-t", "--tsv", dest = 'tsvdir', help="Adaptive Biotechnologies data directory that contains tsv files (each file represents a sample). If specified, will run the preprocess steps that convert the tsv files to fasta files to seqdir. Default=%default\n")
    parser.add_option("-o", "--outdir", dest = 'outdir', default='immunoseqOut', help="Required argument. Output directory. Default=%default\n")
    #parser.add_option('-c', '--count', dest='minReadCount', default=1, type='int', help='Minimum read count cutoff. Default=%default')
    parser.add_option('--sampling', dest='sampling', type='int', help='Sampling size if wish to normalize the read counts of all samples (randomly select the specified size of reads from each sample). Default is no sampling.')
    parser.add_option("--group2samples", dest= 'group2samples', help='File specifying which samples belong to which group. Format: <group> <comma-separted-list-of-samples>')
    parser.add_option("--minReadCount", dest='minReadCount', type='int', default=1, help='Minimum read count: clones with < this number will be filtered out. Default=%default')
    parser.add_option("--minReadPercentage", dest='minReadPercentage', default='0,0.001,0.01,0.05,0.1,0.5,1,5', help='Comma separated list of read percentage cutoffs. Clone must have at least this percentage of total reads to be included. Default=%default')


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
    #from immunoseq import *
    from immunoseq.src.immunoSeq import *
    main()

