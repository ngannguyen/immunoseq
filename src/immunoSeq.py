#!/usr/bin/env python2.6

#nknguyen at soe ucsc edu
#Tue Jul 17 09:09:24 PDT 2012
#Immuno-seq pipeline

'''
Immuno seq pipepline
Inputs: 
    1/ Directory contains fasta files of productive CDR3 sequences (1 file per sample/patient)  in the following format:
        >sampleName;userDefinedCloneID|Vgenes|Jgenes|Dgenes;size=#
        Sequence
    Vgenes, Jgenes and D genes are comma separated, size = clone size (number of reads)
    Example: 
        >SBC1;Clusters6438|TRBV7-2|TRBJ2-3|TRBD1-1;size=1
        CASSFTQGETDTQYF

    2/ File that maps sample to group, with format: <Group> <sample1,sample2,...>
    Example:
        patients as1R,as2R,as3R
        controls c1,c2,c3

Outputs/Functions:
    1/ Gene usage analysis: V, J, D, DJ, VJ, VDJ usage
    2/ CDR3 length distribution
    3/ Amino acid usage
    4/ Nucleotide to amino acid convergence
       Nucleotide/aa to VDJ or VJ
    5/ Clonesize distribution (#reads/clone vs #clones)
    
    #(Differential expression of clones: heatmap where columns = samples and rows = clones.) 
    #Clustering...
'''

import os, sys, re, time, gzip
import cPickle as pickle
from optparse import OptionParser
#import xml.etree.ElementTree as ET

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import logger
from sonLib.bioio import system
#from sonLib.bioio import getTempDirectory
#from sonLib.bioio import setLogLevel

#import immunoseq.lib.immunoseqLib as iseqlib
import immunoseqLib as iseqlib
from immunoseq.lib.geneUsageLib import *
from immunoseq.lib.lendistLib import *
from immunoseq.lib.aausageLib import *
import numpy as np

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
    '''Start of the main pipeline.
       Firstly, read input fasta files
       When done, call following steps to do sampling if needed and all the analyses
    '''
    def __init__(self, indir, options):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.options = options
    
    def run(self):
        #read input fasta files:
        ext = 'fa'
        files = iseqlib.getfiles(self.indir, ext)
        globalTempDir = self.getGlobalTempDir()
        for file in files:
            filepath = os.path.join(self.indir, file)
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
        ext = "pickle"
        files = iseqlib.getfiles(self.indir, ext)
        samples = [file.split('.')[0] for file in files]

        for sample in samples:
            samplefile = os.path.join(self.indir, "%s.%s"%(sample, ext))
            sampledir = os.path.join(globalTempDir, sample)
            system("mkdir -p %s" %sampledir)
            if self.options.sampling:
                for i in xrange(self.options.numsam): #sampling a number of times
                    samplingdir = os.path.join(sampledir, "%d" %i)
                    system("mkdir -p %s" %samplingdir)
                    self.addChildTarget( Sampling(samplefile, samplingdir, self.options) )
            else:
                tempoutdir = os.path.join(sampledir, "0")
                system("mkdir -p %s" %tempoutdir)
                
                #filtering if selected Vs and/or selected Js were specified
                if self.options.vs or self.options.js:
                    sampleObj = pickle.load( gzip.open(samplefile, "rb") )
                    subsample = iseqlib.filterSampleByGenes(sampleObj, self.options.vs, self.options.js)
                    system("rm %s" %samplefile)
                    pickle.dump( subsample, gzip.open(samplefile, "wb") )

                self.addChildTarget( Analyses(samplefile, tempoutdir, self.options) )
        #Calculate means & standard deviations of samplings
        self.setFollowOnTarget( AverageResults(globalTempDir, self.options) )

class Sampling(Target):
    '''Sampling and and call analyses
    '''
    def __init__(self, samplefile, outdir, options):
        Target.__init__(self)
        self.samplefile = samplefile
        self.outdir = outdir
        self.options = options
        self.size = options.sampling

    def run(self):
        globalTempDir = self.getGlobalTempDir()
        #sampling
        sample = pickle.load( gzip.open(self.samplefile, "rb") )
        if self.options.uniq:
            subsample = iseqlib.samplingSample_weightedUniq(sample, size)
        else:
            subsample = iseqlib.samplingSample(sample, self.size)
        
        #filtering if selected Vs and/or selected Js were specified
        subsample = iseqlib.filterSampleByGenes(subsample, self.options.vs, self.options.js)

        picklefile = os.path.join(globalTempDir, "%s.pickle" %sample.name)
        pickle.dump(subsample, gzip.open(picklefile, "wb"))
        self.addChildTarget( Analyses(picklefile, self.outdir, self.options) )

class  Analyses(Target):
    '''Call jobs to compute different analyses
    '''
    def __init__(self, infile, outdir, options):
        Target.__init__(self)
        self.infile = infile
        self.outdir = outdir
        self.options = options

    def run(self):
        #analyses = ['geneusage', 'cdr3len', 'aausage', 'seq2vdj', 'clonesize']
        for analysis in self.options.analyses:
            outdir = os.path.join(self.outdir, analysis)
            system("mkdir -p %s" %outdir)
            self.addChildTarget( Analysis(self.infile, outdir, analysis) )

####################################
############# ANALYSES #############
####################################
class Analysis(Target):
    '''
    '''
    def __init__(self, infile, outdir, type):
        Target.__init__(self)
        self.infile = infile
        self.outdir = outdir
        self.type = type

    def run(self):
        #analyses = ['geneusage', 'cdr3len', 'aausage', 'clonesize']
        sample = pickle.load( gzip.open(self.infile, "rb") )
        if len(sample.seqs) == 0:
            #raise ValueError("Sample %s has zero sequence!. Pickle file is %s" %(sample.name, self.infile))
            return


        if self.type == 'geneusage':#V, D, J, DJ, VJ, VDJ usage:
            obj = getGene2count(sample.seqs) #type2gene2count
        elif self.type == 'cdr3len':
            obj =  getLenDist(sample.seqs) #len2freq
        elif self.type == 'aausage':
            sizes = [1, 2, 3] #single-aa, di-aa, tri-aa, 4-aa
            #sizes = [1, 2, 3, 4] #single-aa, di-aa, tri-aa, 4-aa
            #sizes = [1]
            uniq = True
            obj_uniq =  len2aaSampleProfile(sample, sizes, uniq) #len2size2pos2aa2usage
            obj_read =  len2aaSampleProfile(sample, sizes, not uniq) #len2size2pos2aa2usage
            obj = [obj_uniq, obj_read]
        #elif self.type == 'clonesize':
        #    obj = 
        else:
            raise ValueError("%s is an invalid analysis type" %self.type)

        pickleFile = os.path.join(self.outdir, "%s.pickle" %sample.name)
        pickle.dump(obj, gzip.open(pickleFile, 'wb'))

##########################################################
####### COMPUTE MEANS & STDS OF THE SAMPLINGS ############
##########################################################
class AverageResults(Target):
    '''
    '''
    def __init__(self, indir, options):
        Target.__init__(self)
        self.indir = indir
        self.options = options

    def run(self):
        analyses = self.options.analyses
        #analyses = ['geneusage', 'cdr3len', 'aausage', 'clonesize']
        if 'geneusage' in analyses:
            self.addChildTarget( GeneUsage(self.indir, self.options) )
        if 'cdr3len' in analyses:
            self.addChildTarget( Cdr3len(self.indir, self.options) )
        if 'aausage' in analyses:
            uniq = True
            self.addChildTarget( Aausage(self.indir, self.options, uniq) )
            self.addChildTarget( Aausage(self.indir, self.options, not uniq) )

########################################
######## GENE USAGE ####################
########################################
class GeneUsage(Target):
    '''
    '''
    def __init__(self, indir, options):
        Target.__init__(self)
        self.indir = indir
        self.options = options

    def run(self):
        type2genelist = getInputGenelists(self.options)
        outdir = os.path.join(self.options.outdir, "geneusage")
        alldir = os.path.join(outdir, 'all')
        geneuseddir = os.path.join(alldir, "geneused")
        system("mkdir -p %s" % geneuseddir)
        #types = ['v', 'j', 'd', 'dj', 'vdj']
        types = ['v', 'j']
        if not self.options.noDs:
            types.append('d')
        
        for type in types:
            typedir = os.path.join(alldir, type)
            system("mkdir -p %s" %typedir)   #example: outdir/geneusage/all/v/
        
        if self.options.noDs:
            combtypes = ['vj']
        else:
            combtypes = ['vj', 'dj', 'vdj']
        vjtypes = ['abs', 'rel', 'absuniq', 'reluniq'] #ex: outdir/geneusage/all/vj/absuniq
        for comb in combtypes:
            for t in vjtypes:
                system("mkdir -p %s" % os.path.join(alldir, comb, t)) #ex: outdir/geneusage/all/dj/uniq

        samples = os.listdir(self.indir)
        samples = [s.split('.')[0] for s in samples]
        sample2mean = {}
        for sample in samples:
            sampledir = os.path.join(self.indir, sample)
            aggType2gene2count = {}
            
            #get stats for all samplings
            for i in xrange(self.options.numsam):
                picklefile = os.path.join(sampledir, "%d" %i, "geneusage", "%s.pickle" %sample)
                if not os.path.exists(picklefile):
                    type2gene2count = {}
                else:
                    type2gene2count = pickle.load( gzip.open(picklefile, 'rb') )

                #update aggType2gene2count
                addSamplingStats(type2gene2count, aggType2gene2count, i)
            
            #Average stats of the samplings:
            avrstats, stdstats = avrSamplingStats(aggType2gene2count)
            
            sample2mean[sample] = avrstats
            if avrstats:
                #Print v, j, d, dj, vdj
                usageTab(types, sample, avrstats, stdstats, type2genelist, alldir)
                #Percentage of total genes used/observed:
                outfile = os.path.join(geneuseddir, "%s.txt" %sample)
                geneUsed(avrstats, type2genelist, outfile)
                #Print vj, dj, vdj
                for comb in combtypes:
                    combdir = os.path.join(alldir, comb)
                    rowtype = comb[0]
                    coltype = comb[1:]
                    #print avrstats[comb]
                    if self.options.numsam == 1:
                        getVJusageSample(sample, rowtype, coltype, avrstats, None, type2genelist, combdir)
                    else:
                        getVJusageSample(sample, rowtype, coltype, avrstats, stdstats, type2genelist, combdir)
            
        #Combine geneused:
        abs = True
        geneusedfile = os.path.join(geneuseddir, 'summary.txt')
        geneUsedSummary(sample2mean, type2genelist, self.options.group2samples, geneusedfile, not abs)
        geneusedAbsfile = os.path.join(geneuseddir, 'summary-abs.txt')
        geneUsedSummary(sample2mean, type2genelist, self.options.group2samples, geneusedAbsfile, abs)

        #GROUP AVERAGE
        avrdir = os.path.join(outdir, 'avr')
        system("mkdir -p %s" %avrdir)
        for type in types:
            typedir = os.path.join(avrdir, type)
            system("mkdir -p %s" %typedir)
        for comb in combtypes:
            for t in vjtypes:
                system("mkdir -p %s" % os.path.join(avrdir, comb, t)) #ex: outdir/geneusage/avr/dj/uniq

        for group, groupsamples in self.options.group2samples.iteritems():
            groupstats = {}
            for i, s in enumerate(groupsamples):
                type2gene2count = sample2mean[s]
                addSamplingStats(type2gene2count, groupstats, i)
            avrstats, stdstats = avrSamplingStats(groupstats)
            #Print v, j, d, dj, vdj
            usageTab(types, group, avrstats, stdstats, type2genelist, avrdir)
            #Print vj, dj, vdj
            for comb in combtypes:
                combdir = os.path.join(avrdir, comb)
                rowtype = comb[0]
                coltype = comb[1:]
                if len(groupsamples) == 1:
                    getVJusageSample(group, rowtype, coltype, avrstats, None, type2genelist, combdir)
                else:
                    getVJusageSample(group, rowtype, coltype, avrstats, stdstats, type2genelist, combdir)

        #Call to draw plots:
        if not self.options.vs and not self.options.js:
            self.setFollowOnTarget(GeneUsagePlot(outdir, self.options))

class GeneUsagePlot(Target):
    '''Plot V usage, J usage, VJ usage, and print some summary
    '''
    def __init__(self, indir, options):
        Target.__init__(self)
        self.indir = indir
        self.group2samplesFile = options.group2samplesFile
        self.targetGroup = options.ttestTargetGroup
        self.pval = options.pval
        self.sampleOrder = options.sampleOrder
        self.noDs = options.noDs

    def run(self):
        outdir = os.path.join(self.indir, "plots")
        system("mkdir -p %s" %outdir)
        levels = ['all', 'avr']
        if self.noDs:
            combs = ['vj']
        else:
            combs = ['vj', 'dj', 'vdj']
        
        for level in levels:
            for comb in combs:
                vjtypes = ['abs', 'rel', 'absuniq', 'reluniq']
                for type in vjtypes:
                    typedir = os.path.join(self.indir, level, comb, type) #ex: indir/all/vj/abs/
                    typeoutdir = os.path.join(outdir, level, comb, type) #ex: plots/all/vj/abs/
                    system("mkdir -p %s" %typeoutdir)
                    cmd = 'usagePlot.py -i %s -o %s --group2samples %s --vj -p %f ' %(typedir, typeoutdir, self.group2samplesFile, self.pval)
                    if self.sampleOrder:
                        cmd += '--sampleOrder %s ' % self.sampleOrder
                    if level == 'all':
                        cmd += ' --ttest '
                        if self.targetGroup:
                            cmd += ' --ttestTargetGroup %s ' %self.targetGroup
                    if type == 'abs' or type == 'absuniq':
                        cmd += ' -a '
                    system(cmd)

##########################################
########## CDR3 LENGTH DISTRIBUTION ######
##########################################
class Cdr3len(Target):
    '''Length distribution
    '''
    def __init__(self, indir, options):
        Target.__init__(self)
        self.indir = indir
        self.options = options

    def run(self):
        sample2avr = {}
        sample2std = {}
        samples = os.listdir(self.indir)
        samples = [s.split('.')[0] for s in samples]
        for sample in samples:
            sampledir = os.path.join(self.indir, sample)
            agglen2freq = {}
            #get stats for all samplings
            for i in xrange(self.options.numsam):
                picklefile = os.path.join(sampledir, "%d" %i, "cdr3len", "%s.pickle" %sample)
                if not os.path.exists(picklefile):
                    len2freq = {}
                else:
                    len2freq = pickle.load( gzip.open(picklefile, 'rb') )
                #update aggType2gene2count
                addSamplingLenDist(len2freq, agglen2freq, i)
            
            #Average stats of the samplings:
            l2avr, l2std = avrSamplingLenDist(agglen2freq)
            sample2avr[sample] = l2avr
            sample2std[sample] = l2std
        #Get union lendist:
        getUnionLenDist(sample2avr)
        getUnionLenDist(sample2std)

        #Draw lendist:
        outdir = os.path.join(self.options.outdir, "lendist")
        system("mkdir -p %s" %outdir)
        uniq = True
        uniqoutfile = os.path.join(outdir, 'lendist_uniq')
        drawLenDist2(sample2avr, uniqoutfile, uniq)
        readoutfile = os.path.join(outdir, 'lendist_read')
        drawLenDist2(sample2avr, readoutfile, not uniq)

        #Print tables:
        taboutdir = os.path.join(outdir, 'tables')
        system("mkdir -p %s" %taboutdir)
        printLenDistTab(sample2avr, sample2std, taboutdir)

        #t-test
        if self.options.group2samples:
            uniq = True
            ttestoutfile = os.path.join(outdir, 'ttests.txt')
            lenDistTtests(sample2avr, self.options.group2samples, ttestoutfile, self.options.ttestTargetGroup, self.options.pval)

        #clonesize distribution for each length

##############################################
######### AMINO ACID USAGE ###################
##############################################
class Aausage(Target):
    '''Amino acid usage
    '''
    def __init__(self, indir, options, uniq):
        Target.__init__(self)
        self.indir = indir
        self.options = options
        self.uniq = uniq

    def run(self):
        #len2size2pos2aa2usage
        sample2avr = {}
        sample2std = {}
        samples = os.listdir(self.indir)
        samples = [s.split('.')[0] for s in samples]
        for sample in samples:
            sampledir = os.path.join(self.indir, sample)
            aggstats = {}
            #get stats for all samplings
            for i in xrange(self.options.numsam):
                picklefile = os.path.join(sampledir, "%d" %i, "aausage", "%s.pickle" %sample)
                if not os.path.exists(picklefile):
                    stats = {}
                else:
                    statsList = pickle.load( gzip.open(picklefile, 'rb') )
                    if self.uniq:
                        stats = statsList[0]
                    else:
                        stats = statsList[1]
                #update aggType2gene2count
                addSamplingAausage(stats, aggstats, i)
            
            #Average stats of the samplings:
            avrstats, stdstats = avrSamplingAausage(aggstats)
            sample2avr[sample] = avrstats
            sample2std[sample] = stdstats

        #print Transfac format:
        if self.uniq:
            outdir = os.path.join(self.options.outdir, 'aausage', "uniq")
        else:
            outdir = os.path.join(self.options.outdir, 'aausage', "read")
        system('mkdir -p %s' %outdir)
        logoTextdir = os.path.join(outdir, 'logotextfiles')
        system('mkdir -p %s' %logoTextdir)
        logodir = os.path.join(outdir, 'logo')
        for sample in samples:
            len2size2pos2aausage = sample2avr[sample]
            for l in len2size2pos2aausage:
                textdir = os.path.join(logoTextdir, "%d" %l)
                system('mkdir -p %s' %textdir)
                if 1 in len2size2pos2aausage[l]:
                    pos2aausage = len2size2pos2aausage[l][1]
                    outfile = os.path.join(textdir, '%s.txt' %sample)
                    printTransfacFormat(outfile, sample, pos2aausage)
                    plotdir = os.path.join(logodir, "%d" %l)
                    system('mkdir -p %s' %plotdir)
                    logofile = os.path.join(plotdir, '%s.pdf' %sample)
                    title = "%s" %sample
                    xlabel = "Position"
                    ylabel = "Bits"
                    cmd = '/cluster/software/bin/weblogo -F pdf -D transfac -s large -t "%s" -x "%s" -y "%s" --errorbars NO < %s > %s' %(title, xlabel, ylabel, outfile, logofile)
                    system(cmd)

        #print matrix where columns = samples, rows = position+aa
        matrixdir = os.path.join(outdir, 'tables')
        system('mkdir -p %s' % matrixdir)
        getAaUsage(sample2avr, matrixdir)

        #ttests:
        if self.options.ttestTargetGroup:
            ttestdir = os.path.join(outdir, 'ttests')
            system("mkdir -p %s" %ttestdir)
            targetSamples = self.options.group2samples[self.options.ttestTargetGroup]
            aausage_ttests(ttestdir, sample2avr, targetSamples, self.options.pval)

def getInputGenelists(options):
    type2list = {}
    if options.vgenes:
        type2list['v'] = options.vgenes
    if options.jgenes:
        type2list['j'] = options.jgenes
    if options.dgenes:
        type2list['d'] = options.dgenes
    if options.dgenes and options.jgenes:
        dj = []
        for d in options.dgenes:
            for j in options.jgenes:
                dj.append("%s|%s" %(d, j))
        type2list['dj'] = dj
    if options.vgenes and options.jgenes:
        vj = []
        for v in options.vgenes:
            for j in options.jgenes:
                vj.append("%s|%s" %(v, j))
        type2list['vj'] = vj
    if options.vgenes and options.jgenes and options.dgenes:
        vdj = []
        for v in options.vgenes:
            for d in options.dgenes:
                for j in options.jgenes:
                    vdj.append("%s|%s|%s" %(v, d, j))
        type2list['vdj'] = vdj
    return type2list

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
    analyses = ['geneusage', 'cdr3len', 'aausage']
    if options.analyses == 'all':
        options.analyses = analyses
    else:
        options.analyses = options.analyses.split(',')
        for analysis in options.analyses:
            if analysis not in analyses:
                raise ValueError("Option for analysis %s is not valid." %analysis)
    if options.vgenes:
        options.vgenes = iseqlib.readList(options.vgenes, ',')
    if options.jgenes:
        options.jgenes = iseqlib.readList(options.jgenes, ',')
    if options.dgenes:
        options.dgenes = iseqlib.readList(options.dgenes, ',')
    if options.group2samplesFile:
        options.group2samples, options.sample2group = iseqlib.readGroup2samples(options.group2samplesFile)
    if options.vs:
        options.vs = options.vs.split(',')
    if options.js:
        options.js = options.js.split(',')
    #if options.sampleOrder:
    #    options.sampleOrder = options.sampleOrder.split(',')

def addOptions( parser ):
    parser.add_option("-s", "--sequences", dest = 'seqdir', help="Required argument. Input directory that contains fasta files of productive CDR3 sequences (each file represents a sample)\n")
    parser.add_option("-t", "--tsv", dest = 'tsvdir', help="Adaptive Biotechnologies data directory that contains tsv files (each file represents a sample). If specified, will run the preprocess steps that convert the tsv files to fasta files to seqdir. Default=%default\n")
    parser.add_option("-o", "--outdir", dest = 'outdir', default='immunoseqOut', help="Required argument. Output directory. Default=%default\n")
    parser.add_option('-u', '--uniq', dest='uniq', action='store_true', default=False, help='If specified, will perform the weighted-unique (clones) samplings insteads of reads samplings')
    parser.add_option('--sampling', dest='sampling', type='int', help='Sampling size if wish to normalize the read counts of all samples (randomly select the specified size of reads from each sample). Default is no sampling.')
    parser.add_option("--numsam", dest='numsam', type='int', default=1, help='If sampling is specified, this is the number of samplings to perform for each sample')
    parser.add_option("--group2samples", dest= 'group2samplesFile', help='File specifying which samples belong to which group. Format: <group> <comma-separted-list-of-samples>')
    parser.add_option('--ttestTargetGroup', dest='ttestTargetGroup', help='If specified, will compare samples belong to this group versus all other samples')
    parser.add_option("--minReadCount", dest='minReadCount', type='int', default=1, help='Minimum read count: clones with < this number will be filtered out. Default=%default')
    parser.add_option("--minReadPercentage", dest='minReadPercentage', default='0,0.001,0.01,0.05,0.1,0.5,1,5', help='Comma separated list of read percentage cutoffs. Clone must have at least this percentage of total reads to be included. Default=%default')
    parser.add_option("--analyses", dest='analyses', default='all', help='Comma separated list of analyses to do. Valid values are[geneusage, cdr3len, aausage, clonesize]. Default="all", meaning all analyses')
    parser.add_option("--vgenes", dest='vgenes', default='/hive/users/nknguyen/immuno/imgt/TRB/v2.lst', help='File contain known vgenes' )
    parser.add_option("--jgenes", dest='jgenes', default='/hive/users/nknguyen/immuno/imgt/TRB/j.lst', help='File contain known jgenes' )
    parser.add_option("--dgenes", dest='dgenes', default='/hive/users/nknguyen/immuno/imgt/TRB/d.lst', help='File contain known dgenes' )
    parser.add_option('-p', '--pval', dest='pval', default=1.0, type='float', help='pvalue cutoff. Only report ttests with pval <= this value')
    parser.add_option('--sampleOrder', dest='sampleOrder', type='string', help='Comma seperated list of samples in the order of interest')
    #Filtering options (only look at clones with at least 1 V and at least 1 J in the specified lists
    parser.add_option('--vs', dest='vs', help='Comma separated list of selected V genes')
    parser.add_option('--js', dest='js', help='Comma separated list of selected J genes')
    parser.add_option('--noDs', dest='noDs', action='store_true', default=False, help='If specified, ignore analyses involved in D genes')

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
    from immunoseq.src.immunoSeq import *
    main()

