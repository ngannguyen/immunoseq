#nknguyen soe ucsc edu
#Tue Jul 17 12:48:10 PDT 2012
#

import sys, os, re, time
from optparse import OptionParser
from sonLib.bioio import system
from immunoseq.lib.immunoseqLib import *
import numpy as np
from scipy.stats.stats import pearsonr
from scipy.stats import ttest_ind

def getLen2seqsSample(sample):
    len2seqs = {}
    for header, seq in sample.seqs.iteritems():
        l = len(seq.seq)
        if l not in len2seqs:
            len2seqs[l] = {header:seq}
        else:
            len2seqs[l][header] = seq
    #sample.len2seqs = len2seqs
    return len2seqs

def setLen2seqs(samples):
    for sample in samples:
        len2seqs = getLen2seqsSample(sample)
        sample.len2seqs = len2seqs

def aaSampleProfile(seqs, length, sizes, uniq):
    size2pos2aa2usage = {}
    if uniq:
        totalCount = len(seqs)
    else:
        totalCount = sum([seq.count for seq in seqs.values()])
    #sizes = [1, 2, 3, 4] #single-aa, di-aa, tri-aa, 4-aa
    
    for size in sizes:
        size2pos2aa2usage[size] = {}
        for i in xrange(0, length - size + 1): #position
            aa2usage = {}
	    for seq in seqs.values():
                aa = seq.seq[i:i+size] #aa of current sequence
                count = 1.0
                if not uniq:
                    count = seq.count
                if aa not in aa2usage:
                    aa2usage[aa] = count
                else:
                    aa2usage[aa] += count
            #Normalize:
            for aa in aa2usage:
                aa2usage[aa] /= float(totalCount)
            size2pos2aa2usage[size][ '-'.join( [str(p) for p in xrange(i, i+size)] ) ] = aa2usage 
    
    return size2pos2aa2usage

def len2aaSampleProfile(sample, sizes, uniq):
    len2size2pos2aa2usage = {}
    len2seqs = getLen2seqsSample(sample)
    for l in len2seqs.keys():
        size2pos2aa2usage = aaSampleProfile(len2seqs[l], l, sizes, uniq)
        len2size2pos2aa2usage[l] = size2pos2aa2usage
    return len2size2pos2aa2usage

def addSamplingAausage(len2size2pos2aa2usage, aggstats, i):
    for l, size2stats in len2size2pos2aa2usage.iteritems():
        if l not in aggstats:
            aggstats[l] = {}
            for size, pos2aausage in size2stats.iteritems():
                aggstats[l][size] = {}
                for pos, aausage in pos2aausage.iteritems():
                    aggstats[l][size][pos] = {}
                    for aa, count in aausage.iteritems():
                        aggstats[l][size][pos][aa] = [0]*i + [count]
        else:
            for size, pos2aausage in size2stats.iteritems():
                for pos, aausage in pos2aausage.iteritems():
                    for aa, count in aausage.iteritems():
                        if aa not in aggstats[l][size][pos]:
                            aggstats[l][size][pos][aa] = [0]*i + [count]
                        else:
                            aggstats[l][size][pos][aa].append(count)

def avrSamplingAausage(aggstats):
    avrstats = {}
    stdstats = {}
    for l, size2stats in aggstats.iteritems():
        avrstats[l] = {}
        stdstats[l] = {}
        for size, pos2aausage in size2stats.iteritems():
            avrstats[l][size] = {}
            stdstats[l][size] = {}
            for pos, aausage in pos2aausage.iteritems():
                avrstats[l][size][pos] = {}
                stdstats[l][size][pos] = {}
                for aa, countVec in aausage.iteritems():
                    avrstats[l][size][pos][aa] = np.mean(countVec)
                    stdstats[l][size][pos][aa] = np.std(countVec)
    #DEBUG
    #for l, size2stats in avrstats.iteritems():
    #    print l
    #    for size, pos2aausage in size2stats.iteritems():
    #        print size
    #        for pos, aausage in pos2aausage.iteritems():
    #            print pos
    #            for aa, usage in aausage.iteritems():
    #                print aa, usage
    #raise ValueError("DEBUGGING...")
    #END DEBUG
    return avrstats, stdstats

def printTransfacFormat(outfile, sample, pos2aa2count):
    f = open(outfile, 'w')
    f.write("ID\t%s\n" %sample)
    f.write("BF\t%s\n" %sample)
    aa = 'ACDEFGHIKLMNPQRSTVWY'
    aaletters = [l for l in aa]
    f.write("P0\t%s\n" %('\t'.join(aaletters)) )
    positions = sorted( [int(p) for p in pos2aa2count.keys()] )
    for pos in positions:
        aa2count = pos2aa2count[ str(pos) ]
        f.write("%d" % (pos + 1) ) #convert pos to base 1
        for a in aaletters:
            count = 0
            if a in aa2count:
                count = aa2count[a]
            f.write( "\t%f" % (count*10000) )
        f.write("\n")
    f.write("XX\n")
    f.write("//\n")
    f.close()

def aaProfilePerLen(sample2len2size2pos2aausage, length):
    size2pos2aa2sample2usage = {}
    sizes = [1, 2, 3, 4] #single-aa, di-aa, tri-aa, 4-aa
    for sample, stats in sample2len2size2pos2aausage.iteritems():
        if length not in stats:
            continue

        #seqs = sample.len2seqs[length]
        #size2pos2aa2usage =  aaSampleProfile(seqs, length, sizes)
        size2pos2aa2usage = stats[length]

        #add to pos2aa2sample2usage
        for size, pos2aa2usage in size2pos2aa2usage.iteritems():
            if size not in size2pos2aa2sample2usage:
                size2pos2aa2sample2usage[size] = {}
            pos2aa2sample2usage = size2pos2aa2sample2usage[size]
            for pos, aa2usage in pos2aa2usage.iteritems():
                if pos not in pos2aa2sample2usage:
                    pos2aa2sample2usage[pos] = {}
                    for aa, usage in aa2usage.iteritems():
                        pos2aa2sample2usage[pos][aa] = {sample: usage}
                else:
                    for aa, usage in aa2usage.iteritems():
                        if aa not in pos2aa2sample2usage[pos]:
                            pos2aa2sample2usage[pos][aa] = {sample: usage}
                        else:
                            pos2aa2sample2usage[pos][aa][sample] = usage
            #size2pos2aa2sample2usage[size] = pos2aa2sample2usage

    return size2pos2aa2sample2usage

def printProfile(size2pos2aa2sample2usage, sampleNames, outdir):
    for size, pos2aa2sample2usage in size2pos2aa2sample2usage.iteritems():
        outfile = os.path.join(outdir, "%d.txt" %size)
        f = open(outfile, 'w')
        f.write("#Pos-AA\t%s\n" %( "\t".join(sampleNames) ))
        for pos, aa2sample2usage in pos2aa2sample2usage.iteritems():
            for aa, sample2usage in aa2sample2usage.iteritems():
                rowname = "%s_%s" %(pos, aa)
                usage = []
                for sample in sampleNames:
                    if sample not in sample2usage:
                        usage.append('0')
                    else:
                        usage.append( str( sample2usage[sample] ) )
                f.write("%s\t%s\n" %(rowname, "\t".join(usage) ) )
        f.close()

def getAaUsage(sample2stats, outdir):
    lens = xrange(12, 18)
    sampleNames = sorted( sample2stats.keys() )

    for l in lens:
        lendir = os.path.join(outdir, "%d" %l)
        system("mkdir -p %s" % lendir)
        size2pos2aa2sample2usage = aaProfilePerLen(sample2stats, l)
        printProfile(size2pos2aa2sample2usage, sampleNames, lendir)
    return

def getaaUsageMultiSamples(indir, outdir, sizes, uniq):
    samples = readfiles( indir, 1, 1 )
    sample2stats = {}
    for sample in samples:
        stats = len2aaSampleProfile(sample, sizes, uniq)
        sample2stats[sample.name] = stats
    getAaUsage(sample2stats, outdir)

############# TTESTS #########################
def getSamplePosUsage(aas, sample, aa2sample2usage):
    vec = []
    for aa in aas:
        if aa not in aa2sample2usage:
            raise ValueError('Unknown aa %s, known aas include: %s' %(aa, ','.join(aa2sample2usage.keys())))
        if sample not in aa2sample2usage[aa]:
            vec.append(0.0)
        else:
            vec.append(aa2sample2usage[aa][sample])
    return vec

def aausage_valTtestsPos(samples, aa2sample2usage, targetSamples):
    aas = sorted(aa2sample2usage.keys())
    aa2stats = {}

    sample2total = {}
    for sample2usage in aa2sample2usage.values():
        for sample, usage in sample2usage.iteritems():
            if sample not in sample2total:
                sample2total[sample] = usage
            else:
                sample2total[sample] += usage

    for aa in aas:
        targetVals = []
        otherVals = []
        sample2usage = aa2sample2usage[aa]
        for s in samples:
            if s not in sample2total or sample2total[s] == 0:#skip empty sample
                continue

            v = 0.0
            if s in sample2usage:
                v = sample2usage[s]
            if s in targetSamples:
                targetVals.append(v)
            else:
                otherVals.append(v)
        tval, pval = ttest_ind(targetVals, otherVals)
        aa2stats[aa] = [tval, pval, np.mean(targetVals), np.std(targetVals), np.mean(otherVals), np.std(otherVals)]
    return aa2stats

def aausage_corrTtestsPos(samples, aa2sample2usage, targetSamples):
    targetPairs = []
    otherPairs = []
    
    aas = sorted( aa2sample2usage.keys() )
    #targetGroup pairs:
    for i in xrange( len(targetSamples) - 1 ):
        s1 = targetSamples[i]
        vec1 = getSamplePosUsage(aas, s1, aa2sample2usage)
        if sum(vec1) == 0: #empty sample
            continue

        for j in xrange( i+1, len(targetSamples) ):
            s2 = targetSamples[j]
            vec2 = getSamplePosUsage(aas, s2, aa2sample2usage)
            if sum(vec2) == 0:
                continue

            corr, pval = pearsonr(vec1, vec2)
            targetPairs.append(corr)
    
    #Diff pairs:
    for s1 in samples:
        if s1 in targetSamples:
            continue
        vec1 = getSamplePosUsage(aas, s1, aa2sample2usage)
        if sum(vec1) == 0:
            continue
        for s2 in targetSamples:
            vec2 = getSamplePosUsage(aas, s2, aa2sample2usage)
            if sum(vec2) == 0:
                continue
            corr, pval = pearsonr(vec1, vec2)
            otherPairs.append(corr)

    #t-test:
    tval, pval = ttest_ind(targetPairs, otherPairs)
    return tval, pval, np.mean(targetPairs), np.std(targetPairs), np.mean(otherPairs), np.std(otherPairs)

def aausage_ttestsSize(type, outfile, samples, pos2aa2sample2usage, targetSamples, pvalcutoff):
    f = open(outfile, 'w')
    for pos, aa2sample2usage in pos2aa2sample2usage.iteritems():
        if type == 'corr':
            ttestStats = aausage_corrTtestsPos(samples, aa2sample2usage, targetSamples)
            if ttestStats[1] <= pvalcutoff:
                printTtestResults(f, str(pos), ttestStats)
        elif type == 'val':
            aa2ttestStats = aausage_valTtestsPos(samples, aa2sample2usage, targetSamples)
            for aa, ttestStats in aa2ttestStats.iteritems():
                if ttestStats[1] <= pvalcutoff:
                    name = '%s-%s' %(pos, aa)
                    printTtestResults(f, name, ttestStats)
        else:
            raise ValueError("Unknown ttest type: %s. Valid values are ['corr', 'val']" %type)
    f.close()

def aausage_ttests(outdir, sample2stats, targetSamples, pvalcutoff):
    lens = xrange(10, 20)
    sampleNames = sorted( sample2stats.keys() )

    for l in lens: #each length
        lendir = os.path.join(outdir, "%d" %l)
        size2pos2aa2sample2usage = aaProfilePerLen(sample2stats, l)
        for size, pos2aa2sample2usage in size2pos2aa2sample2usage.iteritems():
            sizedir = os.path.join(lendir, "%d" %size)
            system("mkdir -p %s" % sizedir)
            #Hypothesis: the amino acid usage at each position of the targetGroup's samples are more correlated with each other than with other samples
            corrttestfile = os.path.join(sizedir, 'corrTtests.txt')
            aausage_ttestsSize('corr', corrttestfile, sampleNames, pos2aa2sample2usage, targetSamples, pvalcutoff)

            #Hypothesis: at each position, certain amino acid is significantly overused/underused in targetGroup samples than in other samples
            valttestfile = os.path.join(sizedir, 'valTtests.txt')
            aausage_ttestsSize('val', valttestfile, sampleNames, pos2aa2sample2usage, targetSamples, pvalcutoff)

def printTtestResults(fh, name, stats):
    #fh: file handle, name=rowname (name of the ttest), stats=[tval, pval, sameMean, sameStd, diffMean, diffStd]
    #f.write("Genetype\tPval\tSameMean\tSameStd\tDiffMean\tDiffStd\n")
    fh.write("%s\t%s\n" %(name, '\t'.join([str(c) for c in stats[1:]])) )


















