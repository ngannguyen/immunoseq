#nknguyen soe ucsc edu
#Tue Jul 17 12:41:31 PDT 2012
#Library of functions used to compute sequence length distribution

import os, sys, re
import matplotlib.pyplot as pyplot
from immunoseq.lib.immunoseqLib import *
import numpy as np
from scipy.stats.stats import pearsonr
from scipy.stats import ttest_ind


def getLenDist(seqs):
    len2freq = {} #key = len, freq = [readFreq, uniqFreq] where readFreq is %of total read with this len, and uniqFreq is percentage of total clones with this len
    
    for s in seqs.values():
        l = len(s.seq)
        if l not in len2freq:
            len2freq[l] = [s.count, 1]
        else:
            counts = len2freq[l]
            len2freq[l] = [counts[0] + s.count, counts[1] + 1]
    
    totalUniqs = len(seqs)
    totalReads = sum([s.count for s in seqs.values()])
    if totalUniqs > 0 and totalReads >0:
        for l, counts in len2freq.iteritems():
            len2freq[l] = [ counts[0]*100.0/totalReads, counts[1]*100.0/totalUniqs ]
    return len2freq

def addSamplingLenDist(len2freq, agglen2freq, i):
    for l, f in len2freq.iteritems():
        if l not in agglen2freq:
            agglen2freq[l] = [ [0.0]*i + [f[0]], [0.0]*i + [f[1]] ]
        else:
            agglen2freq[l][0].append(f[0])
            agglen2freq[l][1].append(f[1])

def avrSamplingLenDist(agglen2freq):
    l2mean = {}
    l2std = {}
    for l, f in agglen2freq.iteritems():
        l2mean[l] = [np.mean(f[0]), np.mean(f[1])]
        l2std[l] = [np.std(f[0]), np.std(f[1])]
    return l2mean, l2std

def getSampleLenDist(samples):
    sample2len2freq = {}
    lens = []
    for sample in samples:
        #len2freq = getLenDist( sample.seqs )
        len2freq = getLenDist( sample.seqs )
        sample2len2freq[sample.name] = len2freq
        for l in len2freq.keys():
            if l not in lens:
                lens.append(l)

    for s, l2f in sample2len2freq.iteritems():
        for l in lens:
            if l not in l2f:
                l2f[l] = [0.0, 0.0]
    return sample2len2freq, sorted(lens)

def getUnionLenDist(sample2len2freq):
    lens = []
    for len2freq in sample2len2freq.values():
        for l in len2freq.keys():
            if l not in lens:
                lens.append(l)
    for s, l2f in sample2len2freq.iteritems():
        for l in lens:
            if l not in l2f:
                l2f[l] = [0.0, 0.0]
    return sorted(lens)

def printLenDistTabSample(len2avr, len2std, outfile):
    f = open(outfile, 'w')
    f.write("Length\tUniqMean\tReadMean\tUniqStd\tReadStd\n")
    for l in sorted(len2avr.keys()):
        avr = len2avr[l]
        std = len2std[l]
        f.write("%d\t%f\t%f\t%f\t%f\n" %(l, avr[1], avr[0], std[1], std[0]))
    f.close()

def printLenDistTab(sample2avr, sample2std, outdir):
    #lens = getUnionLenDist(sample2avr)
    #getUnionLenDist(sample2std)
    for sample, len2avr in sample2avr.iteritems():
        outfile = os.path.join(outdir, '%s.txt' %sample)
        len2std = sample2std[sample]
        printLenDistTabSample(len2avr, len2std, outfile)

def drawLenDistData(axes, sample2lendist, uniq):
    #if uniq is true, return the distribution of sequence length versus percentage of total uniq sequences (clones)
    #if uniq is false, the yaxis is percentage of total reads instead
    #sample2lendist, xdata = getSampleLenDist( samples )
    #xdata = getUnionLenDist(sample2lendist)
    xdata = sorted( sample2lendist.values()[0].keys() )
    if not xdata:
        return

    lines = []
    labels = []
    colors = getColors6()
    #HACK:
    if len(sample2lendist) > len(colors):
        colors.extend( getColors6light() )
    if len(sample2lendist) > len(colors):
        colors.extend( getColors6dark() )

    #HACK color as vs healthy:
    #colors = ["#E31A1C", "#377EB8"] #red, blue
    #patients = ['as1D', 'as11D', 'as8D', 'as15D', 'as16D', 'patients']
    s2c, s2m, s2cLight = sample2color(colors) 
    #END HACK color

    #lightColors = getColors6light()
    markersize = 10.0
    c = -1
    name2line = {}
    for s in sorted( sample2lendist.keys() ):
        len2freq = sample2lendist[s]
        c += 1
        ydata = []
        for l in xdata:
            if uniq:
                ydata.append( len2freq[l][1] )
            else:
                ydata.append( len2freq[l][0] )
        
        #color = colors[c]
        #HACK COLOR
        #color = colors[1]
        #if s in patients:
        #    color = colors[0]
        #END HACK COLOR

        line = axes.plot(xdata, ydata, color=s2c[s], marker=s2m[s], markeredgecolor=s2c[s], markersize = markersize, linestyle='-', linewidth=2)
        #line = axes.plot(xdata, ydata, color=color, marker='o', markeredgecolor=color, markersize = markersize, linestyle='-', linewidth=2)
        #axes.plot(xdata, ydata, color=lightColors[c], linestyle='-', linewidth=0.5)
        lines.append(line)
        labels.append( s)
        name2line[s] = line
    
    axes.set_xlim(8.5, 23.5)
    xticks = xrange(8, 24, 2)
    xticklabels = [ str(x) for x in xticks]
    axes.xaxis.set_ticks(xticks)
    axes.xaxis.set_ticklabels( xticklabels )

    editSpine( axes )
    axes.set_xlabel("Length (amino acids)", size='x-large')
    if uniq:
        axes.set_title('CDR3 Unique Sequence Length Distribution', size="xx-large")
        axes.set_ylabel("Frequency (% of total clones)", size='x-large')
    else:
        axes.set_title('CDR3 Sequence Length Distribution', size="xx-large")
        axes.set_ylabel("Frequency (% of total sequences)", size='x-large')
    
    #HACK: Sort labels:
    labelorder = ['as15D', 'as20D', 'asBD', 'as8D', 'as16D', 'as1D', 'as11D']
    currlabels = []
    currlines = []
    for label in labelorder:
        if label in labels:
            currlabels.append( properName(label) )
            currlines.append(name2line[label])
    legend = axes.legend( currlines, currlabels, numpoints=1, loc="best", ncol=1 )
    #END HACK

    #legend = axes.legend( lines, labels, numpoints=1, loc='best', ncol=1)
    legend.__drawFrame = False
    axes.yaxis.grid(b=True, color="#BDBDBD", linestyle='-', linewidth=0.005)
    axes.xaxis.grid(b=True, color="#BDBDBD", linestyle='-', linewidth=0.005)

def drawLenDist(samples, outfile, uniq):
    dpi = 300
    outformat = 'pdf'
    fig, pdf = initImage2(10.0, 10.0, outformat, outfile, dpi)
    axes = setAxes(fig)
    sample2lendist, xdata = getSampleLenDist( samples )
    drawLenDistData(axes, sample2lendist, uniq)
    writeImage2(fig, pdf, outformat, outfile, dpi)

def drawLenDist2(sample2len2freq, outfile, uniq):
    dpi = 300
    outformat = 'pdf'
    fig, pdf = initImage2(10.0, 10.0, outformat, outfile, dpi)
    axes = setAxes(fig)
    drawLenDistData(axes, sample2len2freq, uniq)
    writeImage2(fig, pdf, outformat, outfile, dpi)

###### TTEST #######
def getFreqVec(lendist, uniq):
    vec = []
    for l in sorted(lendist.keys()):
        freqs = lendist[l]
        if uniq:
            vec.append(freqs[1])
        else:
            vec.append(freqs[0])
    return vec

def lenDistTtest(sample2len2freq, group2samples, uniq, targetGroup):
    samePairs = []
    diffPairs = []

    #same pairs
    for g, gsamples in group2samples.iteritems():
        issame = True
        if targetGroup and g != targetGroup:
            issame = False
        for i in xrange( len(gsamples) - 1 ):
            if gsamples[i] not in sample2len2freq:
                raise KeyError("Could not found sample %s in sample2len2freq\n" %gsamples[i])
            lendist1 = sample2len2freq[ gsamples[i] ]
            vec1 = getFreqVec(lendist1, uniq)
            if sum(vec1) == 0:
                continue

            for j in xrange( i+1, len(gsamples) ):
                if gsamples[j] not in sample2len2freq:
                    raise KeyError("Could not found sample %s in the sample2len2freq\n" %gsamples[j])
                lendist2 = sample2len2freq[ gsamples[j] ]
                vec2 = getFreqVec(lendist2, uniq)
                if sum(vec2) == 0:
                    continue

                corr, pval = pearsonr(vec1, vec2)
                if issame:
                    samePairs.append( corr )
                #else:
                #    diffPairs.append( corr )
    #Diff pairs:
    groups = group2samples.keys()
    for i in xrange( len(groups) - 1 ):
        for gs1 in group2samples[ groups[i] ]:
            lendist1 = sample2len2freq[ gs1 ]
            vec1 = getFreqVec(lendist1, uniq)

            if sum(vec1) == 0: #skip empty sample
                continue
            
            for j in xrange( i+1, len(groups) ):
                for gs2 in group2samples[ groups[j] ]:
                    lendist2 = sample2len2freq[ gs2 ]
                    vec2 = getFreqVec(lendist2, uniq)
                    if sum(vec2) == 0:
                        continue

                    corr, pval = pearsonr(vec1, vec2)
                    diffPairs.append( corr )

    #t-test:
    if len(samePairs) == 0 or len(diffPairs) == 0:
        raise ValueError ('lenDistTtest, one of the vector has zero length.\nVec1: %s.\nVec2: %s\n' %(','.join(vec1), ','.join(vec2)))
    tval, pval = ttest_ind(samePairs, diffPairs)
    return tval, pval, np.mean(samePairs), np.std(samePairs), np.mean(diffPairs), np.std(diffPairs)

def getLens(sample2len2freq):
    lens = []
    for sample, l2f in sample2len2freq.iteritems():
        for l in l2f.keys():
            if l not in lens:
                lens.append(l)
    return lens

def lenDistTtest_singleLength(sample2len2freq, targetSamples, uniq):
    lens = getLens(sample2len2freq)
    len2stats = {}

    emptySamples = []
    for sample, l2f in sample2len2freq.iteritems():
        if not l2f or sum([f[0] for f in l2f.values()]) == 0:
            emptySamples.append(sample)

    for l in lens:
        targetVals = []
        otherVals = []
        for sample, l2f in sample2len2freq.iteritems():
            if sample in emptySamples:#Skip empty samples
                continue

            val = 0.0
            if l in l2f:
                freqs = l2f[l]
                if uniq:
                    val = freqs[1]
                else:
                    val = freqs[0]
            if sample in targetSamples:
                targetVals.append(val)
            else:
                otherVals.append(val)
        tval, pval = ttest_ind(targetVals, otherVals)
        len2stats[l] = [tval, pval, np.mean(targetVals), np.std(targetVals), np.mean(otherVals), np.std(otherVals)]
    return len2stats

def lenDistTtests(sample2len2freq, group2samples, outfile, targetGroup, pvalcutoff):
    uniq = True
    f = open(outfile, 'w')
    f.write("Type\tPval\tSameMean\tSameStd\tDiffMean\tDiffStd\n")
    #uniq = True
    tval, pval, sameMean, sameStd, diffMean, diffStd = lenDistTtest(sample2len2freq, group2samples, uniq, targetGroup)
    f.write("Uniq\t%f\t%f\t%f\t%f\t%f\n" %(pval, sameMean, sameStd, diffMean, diffStd) )
    tval, pval, sameMean, sameStd, diffMean, diffStd = lenDistTtest(sample2len2freq, group2samples, not uniq, targetGroup)
    f.write("Reads\t%f\t%f\t%f\t%f\t%f\n" %(pval, sameMean, sameStd, diffMean, diffStd) )
    
    #ttest for each length:
    if targetGroup and targetGroup in group2samples:
        targetSamples = group2samples[targetGroup]
        f.write("###\n")
        f.write("###\n")
        f.write("#Uniq\n")
        len2results = lenDistTtest_singleLength(sample2len2freq, targetSamples, uniq)
        for l, results in len2results.iteritems():
            if results[1] <= pvalcutoff:
                f.write("%d\t%s\n" %(l, '\t'.join([str(c) for c in results[1:]])))
        
        f.write("###\n")
        f.write("###\n")
        f.write("#Reads\n")
        len2results = lenDistTtest_singleLength(sample2len2freq, targetSamples, not uniq)
        for l, results in len2results.iteritems():
            if results[1] <= pvalcutoff:
                f.write("%d\t%s\n" %(l, '\t'.join([str(c) for c in results[1:]])))
    
    f.close()














