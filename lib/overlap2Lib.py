import os, sys, re, copy
import matplotlib.pyplot as pyplot
import immunoseq.lib.immunoseqLib as iseqlib
from sonLib.bioio import system

def getSubSamples(seqs):
    samples = []
    for header, name2seq in seqs.iteritems():
        for name in name2seq.keys():
            if name not in samples:
                samples.append(name)
    return samples

def sharedSeqFreq(sample):
    numsam2header2freqs = {}
    for header, name2seq in sample.seqs.iteritems():
        numsam = len( name2seq.keys() )
        #freqs = sorted( [ seq.freq for seq in name2seq.values() ] )
        n2freq = {}
        for n, s in name2seq.iteritems():
            #n2freq[n] = s.freq
            n2freq[n] = s.count #HACK

        if numsam not in numsam2header2freqs:
            numsam2header2freqs[numsam] = {header:n2freq}
        else:
            numsam2header2freqs[numsam][header] = n2freq
    return numsam2header2freqs

def printSharedSeqFreq(sample, minsam, outfile):
    f = open(outfile, 'w')
    numsam2header2freqs = sharedSeqFreq(sample)
    for numsam, header2freqs in numsam2header2freqs.iteritems():
        if numsam >= minsam:
            for header, n2f in header2freqs.iteritems():
                #freqs = n2f.values()
                #f.write("%s\t%s\n" %(header, ','.join( [str(freq) for freq in freqs] )) )
                names = []
                freqs = []
                for n, freq in n2f.iteritems():
                    names.append(n)
                    freqs.append( str(freq) )
                f.write("%s\t%s\t%s\n" % (header, ','.join(names) , ','.join(freqs) ))

    f.close()
    return

def printSharedSeqFreqAll( samples, minsam, outdir ):
    for s in samples:
        file = os.path.join( outdir, "%s-freqs.txt" %s.name )
        printSharedSeqFreq(s, minsam, file)

def filterByNumSample(sample, minsam, outfile, minfreq):
    #sample = calcFreq(sample) #HACK
    f = open(outfile, 'w')
    for header, name2seq in sample.seqs.iteritems():
        if len(name2seq.keys()) >= minsam:
            for name, seq in name2seq.iteritems():
                if seq.freq < minfreq:
                    continue
                h = seq.getFastaHeader()
                f.write("%s\n" %h)
                f.write("%s\n" %seq.seq)
    f.close()
    return

def filterByNumSampleAll(samples, minsam, outdir, minfreq):
    for s in samples:
        outfile = os.path.join(outdir, s.name)
        filterByNumSample(s, minsam, outfile, minfreq)

#========== DRAWING ===
def combineSamples(samples, name):
    #combine iseqlib samples into One sample where seqs = {header: {subSampleName: Seq}}
    sample = iseqlib.Sample(name)
    seqs = {}
    for s in samples:
        for header, seq in s.seqs.iteritems():
            if header not in seqs:
                seqs[header] = {s.name: seq}
            else:
                seqs[header][s.name] = seq
    sample.seqs = seqs
    return sample

def getSharedSeqDist(samples, uniq):
    sample2dist = {}
    for sample in samples:
        #sample = calcFreq(sample) #HACK
        numsam2count = {}
        for header, name2seq in sample.seqs.iteritems():
            numsam = len( name2seq.keys() ) #number of samples having this sequence
            if uniq:
                if numsam not in numsam2count:
                    numsam2count[numsam] = 1
                else:
                    numsam2count[numsam] += 1
            else:
                reads = sum([seq.count for seq in name2seq.values()])
                if numsam not in numsam2count:
                    numsam2count[numsam] = reads
                else:
                    numsam2count[numsam] += reads
        sample2dist[sample.name] = numsam2count
    return sample2dist

def printNumsamVsClonesDist(outfile, sample2dist, relative):
    #Get the union list of numsam:
    nums = []
    for s, num2count in sample2dist.iteritems():
        for num in num2count:
            if num not in nums:
                nums.append(num)
    nums.sort()

    #Convert to percentage if option 'relative' is True
    #if relative:
    #    for s, n2c in sample2dist.iteritems():
    #        total = sum(n2c.values())
    #        if total >0:
    #            for n, c in n2c.iteritems():
    #                n2c[n] = c*100.0/total

    f = open(outfile, 'w')
    f.write("Group\t%s\n" %( '\t'.join([str(n) for n in nums]) ))
    for s, num2count in sample2dist.iteritems():
        f.write("%s" %s)
        total = sum(num2count.values())
        for n in nums:
            if n not in num2count or total == 0:
                f.write("\t0")
            else:
                if not relative:
                    f.write("\t%f" %num2count[n])
                else:
                    f.write("\t%f" %( num2count[n]*100.0/total ))
        f.write("\n")
    f.close()

def drawDistData(axes, sample2dist):
    #sample2dist = getSharedSeqDist(samples, uniq)
    
    lines = []
    labels = []
    colors = iseqlib.getColors6()
    #lightColors = getColors6light()
    markersize = 10.0
    c = -1
    xmax = 0
    for s in sorted( sample2dist.keys() ):
        numsam2count = sample2dist[s]
        
        c += 1
        xdata = sorted( numsam2count.keys() )
        xmax = max([xmax, max(xdata)])
        ydata = [ numsam2count[x] for x in xdata ]
        totaly = sum(ydata)
        pcydata = [(100.0*y)/totaly for y in ydata]
        
        #line = axes.plot(xdata, ydata, color=colors[c], marker='o', markeredgecolor=colors[c], markersize = markersize, linestyle='-', linewidth=2)
        line = axes.plot(xdata, pcydata, color=colors[c], marker='o', markeredgecolor=colors[c], markersize = markersize, linestyle='-', linewidth=2)
        #axes.plot(xdata, ydata, color=lightColors[c], linestyle='-', linewidth=0.5)
        lines.append(line)
        labels.append(s)
        print s
        print xdata
        print ydata
        print pcydata
    
    axes.set_yscale('log')
    axes.set_xlim(0.8, xmax + 0.2)
    xticks = xrange(1, xmax + 1)
    xticklabels = [ str(x) for x in xticks ]
    axes.xaxis.set_ticks(xticks)
    axes.xaxis.set_ticklabels( xticklabels )

    axes.set_title('Shared sequences', size="xx-large")
    iseqlib.editSpine( axes )
    axes.set_xlabel("Number of samples", size='x-large')
    axes.set_ylabel("Number of clones", size='x-large')
    legend = axes.legend( lines, labels, numpoints=1, loc='best', ncol=1)
    legend.__drawFrame = False
    axes.yaxis.grid(b=True, color="#BDBDBD", linestyle='-', linewidth=0.005)
    axes.xaxis.grid(b=True, color="#BDBDBD", linestyle='-', linewidth=0.005)

def drawDist(samples, outfile, uniq):
    dpi = 300
    outformat = 'pdf'
    fig, pdf = iseqlib.initImage2(10.0, 10.0, outformat, outfile, dpi)
    axes = iseqlib.setAxes(fig)
    sample2dist = getSharedSeqDist(samples, uniq)
    drawDistData(axes, sample2dist)
    iseqlib.writeImage2(fig, pdf, outformat, outfile, dpi)

def drawNumsamVsClonesDist2(sample2dist, outfile):
    dpi = 300
    outformat = 'pdf'
    fig, pdf = iseqlib.initImage2(10.0, 10.0, outformat, outfile, dpi)
    axes = iseqlib.setAxes(fig)
    drawDistData(axes, sample2dist)
    iseqlib.writeImage2(fig, pdf, outformat, outfile, dpi)

######### CLONE MATRIX ###############
def calcFreq(sample):
    newsample = copy.copy(sample)
    name2total = {}
    for header, name2seq in sample.seqs.iteritems():
        for name, seq in name2seq.iteritems():
            if name in name2total:
                name2total[name] += seq.count
            else:
                name2total[name] = seq.count
    for header, name2seq in sample.seqs.iteritems():
        for name, seq in name2seq.iteritems():
            seq.setFreq( name2total[name] )
    return newsample

def printCloneMatrix(outfile, sample, minsam, freqTransform):
    sample = calcFreq(sample)
    f = open(outfile, 'w')
    f.write( "Clone\t%s\n" % '\t'.join(sample.subsamples) )
    for header, name2seq in sample.seqs.iteritems():
        rowvals = []
        if len(name2seq) < minsam:
            continue
        for name in sample.subsamples:
            count = 0
            if name in name2seq:
                count = name2seq[name].count
                if freqTransform:
                    count = int ( name2seq[name].freq*freqTransform + 0.5 )
            rowvals.append( str(count) )
        f.write("%s\t%s\n" %(header, "\t".join(rowvals)))
    f.close()

def printCloneMatrixAll(outdir, samples, minsam, freqTransform):
    for sample in samples:
        outfile = os.path.join(outdir, "cloneM-%s.txt" %sample.name)
        printCloneMatrix(outfile, sample, minsam, freqTransform)

#def main():
#    samples = readfiles(options.indir, options.count)
#    if options.drawdist:
#        drawDist(samples, options) 
#    #printSharedSeqFreqAll( samples, minsam, "counts-atleast3sams" )
#    if options.fasta:
#        faoutdir = os.path.join(options.outdir, "fasta-atleast%dsams" %options.minsam)
#        system("mkdir -p %s" %faoutdir)
#        filterByNumSampleAll(samples, options.minsam, faoutdir, options.freq)
#    if options.clonematrix:
#        printCloneMatrixAll(options.outdir, samples, options.minsam, options.freqTransform)


