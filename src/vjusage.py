#!/usr/bin/env python

#
#nknguyen at soe ucsc edu
#Aug 05 2011
#Reads in fasta file of cdr3 aa sequences
#The header of each sequence must be in the format:
#>id|v1[,v2,...]|j1[,j2,...];size=#
#Note: if there are tie calls of Vs or Js, the counts are split evenly between the ties.
#For example, if sequence S mapped to V1, V2, V3 and has 90 read counts, then 30 counts go to each V1, V2, V3
#Output: Tabs and plots of v usage, j usage, and vj usage
#

import sys, re, os, random, copy
from optparse import OptionParser
from scipy.stats.stats import pearsonr, spearmanr, kendalltau

from sonLib.bioio import system
import numpy as np
import immunoseq.lib.immunoseqLib as iseqlib 

class Sample:
    def __init__(self, name):
        self.name = name
        self.seqs = {} #key = header, val = Sequence()
        #self.seqs = iseqlib.Seqs()
    
    def getVJusage(self):
        self.usage = {}
        self.usage['v'], self.usage['j'], self.usage['vj'] = getSampleVJusage(self.seqs)

    def setCounts(self):
        #Get total read counts and total uniq sequence count
        self.totalReads = 0
        self.uniqueSeqs = 0
        for k, v in self.usage['vj'].iteritems():
            self.totalReads += v[0]
            self.uniqueSeqs += v[1]

class Sequence:
    def __init__(self, header, sequence):
        l = header.split(';')
        if(len(l) < 2):
            sys.stderr.write("Unvalid sequence header: %s. Sequence must end with ';size=#'\n" %header)
            exit(1)

        #self.header = header
        self.count = int( l[ len(l) -1 ].lstrip('size=') )
        
        clusters = l[1].split(',,')
        vs = []
        js = []
        for cluster in clusters:
            items = cluster.split('|')
            if len(items) < 3:
                sys.stderr.write("Sequence has wrong header format. Should be 'sampleName;id|v1[,v2,...]|j1[,j2,...];size=count'.\n'%s' was given\n" %header)
                exit(1)
            #self.id = items[0]
            currVs = items[1].split(',')
            currJs = items[2].split(',')
            for v in currVs:
                if v not in vs:
                    vs.append(v)
            for j in currJs:
                if j not in js:
                    js.append(j)
        self.v = ','.join( sorted(vs) )
        self.j = ','.join( sorted(js) )
        self.aa = sequence
        self.header = "|".join([self.aa, self.v, self.j])

    def __cmp__(self, other):
        return cmp(self.header, other.header)

    def setCount(self, count):
        self.count = count

    def updateCount(self, count):
        self.count += count

def addAvrSample( samples ):
    ''' Add the average and standardDev of all the samples '''
    if len(samples) == 0:
        return
    avrusage = {'v':{}, 'j':{}, 'vj':{}} #'v':{ 'vgene':[totalreads, uniqseqs] }
    stdusage = {'v':{}, 'j':{}, 'vj':{}} #'v':{ 'vgene':[totalreads, uniqseqs] }

    #get accumulate count across samples:
    for s in samples:
        for type in avrusage:
            g2c = s.usage[type]
            typeusage = avrusage[type]
            for g in g2c:
                if g not in typeusage:
                    typeusage[g] = [ g2c[g] ]
                else:
                    typeusage[g].append( g2c[g] )
                    #typeusage[g][1] += g2c[g][1]
    #average:
    avrsample = Sample('average')
    stdsample = Sample('std')
    for type in avrusage:
        for g in avrusage[type]:
            totalreads = [ sample[0] for sample in avrusage[type][g] ]
            uniqseqs = [ sample[1] for sample in avrusage[type][g] ]
            avrusage[type][g] = [np.mean(totalreads), np.mean(uniqseqs)]
            stdusage[type][g] = [np.std(totalreads), np.std(uniqseqs)]
            
    avrsample.usage = avrusage
    avrsample.setCounts()
    stdsample.usage = stdusage
    stdsample.setCounts()
    
    samples.append(avrsample)
    samples.append(stdsample)

def readFiles( indir ):
    samples = []
    for f in os.listdir( indir ):
        if os.path.isdir( os.path.join(indir, f) ):
            continue
        fh = open( os.path.join(indir, f), 'r')
        s = Sample( f.split('.')[0] )
        
        header = ""
        seq = ""
        for line in fh.readlines():
            if not re.search('\w', line):
                continue
            line = line.strip()
            if line[0] == '>': #header
                if header != "":
                    sequence = Sequence(header, seq)
                    s.seqs[sequence.header] = sequence
                    #s.seqs.add( Sequence(header, seq) )
                header = line.strip().lstrip('>')
                seq = ""
            else:
                seq += line
        fh.close()

        if header != "":
            sequence = Sequence(header, seq)
            s.seqs[sequence.header] = sequence
            #s.seqs.add( Sequence(header, seq) )
        samples.append(s)
    return samples


#============== SAMPLING =================
def samplingSample(sample, size):
    seqs = sample.seqs
    newseqs = {}
    indexList = []#index list represent all the sequences of the sample
    for i, s in enumerate(seqs.values()):
        indexList.extend( [i for j in xrange(s.count)] )

    sampleTotalReads = sum( [s.count for s in seqs.values()] )
    if sampleTotalReads < size:
        return 
    
    #sampling "size" indices from the indexList
    chosenIndices = random.sample( indexList, size )
    headers = seqs.keys()
    #for k in xrange(size):
    #    j = random.randint( 0, len(indexList) -1 ) #randomly pick one sequence
    #    i = indexList.pop(j) #remove the select sequence out of the indexlist
    for i in chosenIndices:
        header = headers[i]
        if header in newseqs:
            newseqs[header].updateCount(1)
        else:
            newseqs[header] = copy.copy( seqs[header] )
            newseqs[header].setCount(1)
    #Update sample with new sequences
    sample.seqs = newseqs
    return

def sampling(samples, size):
    '''Randomly pick "size" sequences from each sample'''
    for sample in samples:
        samplingSample(sample, size)
        sys.stderr.write("Done sampling sample %s\n" %sample.name)
    return

#================ V,J usage ============================
def sortDictByValue( dictionary ):
    items = [ (v, k) for k,v in dictionary.items() ]
    items.sort()
    items.reverse()
    return [ (k, v) for v, k in items ]

def getSampleVJusage(seqs):
    v2count = {} #key = TRBV geneName, val = [totalReads, uniqueSeqs]
    j2count = {} #key = TRBJ geneName, val = [totalReads, uniqueSeqs]
    vj2count = {} #key = 'TRBVgeneName,TRBJgeneName', val = [totalReads, uniqueSeqs]
    for seq in seqs.values():
        vs = seq.v.split(',')
        js = seq.j.split(',')
        vcount = seq.count/len(vs) #average the counts over all tied-genes
        jcount = seq.count/len(js)
        vjcount = seq.count/(len(vs)*len(js))
        for v in vs:
            if v not in v2count:
                v2count[v] = [vcount, 1.0/len(vs)]
            else:
                currvcount = v2count[v]
                v2count[v] = [ vcount + currvcount[0], 1.0/len(vs) + currvcount[1] ]

        for j in js:
            if j not in j2count:
                j2count[j] = [jcount, 1.0/len(js)]
            else:
                currjcount = j2count[j]
                j2count[j] = [ jcount + currjcount[0], 1.0/len(js) + currjcount[1] ]

        for v in vs:
            for j in js:
                vj = '|'.join([v, j])
                if vj not in vj2count:
                    vj2count[vj] = [vjcount, 1.0/(len(vs)*len(js))]
                else:
                    currvjcount = vj2count[vj]
                    vj2count[vj] = [ vjcount + currvjcount[0], 1.0/(len(vs)*len(js))+ currvjcount[1] ]
            
    return v2count, j2count, vj2count

def getUnionGeneList(samples, type):
    #Get the union of vgenes lists from all samples. 
    genes = []
    for s in samples:
        #print s.usage[type].keys()
        for g in s.usage[type].keys():
            if g not in genes:
                genes.append(g)
    #print genes
    #If a sample doesn't have a vgene, put the count of that vgene to 0
    genes.sort()
    for g in genes:
        for s in samples:
            if g not in s.usage[type].keys():
                s.usage[type][g] = [0,0]
    
    return genes

def getUsage(samples, outdir, type):
    genes = getUnionGeneList(samples, type)
    sys.stderr.write("Done getting uniqGeneList\n")
    #Print out usage table for each sample:
    for s in samples:
        g2c = s.usage[type]
        tabfile = os.path.join( outdir, "%s-%s.txt" %(s.name, type) )
        f = open( tabfile, 'w') 
        f.write("Gene\tTotal\tUniq\n")
        for g in genes:
            f.write( "%s\t%d\t%d\n" %(g, g2c[g][0], g2c[g][1]) )
        f.close()

def getVJusage(samples, outdir, abs, uniq):
    #If abs is True: print absolute count, otherwise print frequencies.
    #If uniq is True: using the Uniq sequence Count as the unit, otherwise, use read count
    if abs:
        if uniq:
            outdir = os.path.join(outdir, "absuniq")
        else:
            outdir = os.path.join(outdir, "abs")
    else:
        if uniq:
            outdir = os.path.join(outdir, "reluniq")
        else:
            outdir = os.path.join(outdir, "rel")
    system("mkdir -p %s" %outdir)
    

    for s in samples:
        v2c = s.usage['v']
        j2c = s.usage['j']
        vj2c = s.usage['vj']

        file = os.path.join( outdir, "%s" %s.name )
        if abs:
            file += "_abs"
        if uniq:
            file += "_uniq"
        file += "-vj.txt"

        f = open(file, 'w')
        f.write( "\t%s\n" %( '\t'.join( [j for j in sorted(j2c.keys())] ) ) )
        
        for v in sorted( v2c.keys() ):
            if v == '' or re.search('undefined', v):
                continue
            f.write( "%s" %v )
            for j in sorted( j2c.keys() ):
                vj = '|'.join([v, j])
                if vj not in vj2c:
                    f.write("\t0")
                else:
                    if uniq:#uniq seq count
                        count = vj2c[vj][1]
                    else:#read count
                        count = vj2c[vj][0]

                    if abs:
                        f.write("\t%d" %count)
                    else:#relative
                        if uniq:
                            count = float(count)/s.uniqueSeqs
                        else:
                            count = float(count)/s.totalReads
                        f.write("\t%f" %count)
            f.write("\n")
        f.close()

def topVJusage(samples, outdir, num):
    sample2topvjs = {}
    for sample in samples:
        if sample.name == 'average' or sample.name == 'std':
            continue
        vj2c = sample.usage['vj']
        sortedvjs = sorted( [(vj, counts[0]) for vj, counts in vj2c.iteritems()], key=lambda item:item[1], reverse=True )
        if num > len(sortedvjs):
            topvjs = sortedvjs
        else:
            topvjs = sortedvjs[:num]
        sample2topvjs[sample.name] = topvjs
    
    outfile = os.path.join(outdir, 'topvj.txt')
    f = open(outfile, 'w')
    names = sorted( sample2topvjs.keys() )
    f.write("%s\n" %('\t'.join(names)))
    for i in xrange(0, num):
        row = []
        for n in names:
            #print i, n
            #print sample2topvjs[n][i]
            vj = sample2topvjs[n][i][0]
            row.append( vj )
        f.write("%s\n" %('\t'.join(row) ))
    f.close()

def vjUsage(samples, outdir):
    #types = ['v', 'j', 'vj']
    types = ['v', 'j']
    for type in types:
        typeoutdir = os.path.join(outdir, type)
        system( "mkdir -p %s" %(typeoutdir) )
        getUsage(samples, typeoutdir, type)
        sys.stderr.write("Done %s usage\n" %type)
    return

def checkOptions(parser, args, options):
    if not options.indir:
        parser.error('Input directory was not specified\n')
    if not options.outdir:
        parser.error('Output directory was not specified\n')
    if not os.path.isdir(options.indir):
        parser.error('Input directory %s does not exist or is not a directory.\n' %options.indir)
    if not os.path.isdir(options.outdir):
        parser.error('Output directory %s does not exist or is not a directory.\n' %options.outdir)

    return

def initOptions(parser):
    parser.add_option('-i', '--indir', dest='indir', help='Required argument. Input directory')
    parser.add_option('-o', '--outdir', dest='outdir', help='Required argument. Output directory')
    parser.add_option('-s', '--sampling', dest='sampling', type='int', help='Sampling size, this number should be less than or equal to the size of the smallest sample. If this option is specified, randomly choose this number of sequences from each sample, and compare the V, J, VJ usage from these subsets. This approach, in away, normalize the data.')
    parser.add_option('-t', '--topseqs', dest='numtop', type='int', default=20, help='Number of top sequences to print out. Default = %default')
    #parser.add_option('-f', '--filter', dest='filter', action='store_true', default=False, help='If specified, filter out sequences that mapped unambiguously to multiple V or J genes. Default=%default')

def main():
    usage = "usage: %prog [options]\n"
    parser = OptionParser(usage = usage)
    initOptions(parser)
    options, args = parser.parse_args()
    checkOptions(parser, args, options)

    samples = readFiles( options.indir )
    sys.stderr.write("Done reading input files\n")

    if options.sampling:
        sampling(samples, options.sampling)
        sys.stderr.write("Done sampling\n")

    for sample in samples:
        sample.getVJusage()
        sample.setCounts()
        sys.stderr.write("Done getting usage for %s\n" %sample.name)
    
    #Adding the average of all samples the the sample list
    addAvrSample( samples )
    sys.stderr.write("Done adding average and std sample\n")

    vjUsage(samples, options.outdir)
    sys.stderr.write("Done v usage and j usage\n")

    vjoutdir = os.path.join( options.outdir, "vj")
    system("mkdir -p %s" %vjoutdir)
    #Generate VJ using the uniq sequence count or using the read count, relative or absolute count
    abs = True
    uniq = True
    getVJusage(samples, vjoutdir, abs, not uniq)
    sys.stderr.write("Done vj usage with absolute read count\n")
    getVJusage(samples, vjoutdir, not abs, not uniq)
    sys.stderr.write("Done vj usage with relative read count\n")
    getVJusage(samples, vjoutdir, abs, uniq)
    sys.stderr.write("Done vj usage with absolute uniqSeq count\n")
    getVJusage(samples, vjoutdir, not abs, uniq)
    sys.stderr.write("Done vj usage with relative read count\n")

    #TopVJ usage across samples
    #num = 20
    topVJusage(samples, vjoutdir, options.numtop)

if __name__ == '__main__':
    main()

