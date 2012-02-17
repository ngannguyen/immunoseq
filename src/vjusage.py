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

import sys, re, os
from optparse import OptionParser
from scipy.stats.stats import pearsonr, spearmanr, kendalltau

from sonLib.bioio import system

class Sample:
    def __init__(self, name):
        self.name = name
        self.seqs = []
        self.totalReads = 0
        self.uniqueSeqs = 0
    
    def getVJusage(self):
        self.usage = {}
        self.usage['v'], self.usage['j'], self.usage['vj'] = getSampleVJusage(self.seqs)

class Sequence:
    def __init__(self, header, sequence):
        l = header.split(';')
        if(len(l) < 2):
            sys.stderr.write("Unvalid sequence header: %s. Sequence must end with ';size=#'\n" %header)
            exit(1)

        self.count = int( l[ len(l) -1 ].lstrip('size=') )
        items = l[0].split('|')
        if len(items) < 3:
            sys.stderr.write("Sequence has wrong header format. Should be 'id|v1[,v2,...]|j1[,j2,...];size=count'.\n'%s' was given\n" %header)
            exit(1)
        self.id = items[0]

        vs = sorted( items[1].split(',') )
        self.v = ','.join(vs)

        js = sorted( items[2].split(',') )
        self.j = ','.join(js)
        
        self.aa = sequence

def addAvrSample( samples ):
    ''' Add the average of all the samples '''
    if len(samples) == 0:
        return
    avrsample = Sample('average')
    avrusage = {'v':{}, 'j':{}, 'vj':{}} #'v':{ 'vgene':[totalreads, uniqseqs] }

    #get accumulate count across samples:
    for s in samples:
        for type in avrusage:
            g2c = s.usage[type]
            typeusage = avrusage[type]
            for g in g2c:
                if g not in typeusage:
                    typeusage[g] = g2c[g]
                else:
                    typeusage[g][0] += g2c[g][0] 
                    typeusage[g][1] += g2c[g][1]
    #average:
    for type in avrusage:
        for g in avrusage[type]:
            avrusage[type][g][0] /= len(samples)
            avrusage[type][g][1] /= len(samples)
    avrsample.usage = avrusage
    samples.append(avrsample)

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
                    s.seqs.append( Sequence(header, seq) )
                header = line.strip().lstrip('>')
                seq = ""
            else:
                seq += line
        fh.close()

        if header != "":
            s.seqs.append( Sequence(header, seq) )
        samples.append(s)
    return samples

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
    for seq in seqs:
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
        for g in s.usage[type].keys():
            if g not in genes:
                genes.append(g)

    #If a sample doesn't have a vgene, put the count of that vgene to 0
    genes.sort()
    for g in genes:
        for s in samples:
            if g not in s.usage[type].keys():
                s.usage[type][g] = [0,0]
    
    return genes

def getUsage(samples, outdir, type):
    genes = getUnionGeneList(samples, type)
    #Print out usage table for each sample:
    for s in samples:
        g2c = s.usage[type]
        tabfile = os.path.join( outdir, "%s-%s.txt" %(s.name, type) )
        f = open( tabfile, 'w') 
        for g in genes:
            f.write( "%s\t%d\t%d\n" %(g, g2c[g][0], g2c[g][1]) )
        f.close()

def getVJusage(samples, outdir):
    for s in samples:
        v2c = s.usage['v']
        j2c = s.usage['j']
        vj2c = s.usage['vj']

        file = os.path.join( outdir, "%s-vj.txt" %s.name )
        f = open(file, 'w')
        f.write( "\t%s\n" %( '\t'.join( [j for j in sorted(j2c.keys())] ) ) )
        for v in sorted( v2c.keys() ):
            f.write( "%s" %v )
            for j in sorted( j2c.keys() ):
                vj = '|'.join([v, j])
                if vj not in vj2c:
                    f.write("\t0")
                else:
                    f.write( "\t%d" %( vj2c[vj][0] ) )
            f.write("\n")
        f.close()

def vjUsage(samples, outdir):
    #types = ['v', 'j', 'vj']
    types = ['v', 'j']
    for type in types:
        typeoutdir = os.path.join(outdir, type)
        system( "mkdir -p %s" %(typeoutdir) )
        getUsage(samples, typeoutdir, type)
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

def main():
    usage = "usage: %prog [options]\n"
    parser = OptionParser(usage = usage)
    initOptions(parser)
    options, args = parser.parse_args()
    checkOptions(parser, args, options)

    samples = readFiles( options.indir )

    for sample in samples:
        sample.getVJusage()
    
    #Adding the average of all samples the the sample list
    addAvrSample( samples )
    vjUsage(samples, options.outdir)
    
    vjoutdir = os.path.join( options.outdir, "vj")
    system("mkdir -p %s" %vjoutdir)
    getVJusage(samples, vjoutdir)

if __name__ == '__main__':
    main()

