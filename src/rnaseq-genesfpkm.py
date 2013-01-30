#!/usr/bin/env python

#Jun 11 2012
#
#Input: input directory containing genes.fpkm_tracking of each sample
#       file that maps geneID to name 
#       (list of genes of interest)
#Output: FPKM table where rows = genes and columns = samples
#

import os, sys, re, time
from optparse import OptionParser
from sonLib.bioio import system
from immunoseq.lib.immunoseqLib import *

class SampleGene():
    def __init__(self, fpkm, low, high, status):
        self.fpkm = fpkm
        self.fpkmlow = low
        self.fpkmhigh = high
        self.status = status

class Gene():
    def __init__(self, line):
        #tracking_id    class_code  nearest_ref_id  gene_id gene_short_name tss_id  locus   length  coverage    FPKM    FPKM_conf_lo    FPKM_conf_hi    FPKM_status
        #Q8NH21  -   -   Q8NH21  Q8NH21  -   chr1:69090-70008    -   -   4.94066e-324    4.94066e-324    4.94066e-324    OK
        items = line.strip().split('\t')
        if len(items) < 13:
            raise ValueError('Gene record has less than 13 (only %d) fields:\n%s\n' %(len(items), line))
        self.id = items[4]
        self.name = items[5]
        self.geneid = items[3]
        self.locus = items[6]
        self.samples = [] 
        self.sumfpkm = 0.0
        for i in xrange(9, len(items) - 3, 4):
            fpkm = float(items[i])
            fpkmlow = float(items[i+1])
            fpkmhigh = float(items[i+2])
            status = items[i+3]
            sample = SampleGene( fpkm, fpkmlow, fpkmhigh, status)
            self.samples.append(sample)
            self.sumfpkm += fpkm

def getSamples(header):
    fields = header.split('\t')
    if len(fields) < 13:
        raise ValueError('Gene fpkm_tracking header has less than 13 (only %d) fields:\n%s\n' %(len(fields), header))

    samples = []
    for i in xrange( 9, len(fields), 4):
        name = fields[i].rstrip('_FPKM')
        samples.append(name)
    return samples

def readGenelist(file):
    genes = {}
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        genes[line] = 1
    f.close()
    return genes

def readId2name(file):
    id2name = {}
    f = open(file, 'r')
    for line in f:
        items = line.strip().split('\t')
        id = items[0]
        name = items[1]
        if id not in id2name:
            id2name[id] = name
        else:
            raise ValueError('Redundant records for id %s in file %s\n' %(id, file))
    f.close()
    return id2name

def readGeneFile(file, cutoff):
    id2gene = {} 
    f = open(file, 'r')
    samplenames = []
    for line in f:
        if re.match('tracking', line):
            samplenames = getSamples( line.strip() )
            continue
        gene = Gene(line)

        flag = False
        for s in gene.samples:
            if s.fpkm > cutoff:
                flag = True
                break
        if not flag:
            continue

        #if gene.id not in id2gene and gene.sumfpkm > 100:
        if gene.id not in id2gene:
            id2gene[gene.id] = gene
        #else:
        #    raise ValueError('Gene %s has redundant records in file %s\n' %(gene.id, file))
    f.close()
    return id2gene, samplenames

def printFPKMmatrix( samples, id2gene, id2name, genes, outfile, logtransform ):
    #get a look up table of sample to the index in gene.samples
    sample2index = {}
    for i, sample in enumerate( samples ):
        sample2index[sample] = i

    samples = sorted(samples)
    selectedsamples=['P1', 'P10', 'P11', 'P12', 'P13', 'ControlB']
    #selectedsamples=['P01', 'P10', 'P11', 'P12', 'P13', 'PB']
    #selectedsamples=['P01', 'P10', 'P11', 'P12', 'P13', 'PB', 'P18', 'P19', 'P20']

    matrix = [] #2D array where rows = genes and columns = samples
    nonzeromin = float('inf')
    for id, gene in id2gene.iteritems():
        name = id
        if id2name:
            iditems = id.split(',')
            for iditem in iditems:
                if iditem in id2name:
                    name = id2name[iditem]
                    break
        #elif gene.name != '-':
        #    name = gene.name

        if genes:
            nameitems = name.split(',')
            nameitems = [ item.split('.')[0] for item in nameitems ]
            flag = False
            for item in nameitems:
                if item in genes:
                    flag = True
                    break
            if not flag:
               continue
        
        rowVals = [name]
        #for sample in samples:
        for sample in selectedsamples:
            index = sample2index[sample]
            samplegene = gene.samples[index]
            fpkm = samplegene.fpkm
            rowVals.append(fpkm)
            if fpkm > 0:
                nonzeromin = min([nonzeromin, fpkm])
        matrix.append(rowVals)

    if logtransform:
        for row in matrix:
            for i, c in enumerate(row):
                if i == 0:
                    continue
                if c == 0:
                    row[i] = log10(nonzeromin*0.5)
                else:
                    row[i] = log10(c)
    
    #Write to output file:
    f = open(outfile, 'w')
    #f.write("Gene\t%s\n" %("\t".join(samples)))
    f.write("Gene\t%s\n" %("\t".join(selectedsamples)))
    for row in matrix: 
        f.write("%s" % row[0])
        for v in row[1:]:
            f.write("\t%f" %v)
        f.write("\n")
    f.close()

def main():
    parser = initOptions()
    #parser.add_option('-i', '--indir', dest='indir', default='.', help='Input directory. Default=%default') 
    parser.add_option('-i', '--infile', dest='infile', help='Input fpkm_trackig file. Required argument.') 
    parser.add_option('-o', '--outdir', dest='outdir', default='outputs', help='Output directory. Default=%default')
    parser.add_option('-n', '--name', dest='name', help='File mapping gene/tx id to name.')
    parser.add_option('-g', '--genes', dest='genes', help='File containing the gene list of interest')
    parser.add_option('-l', '--logtransform', dest='logtransform', action='store_true', default=False, help='If specified, perform a log10 transform on the data')
    parser.add_option('-m', '--minfpkm', dest='minfpkm', type='float', default=0.0, help='Only include the transcript in output if at least one sample has value >= this cutoff. Defalut=%default')

    options, args = parser.parse_args()
    #checkOptions(parser, args, options)

    id2gene, samples = readGeneFile( options.infile, options.minfpkm )
    
    if options.name:
        options.name = readId2name(options.name)     
    if options.genes:
        options.genes = readGenelist(options.genes)

    system('mkdir -p %s' %options.outdir)
    outfile = os.path.join(options.outdir, 'fpkmMatrix.txt')
    printFPKMmatrix( samples, id2gene, options.name, options.genes, outfile, options.logtransform )


if __name__ == '__main__':
    main()
