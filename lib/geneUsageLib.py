#nknguyen soe ucsc edu
#Tue Jul 17 10:56:47 PDT 2012
#Library of functions used to compute the gene usage

import sys, re, os, random, copy
from optparse import OptionParser
from scipy.stats.stats import pearsonr, spearmanr, kendalltau
from sonLib.bioio import system
import numpy as np
import immunoseq.lib.immunoseqLib as iseqlib 

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

def getGenes(seq, type):
    if type not in ['v', 'j', 'd']:
        raise ValueError("singleUsage, %s is not a valid genetype. Valid choices are v, d, j" %type)
    if type == 'v':
        return seq.vs
    elif type == 'j':
        return seq.js
    else:
        return seq.ds

def singleUsage(seqs, type):
    gene2count = {} #key = genename, val = [totalreads, uniqseqs]

    for seq in seqs.values():
        genes = getGenes(seq, type)
        
        #filter out unvalid genes:
        if len(genes) == 0 or '(undefined)' in genes or '' in genes:
            continue
        count = float(seq.count)/len(genes)
        for gene in genes:
            if gene not in gene2count:
                gene2count[gene] = [count, 1.0/len(genes)]
            else:
                currcount = gene2count[gene]
                gene2count[gene] = [currcount[0] + count, currcount[1] + 1.0/len(genes)]
    return gene2count

def combinationUsage( seqs, types ):
    comb2count = {} #key = combination of geneNames, val = [totalReads, uniqueSeqs]
    for seq in seqs.values():
        type2genes = {}
        totalCombinations = 1
        
        for type in types:
            genes = getGenes(seq, type)
            type2genes[type] = genes
            totalCombinations *= len(genes)
        
        if totalCombinations == 0:
            continue
        count = float(seq.count)/totalCombinations
       
        combs = type2genes[ types[0] ]
        for i in xrange(1, len(types)):
            type = types[i]
            currcombs = []
            for gene in type2genes[type]:
                for comb in combs:
                    currcombs.append( "|".join([comb, gene]) )
            combs = currcombs

        for comb in combs:
            if comb not in comb2count:
                comb2count[comb] = [count, 1.0/totalCombinations]
            else:
                currcount = comb2count[comb]
                comb2count[comb] = [ currcount[0] + count, currcount[1] + 1.0/totalCombinations ] 
    return comb2count

def getGene2count(seqs):
    #Single:
    type2gene2count = { 'v':{}, 'j':{}, 'd': {}, 'dj':{}, 'vj':{}, 'vdj':{} } 
    singletypes = ['v', 'j', 'd']
    for type in singletypes:
        gene2count = singleUsage(seqs, type)
        type2gene2count[type] = gene2count

    #Combination:
    combs = ['dj', 'vj', 'vdj']
    for comb in combs:
        types = [c for c in comb]
        comb2count = combinationUsage(seqs, types)
        type2gene2count[comb] = comb2count
    
    similarGenes = ['TRBV6-5', 'TRBV6-6']
    combineVgenes(type2gene2count, similarGenes) 

    return type2gene2count

def combineVgenes(type2gene2count, genes):
    '''Combine the genes in 'genes' as one gene
    '''
    newcounts = [0.0, 0.0]
    #Calculate combined counts
    for v, counts in type2gene2count['v'].iteritems():
        if v in genes:
            newcounts[0] += counts[0]
            newcounts[1] += counts[1]
    #Delete single genes
    for g in genes:
        if g in type2gene2count['v']:
            del type2gene2count['v'][g]
    #Add combined newgene
    newgene = '/'.join(genes)
    type2gene2count['v'][newgene] = newcounts
    
    #Combinations: vj, vdj
    combs = ['vj', 'vdj']
    for c in combs:
        if c not in type2gene2count:
            continue
        g2counts = {} #key = j or dj gene(s), val = counts
        delkeys = []
        gene2count = type2gene2count[c]
        #Calculate combined counts
        for g, counts in gene2count.iteritems(): #Each VJ or VDJ combination
            items = g.split('|')
            v = items[0] #current V
            if v in genes:
                delkeys.append(g)
                othergene= '|'.join(items[1:]) #current J or DJ
                if othergene not in g2counts:
                    g2counts[othergene] = [counts[0], counts[1]]
                else:
                    g2counts[othergene][0] += counts[0]
                    g2counts[othergene][1] += counts[1]
        #Delete combinations with single gene in genes
        for k in delkeys:
            del gene2count[k]
        #Add new combinations with new combined gene:
        for othergene, newcounts in g2counts.iteritems():
            newcomb = '|'.join([newgene, othergene])
            gene2count[newcomb] = newcounts
        #print gene2count
    

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

def addSamplingStats(type2gene2count, aggType2gene2count, i):
    #i is the order of the current sampling (base 0), or, it's the number of samplings that have already added to aggStats
    for type, gene2count in type2gene2count.iteritems():
        if type not in aggType2gene2count:
            aggType2gene2count[type] = {}
            for gene, counts in gene2count.iteritems():
                aggType2gene2count[type][gene] = [ [c] for c in counts ]
        else:
            aggGene2count = aggType2gene2count[type]
            for gene, counts in gene2count.iteritems():
                if gene not in aggGene2count:
                    aggGene2count[gene] = [ [0.0]*i + [c] for c in counts] #previous simulation didn't have this gene
                else:
                    aggCounts = aggGene2count[gene]
                    aggCounts[0].append(counts[0])
                    aggCounts[1].append(counts[1])
                    aggType2gene2count[type][gene] = aggCounts 

def avrSamplingStats(aggType2gene2count):
    #Average stats of the samplings:
    avrtype2gene2count = {}
    stdtype2gene2count = {}
    for type, gene2count in aggType2gene2count.iteritems():
        avrtype2gene2count[type] = {}
        stdtype2gene2count[type] = {}
        for gene, counts in gene2count.iteritems():
            meanReads = np.mean(counts[0])
            meanUniqs = np.mean(counts[1])
            avrtype2gene2count[type][gene] = [meanReads, meanUniqs]

            stdReads = np.std(counts[0])
            stdUniqs = np.std(counts[1])
            stdtype2gene2count[type][gene] = [stdReads, stdUniqs]
    return avrtype2gene2count, stdtype2gene2count

def usageTab(types, sample, avrstats, stdstats, type2genelist, outdir):
    for type in types:
        avrgene2count = {}
        stdgene2count = {}
        totalreads = 0
        totaluniqs = 0

        if type in avrstats:
            avrgene2count = avrstats[type]
            stdgene2count = stdstats[type]
            totalreads = sum([counts[0] for counts in avrgene2count.values()])
            totaluniqs = sum([counts[1] for counts in avrgene2count.values()])
            #if totalreads == 0 or totaluniqs == 0:
            #    raise ValueError("sample with zero read/sequence")

        if type in type2genelist:
            genes = type2genelist[type]
        else:
            genes = sorted( avrgene2count.keys() )
        typedir = os.path.join(outdir, type) 
        outfile = os.path.join(typedir, "%s-%s.txt" %(sample, type) )
        f = open(outfile, 'w') 
        f.write("Gene\tReads\t%Reads\tUniq\t%uniq\tStdReads\tStdUniq\n")
        #numpass = 0
        for g in genes:
            if g not in avrgene2count:
                sys.stderr.write("Gene %s is not in avrgene2count %s\n"  %(g, ','.join(avrgene2count.keys()) ))
                 
                avrcounts = [0.0, 0.0]
                stdcounts = [0.0, 0.0]
            else:
                #numpass += 1
                avrcounts = avrgene2count[g]
                stdcounts = stdgene2count[g]
            read = avrcounts[0]
            uniq = avrcounts[1]
            readPc = 0.0
            uniqPc = 0.0
            if totalreads > 0:
                readPc = iseqlib.getPc(read, totalreads)
            if totaluniqs > 0:
                uniqPc = iseqlib.getPc(uniq, totaluniqs)

            f.write("%s\t%d\t%f\t%d\t%f\t%d\t%d\n" %(g, read, readPc, uniq, uniqPc, stdcounts[0], stdcounts[1]))
        f.close()
        #if numpass == 0:
        #    raise ValueError("ERROR\n")

def geneUsedSample(avrstats):
    types = ['v', 'j', 'd', 'dj', 'vj', 'vdj']
    type2count = {} #key = genetype, val = count
    for type in types:
        if type not in avrstats:
            continue
        avrgene2count = avrstats[type]
        used = 0
        for counts in avrgene2count.values():
            if counts[0] > 0:
                used += 1
        type2count[type] = used
    return type2count

def geneUsed(avrstats, type2genelist, outfile):
    types = ['v', 'j', 'd', 'dj', 'vj', 'vdj']
    f = open(outfile, 'w')
    f.write("Genetype\tTotal\tUsed\tPercentage\n")
    type2count = geneUsedSample(avrstats)
    for type in types:
        if type not in type2count:
            continue
        used = type2count[type]
        if type in type2genelist:
            total = len(type2genelist[type])
            if total == 0:
                raise ValueError("Genetype %s has zero genes\n" %type)
            pc = 100.0*used/total
            f.write("%s\t%d\t%d\t%f\n" %(type, total, used, pc))
        else:
            f.write("%s\tNA\t%d\tNA\n" %(type, used))
    f.close()

def geneUsedSummary(sample2stats, type2genelist, group2samples, outfile, abs):
    #Row = sample, column = genetype
    f = open(outfile, 'w')
    types = ['d', 'j', 'v', 'dj', 'vj', 'vdj']
    f.write("Sample\t%s\n" %('\t'.join(types)))
    for group in sorted(group2samples.keys()):
        samples = group2samples[group]
        type2avrcount = {'d':0, 'j':0, 'v':0, 'dj':0, 'vj':0, 'vdj':0}
        #Print stats of each sample
        for sample in samples:
            f.write("%s" %sample)
            type2count = geneUsedSample( sample2stats[sample] )
            for type in types:
                count = 0
                if type in type2count:
                    count = type2count[type]
                
                type2avrcount[type] += count
                if abs:
                    f.write("\t%d" %count)
                else:#calculate percentage
                    if type in type2genelist:
                        total = len( type2genelist[type] )
                        if total == 0:
                            raise ValueError("Genetype %s has zero genes\n" %type)
                        pc = 100.0*count/total
                        f.write("\t%f" %pc)
                    else:
                        f.write("\tNA")
            f.write("\n")
        
        #Group average stats:
        f.write("%s" %group)
        for type in types:
            avrcount = float(type2avrcount[type])/len(samples)
            if abs:
                f.write("\t%d" % avrcount)
            else:
                if type in type2genelist:
                    total = len( type2genelist[type] )
                    if total == 0:
                        raise ValueError("Genetype %s has zero genes\n" %type)
                    pc = 100.0*avrcount/total
                    f.write("\t%f" %pc)
                else:
                    f.write("\tNA")
        f.write("\n")

    f.close()

#def getUsage(samples, outdir, type):
#    genes = getUnionGeneList(samples, type)
#    sys.stderr.write("Done getting uniqGeneList\n")
#    #Print out usage table for each sample:
#    for s in samples:
#        g2c = s.usage[type]
#        tabfile = os.path.join( outdir, "%s-%s.txt" %(s.name, type) )
#        f = open( tabfile, 'w') 
#        f.write("Gene\tTotal\tUniq\n")
#        for g in genes:
#            f.write( "%s\t%d\t%d\n" %(g, g2c[g][0], g2c[g][1]) )
#        f.close()

#def getVJusage(sample, type2gene2count, type2genelist, outdir, abs, uniq, std):
def getVJusage(sample, rowtype, coltype, type2gene2count, type2genelist, outdir, abs, uniq, std):
    #If abs is True: print absolute count, otherwise print frequencies.
    #If uniq is True: using the Uniq sequence Count as the unit, otherwise, use read count
    #rowtype = genetype represented by the rows, coltype = genetype represented by the columns 
    #(For exmple to represent vj recombinations, rows can be Vs and columns can be Js)
    if rowtype not in type2gene2count or coltype not in type2gene2count or (rowtype + coltype) not in type2gene2count:
        return

    v2c = type2gene2count[rowtype]
    j2c = type2gene2count[coltype]
    vj2c = type2gene2count[rowtype + coltype]
    totaluniqs = sum([c[1] for c in vj2c.values() ])
    totalreads = sum([c[0] for c in vj2c.values()])
    if totaluniqs == 0 or totalreads == 0:
        return
        #print vj2c
        #raise ValueError("Sample %s has zero sequence. rowtype: %s, coltype: %s. Totaluniqs: %d, totalreads: %d" %(sample, rowtype, coltype, totaluniqs, totalreads))

    if abs:
        outdirname = 'abs'
    else:
        outdirname = 'rel'
    if uniq:
        outdirname += "uniq"
    outdir = os.path.join(outdir, outdirname)
    if not std:
        file = os.path.join(outdir, "%s-vj.txt" %sample)
    else:
        file = os.path.join(outdir, "std%s-vj.txt" %sample)

    f = open(file, 'w')
    jgenes = [j for j in sorted(j2c.keys())] #column genes
    if coltype in type2genelist:
        jgenes = type2genelist[coltype]
    
    vgenes = [v for v in sorted(v2c.keys())] #row genes
    if rowtype in type2genelist:
        vgenes = type2genelist[rowtype]

    f.write( "\t%s\n" %( '\t'.join(jgenes) ) )
    
    for v in vgenes:
        if v == '' or re.search('undefined', v):
            continue
        f.write( "%s" %v )
        for j in jgenes:
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
                        count = float(count)/totaluniqs
                    else:
                        count = float(count)/totalreads
                    f.write("\t%f" %count)
        f.write("\n")
    f.close()

def getVJusageSample(sample, rowtype, coltype, avrstats, stdstats, type2genelist, outdir):
    #Print vj
    abs = True
    uniq = True
    std = True #if true, print the standard deviation
    getVJusage(sample, rowtype, coltype, avrstats, type2genelist, outdir, abs, uniq, not std)
    getVJusage(sample, rowtype, coltype, avrstats, type2genelist, outdir, abs, not uniq, not std)
    getVJusage(sample, rowtype, coltype, avrstats, type2genelist, outdir, not abs, uniq, not std)
    getVJusage(sample, rowtype, coltype, avrstats, type2genelist, outdir, not abs, not uniq, not std)
    #print stds:
    if stdstats:
        getVJusage(sample, rowtype, coltype, stdstats, type2genelist, outdir, abs, uniq, std)
        getVJusage(sample, rowtype, coltype, stdstats, type2genelist, outdir, abs, not uniq, std)
        getVJusage(sample, rowtype, coltype, stdstats, type2genelist, outdir, not abs, uniq, std)
        getVJusage(sample, rowtype, coltype, stdstats, type2genelist, outdir, not abs, not uniq, std)

#def getGeneUsage(sample, outdir):
#    '''Get V, D, J, VDJ, VJ and DJ usage
#    '''
#    getVJusage()
#        sample.setCounts()
#        sys.stderr.write("Done getting usage for %s\n" %sample.name)
#    
#    #Adding the average of all samples the the sample list
#    addAvrSample( samples )
#    sys.stderr.write("Done adding average and std sample\n")
#
#    vjUsage(samples, options.outdir)
#    sys.stderr.write("Done v usage and j usage\n")
#
#    vjoutdir = os.path.join( options.outdir, "vj")
#    system("mkdir -p %s" %vjoutdir)
#    #Generate VJ using the uniq sequence count or using the read count, relative or absolute count
#    abs = True
#    uniq = True
#    getVJusage(samples, vjoutdir, abs, not uniq)
#    sys.stderr.write("Done vj usage with absolute read count\n")
#    getVJusage(samples, vjoutdir, not abs, not uniq)
#    sys.stderr.write("Done vj usage with relative read count\n")
#    getVJusage(samples, vjoutdir, abs, uniq)
#    sys.stderr.write("Done vj usage with absolute uniqSeq count\n")
#    getVJusage(samples, vjoutdir, not abs, uniq)
#    sys.stderr.write("Done vj usage with relative read count\n")

