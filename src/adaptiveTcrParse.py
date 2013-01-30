#!/usr/bin/env python

#nknguyen at soe ucsc edu
#Jan 10 2011
#Parse adaptive TCR's output tsv file
#
#Input: Directory that contains Adaptive Biotechnologies tsv files:
#sequenceID container nucleotide aminoAcid normalizedFrequency normalizedCopy rawFrequency copy cdr3Length VFamilyName VGeneName VTies DGeneName JGeneName JTies VDeletion d5Deletion d3Deletion JDeletion n2Insertion n1Insertion sequenceStatus VIndex n1Index n2Index DIndex JIndex
#
#Output: Output fasta files of productive sequences to output directory
#   The header of each sequence has the following format: 
#   >sampleName;userDefinedCloneID|Vgenes|Jgenes;size=#

import os, re, sys, copy
import immunoseq.lib.immunoseqLib as iseqlib
from scipy.stats.stats import pearsonr, spearmanr, kendalltau

class Sequence:
    def __init__(self, line):
        items = line.strip().split('\t')
        if len(items) < 27:
            sys.stderr.write('Wrong tsv format. Expected 27 fields, only have %d\n%s\n' %(len(items), line))
            sys.exit(1)
        self.id = items[0]
        self.nuc = items[2]
        self.aa = items[3]
        self.normFreq = -1.0
        self.normCount = -1
        if items[4] != '':
            self.normFreq = float(items[4])
        if items[5] != '':
            self.normCount = int(items[5])
        self.freq = float(items[6])
        self.count = int(items[7])
        self.cdr3len = int(items[8])
        self.vfam = items[9]
        self.vgene = items[10]
        self.vties = items[11]

        vgenes = self.vgene.split('/')
        if len(vgenes) > 1:
            self.vgene = vgenes[0]
            vs = self.vties.split(', ')
            for v in vgenes[1:]:
                if v not in vs:
                    vs.append(v)
            self.vties = ', '.join(vs)

        self.vgene = self.vgene.split('*')[0]
        
        self.dgene = items[12]
        self.jgene = items[13]
        self.jties = items[14]
        
        jgenes = self.jgene.split('/')
        if len(jgenes) > 1:
            self.jgene = jgenes[0]
            js = self.jties.split(', ')
            for j in jgenes[1:]:
                if j not in js:
                    js.append(j)
            self.jties = ', '.join(js)
        self.jgene = self.jgene.split('*')[0]
        
        self.vdel = int(items[15])
        self.d5del = int(items[16])
        self.d3del = int(items[17])
        self.jdel = int(items[18])
        self.n2ins = int(items[19])
        self.n1ins = int(items[20])
        self.status = items[21]
        self.vindex = int(items[22])
        self.n1index = int(items[23])
        self.n2index = int(items[24])
        self.dindex = int(items[25])
        self.jindex = int(items[26])
        
        self.cdr3nuc = self.nuc[self.vindex:self.jindex]
        self.inframenuc = self.nuc[ self.vindex%3: len(self.nuc) - ((len(self.nuc) - self.vindex%3)%3) ]
        self.longaa = iseqlib.nt2aa( self.inframenuc ) 

        self.vpos2snp = None
        self.jpos2snp = None

class Snp:
    def __init__(self, refname, pos, refAllele, altAllele, score):
        self.refname = refname
        self.pos = pos
        self.ref = refAllele
        self.alt = altAllele
        self.score = score
        self.alt2score = {}
        self.alt2score[altAllele] = score #key = alt, val = score of that alternate allele
        self.totalReads = -1

    def updateScore(self, newscore):
        self.score = newscore

    def addAltAllele(self, alt, score):
        alts = self.alt.split(',')
        if alt not in alts:
            self.alt += ",%s" %alt
            self.alt2score[alt] = score
        else:
            self.alt2score[alt] += score
    
    def setTotalReads(self, totalReads):
        self.totalReads = totalReads

    def filterNoiseAltAlleles(self):
        for a in self.alt2score.keys():
            if self.alt2score[a] == 0:
                del self.alt2score[a]
        self.alt = ','.join( self.alt2score.keys() )

    def getStr(self):
        return "%s\t%d\t%s\t%s\t%d" %(self.refname, self.pos, self.ref, self.alt, self.score)

def readFile(file):
    f = sys.stdin
    hasNorm = True
    if file != '-' and os.path.exists(file):
        f = open(file, 'r')
    seqs = []
    for line in f:
        items = line.strip().split('\t')
        if len(items) < 27 or items[0] =='sequenceID':
            continue
        seq = Sequence(line)
        if seq.normFreq < 0 or seq.normCount < 0:
            hasNorm = False
        seqs.append(seq)
        #if seq.status == 'Productive':
        #    seqs.append(seq)
    f.close()
    return seqs, hasNorm

def printSeqs(seqs, outfile, useNorm, options):
    samplename = os.path.basename(outfile).rstrip(".fa")
    f = sys.stdout
    if outfile != '-':
        f = open(outfile, 'w')
    
    i = 0
    numFailedExtends = 0
    for seq in seqs:
        #Filtering
        if options.nonProductive and seq.status == 'Productive':
            continue
        if options.productive and seq.status != 'Productive':
            continue
        if options.uniqv and seq.vties != '':
            continue
        if options.uniqj and seq.jties != '':
            continue
        if seq.vgene == '(undefined)' or seq.jgene == '(undefined)':
            continue

        if options.jgenes:
            extendok = fillInSeq(seq, options.jgenes, options.jgene2pos)
            if extendok == -1:
                numFailedExtends += 1

        vstr = seq.vgene
        if seq.vties != '':
            vstr = ','.join( seq.vties.split(', ') )
        jstr = seq.jgene
        if seq.jties != '':
            jstr = ','.join( seq.jties.split(', ') )
        
        header = ">%s;%s|%s|%s|%s;size=%d" %(samplename, str(i), vstr, jstr, seq.dgene, seq.normCount)
        if not useNorm:
            header = ">%s;%s|%s|%s|%s;size=%d" %(samplename, str(i), vstr, jstr, seq.dgene, seq.count)
        f.write("%s\n" %header)
        if options.nuc:
            #f.write("%s\n" %seq.nuc)
            if options.long:
                f.write("%s\n" %seq.inframenuc)
            else:
                f.write("%s\n" %seq.cdr3nuc)
        else:
            if options.long:
                f.write("%s\n" %seq.longaa)
            else:
                f.write("%s\n" %seq.aa)
        i+=1

    if options.jgenes:
        sys.stderr.write('Number of failed extensions is %d\n' %numFailedExtends)

#====== THIS SECTION IS TO INFER THE RIGHT END OF THE CDR3 SEQUENCES ========
def readGeneFile(file):
    genes = {} #key = name, val = seq
    f = open(file, 'r')
    header = ''
    seq = ''
    for line in f:
        if line[0] == '>':
            if seq != '':
                genes[header] = seq
            header = line.strip().lstrip('>')
            seq = ''
        else:
            seq = line.strip()
    if seq != '' and header != '':
        genes[header] = seq
    f.close()
    return genes

def readGene2pos(file):
    gene2pos = {}
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        items = line.split()
        gene2pos[items[0]] = int(items[1])
    f.close()
    return gene2pos

def fillInSeq(seq, genes, gene2pos):
    '''The earlier version of adaptiveTCR tsv files did not always have the complete CDR3 sequences (because reads only support partial CDR3)
    This function looked up the original J gene sequences that the clone to mapped to, and fill in the its sequences to have complete CDR3.
    '''
    jfrag = seq.nuc[seq.jindex:].upper()
    jseq = genes[seq.jgene].upper()

    pos = gene2pos[seq.jgene] + 3
    #find where the jfrag is in the jseq:
    matches = re.finditer(jfrag, jseq)
    starts = [ m.start() for m in matches ]
    #Now finding all possible matches:
    fullseqs = []
    newcdr3s = []

    for s in starts:
        e = s + len(jfrag)
        if s <= 15 and e >= 10 and e < pos: #the primer should be somewhere to the left (5'end) of the J gene and should cover somewhere after the 10 leftmost nucleotides
            extraseq = jseq[e :pos]
            cdr3len = len(seq.nuc) - seq.vindex + len(extraseq)
            if cdr3len%3 == 0: #inframe
                #Translate to cdr3:
                newNuc = seq.nuc + extraseq
                newaa = iseqlib.nt2aa( newNuc[seq.vindex:] )
                if newaa[ len(newaa) -1 ] == 'F' and '*' not in newaa:
                    if seq.aa not in newaa:
                        sys.stderr.write('Warning: the new infered aa doesnot contain current aa\n')
                    fullseqs.append( newNuc )
                    newcdr3s.append( newaa )
    if len(fullseqs) == 0:
        sys.stderr.write('Attempted to fill in the right side of the nuc and CDR3 sequences. Zero productive matches found. Sequence looked at was: %s, %s, %s\n' %(seq.id, seq.nuc, seq.aa))
        sys.stderr.write("jfrag: %s, jseq: %s, matchesStarts: *%s*\n" %(jfrag, jseq, ','.join(starts)) )
        return -1
    elif len(fullseqs) > 1:
        sys.stderr.write('Attempted to fill in the right side of the nuc and CDR3 sequences. Multiple productive matches found - could not decide. Sequence looked at was: %s, %s, %s\n' %(seq.id, seq.nuc, seq.aa))
        return -1
    else:
        seq.nuc = fullseqs[0]
        seq.aa = newcdr3s[0]
        seq.jindex = seq.vindex + len(seq.aa)*3 
        seq.cdr3nuc = seq.nuc[seq.vindex:seq.jindex]
        seq.inframenuc = seq.nuc[seq.vindex%3: len(seq.nuc) - ((len(seq.nuc) - seq.vindex%3)%3)]
        seq.longaa = iseqlib.nt2aa( seq.inframenuc ) 
    return 1

#========== CALLING SNPs ===========
def findSnps(seqs, useNorm, options):
    #samplename = os.path.basename(outfile).rstrip(".txt")
    gene2pos2snp = {} 
    gene2count = {}

    for seq in seqs:
        #Filtering
        if options.nonProductive and seq.status == 'Productive':
            continue
        if options.productive and seq.status != 'Productive':
            continue
        if options.uniqv and seq.vties != '':
            continue
        if options.uniqj and seq.jties != '':
            continue

        if not options.jgenes or not options.jgene2pos or not options.vgenes or not options.vgene2pos:
            raise ValueError("Cannot find SNP without knowing the reference gene sequences. Please specify all jgenes, jgene2pos, vgene, vgene2pos.\n")

        count = seq.normCount
        if not useNorm:
            count = seq.count
        
        if seq.vgene not in gene2count:
            gene2count[seq.vgene] = count
        else:
            gene2count[seq.vgene] += count
        if seq.jgene not in gene2count:
            gene2count[seq.jgene] = count
        else:
            gene2count[seq.jgene] += count

        #V SNPs
        vpos2snp = getVsnps(seq, options.vgenes, options.vgene2pos) #vpos2snp{ pos : (ref, alt) } (position along the Vgene where the nucleotide is different from the reference allele)
        seq.vpos2snp = vpos2snp
        #print seq.vpos2snp

        #update gene2pos2snp ( {gene: {pos: Snp}} )
        for pos, (ref, alt) in vpos2snp.iteritems():
            snp = Snp(seq.vgene, pos, ref, alt, count)
            if seq.vgene not in gene2pos2snp:
                gene2pos2snp[seq.vgene] = { pos:snp }
            elif pos not in gene2pos2snp[seq.vgene]:
                gene2pos2snp[seq.vgene][pos] = snp
            else:
                snp = gene2pos2snp[seq.vgene][pos]
                if snp.ref != ref:
                    raise RuntimeError("Conflict in reference Allele.\n")
                snp.updateScore( snp.score + count )
                snp.addAltAllele( alt, count )
            
            #HACK - DEBUG
            #if seq.vgene == 'TRBV7-9' and pos == 262:
            #    print snp.ref, snp.alt, count, snp.score

        #J SNPs:
        jpos2snp = getJsnps(seq, options.jgenes, options.jgene2pos)
        seq.jpos2snp = jpos2snp

        #update gene2pos2snp ( {gene: {pos: Snp}} )
        for pos, (ref, alt) in jpos2snp.iteritems():
            snp = Snp(seq.jgene, pos, ref, alt, count)
            if seq.jgene not in gene2pos2snp:
                gene2pos2snp[seq.jgene] = { pos:snp }
            elif pos not in gene2pos2snp[seq.jgene]:
                gene2pos2snp[seq.jgene][pos] = snp
            else:
                snp = gene2pos2snp[seq.jgene][pos]
                if snp.ref != ref:
                    raise RuntimeError("Conflict in reference Allele.\n")
                snp.updateScore( snp.score + count )
                snp.addAltAllele( alt, count )

    #filter out alternative alleles with low count (and most likely to be noise/sequencing errors):
    avrNoiseFreq, gene2pos2noiseAlleles = filterNoiseAlleles(gene2pos2snp, gene2count, options.minAlleleScore, options.minAlleleFreq)
    #sampleName = os.path.basename(outfile).split('.')[0]
    #sys.stderr.write("Sample %s, Average Noise Frequency is %f%%\n" %(sampleName, avrNoiseFreq) )
   
    #DEBUG:
    #for seq in seqs:
    #    if seq.vpos2snp and len(seq.vpos2snp) > 0:
    #        print seq.vpos2snp
    #END DEBUG

    return gene2pos2snp, gene2count, gene2pos2noiseAlleles

def printSnps(outfile, gene2pos2snp, options):
    f = sys.stdout
    if outfile != '-':
        f = open(outfile, 'w')
    
    #print to output file:
    if not options.mapToRef:
        f.write("#Genename\tPosition\tRef\tAlt\tReadCount\tPercentage\n")
        for gene in sorted(gene2pos2snp.keys()):
            pos2snp = gene2pos2snp[gene]
            for pos in sorted(pos2snp.keys()):
                snp = pos2snp[pos]
                f.write( "%s\t%.3f\n" % (snp.getStr(), snp.score*100.0/snp.totalReads) )
   
    else: #PRINT in the PersonalVariant format and in the reference chromosome coordinates to display on the browser
        f.write("#chrom\tchromStart\tchromEnd\talleles\talleleCount\talleleFreq\talleleScores\n")
        for gene in sorted(gene2pos2snp.keys()):
            if gene not in options.gene2offset:
                sys.stderr.write("Gene %s not found in the reference.\n" %gene)
            (goffset, gstrand, gchr) = options.gene2offset[gene]

            pos2snp = gene2pos2snp[gene]
            for pos in sorted(pos2snp.keys()):
                snp = pos2snp[pos]
                refAlleleCount = snp.totalReads - snp.score
                chromStart = goffset + pos
                if gstrand == '-':
                    chromStart = goffset - pos -1
                chromEnd = chromStart + 1
                alleles = snp.alt.strip(',').split(',')
                alleleStr = snp.ref 
                alleleCount = 1
                freqstr = "%d" %refAlleleCount
                for a in alleles:
                    if a != '':
                        alleleCount += 1
                        alleleStr += "/%s" %a
                        freqstr += ",%d" %snp.alt2score[a]

                f.write("%s\t%d\t%d\t%s\t%d\t%s\t%s\n" %(gchr, chromStart, chromEnd, alleleStr, alleleCount, freqstr, freqstr))

def filterNoiseAlleles(gene2pos2snp, gene2count, mincount, minfreq):
    avrNoiseReads = 0
    total = 0
    gene2pos2noiseAlleles = {} #gene:{pos : [list of noise allele]}
    for gene, pos2snp in gene2pos2snp.iteritems():
        totalReads = gene2count[gene]

        for pos in pos2snp.keys():
            snp = pos2snp[pos]
            #DEBUG:
            if snp.score != sum(snp.alt2score.values()):
                sys.stderr.write("Score %d != sum of alt scores: %d\n" %(snp.score, sum(snp.alt2score.values()) ) )
            #END DEBUG

            noiseReadCount = 0
            altNoiseReadCount = 0
            noisealleles = []
            noisecounts = []

            #check ref allele:
            refAlleleCount = totalReads - snp.score
            refAlleleFreq = refAlleleCount*100.0/totalReads
            if refAlleleCount < mincount or refAlleleFreq < minfreq: #noise
                noiseReadCount += refAlleleCount
                noisealleles.append( snp.ref )
                noisecounts.append( refAlleleCount )
            
            #check all the alternative alleles:
            alleles = snp.alt.split(',')
            for a in alleles:
                acount = snp.alt2score[a]
                afreq = acount*100.0/totalReads
                if acount < mincount or afreq < minfreq:
                    snp.alt2score[a] = 0 
                    noiseReadCount += acount
                    altNoiseReadCount += acount
                    if a not in noisealleles:
                        noisealleles.append(a)
                        noisecounts.append(acount)
            
            #Adjust totalReads and score:
            if snp.score == altNoiseReadCount: #all alternative alleles were noise, so Not a SNP, remove from the snp list
                del pos2snp[pos]
            else:
                snp.setTotalReads( totalReads - noiseReadCount )
                snp.updateScore( snp.score - altNoiseReadCount )
                snp.filterNoiseAltAlleles()

            #Update gene2pos2noiseAlleles:
            if len(noisealleles) > 0:
                if gene not in gene2pos2noiseAlleles:
                    gene2pos2noiseAlleles[gene] = { pos: noisealleles }
                else:
                    gene2pos2noiseAlleles[gene][pos] = noisealleles

            #DEBUG
            #print gene, pos, snp.ref, snp.alt, noisealleles, noisecounts, totalReads - noiseReadCount, altNoiseReadCount
            #END DEBUG

            #Noise frequency
            avrNoiseReads += noiseReadCount
            total += totalReads

    #Debugging/testing:
    for gene, pos2snp in gene2pos2snp.iteritems():
        for pos, snp in pos2snp.iteritems():
            if snp.alt == '':
                sys.stderr.write("SNP with no alternative alleles: %s\n" %snp.getStr())
            for alt, score in snp.alt2score.iteritems():
                if score < mincount or score*100.0/gene2count[gene] < minfreq:
                    sys.stderr.write("SNP with noie alleles: %s\n" %snp.getStr())
    #End Debugging/testing:

    return avrNoiseReads*100.0/total, gene2pos2noiseAlleles
    
def filterNoiseReads(seqs, gene2pos2noiseAlleles, noiseFile, options):
    newseqs = []
    f = open(noiseFile, 'w')
    for seq in seqs:
        count = seq.normCount
        if options.raw:
            count = seq.count
        
        if not seq.vpos2snp and not seq.jpos2snp:
            newseqs.append(seq)
            continue
        
        #print "has vpos2snp and jpos2snp"
        flagnoise = False #true if the seq contains sequencing noise, and consequently will be filtered out
        
        #CHECK V GENE
        if seq.vgene in gene2pos2noiseAlleles:
            pos2noiseAlleles = gene2pos2noiseAlleles[seq.vgene]
            #Check to see if the sequence contains 'alternative allele' which is likely to be noise
            for pos, (ref, alt) in seq.vpos2snp.iteritems():
                if pos not in pos2noiseAlleles:
                    continue
                noiseAlleles = pos2noiseAlleles[ pos ]
                if alt in noiseAlleles:
                    flagnoise = True
                    #print "alt %s, count %d, freq %f, gene %s, pos %d, noisealleles %s" %(alt, seq.count, seq.freq, seq.vgene, pos, ','.join(noiseAlleles))
                    break
            if not seq.vpos2snp: #No alternative allele, check to see if the sequence contains the reference allele, but which is likely to be noise
                for pos, noiseAlleles in pos2noiseAlleles.iteritems():
                    seqAllele = seq.nuc[ pos - seq.vrefstart ] 
                    noiseAlleles = pos2noiseAlleles[ pos ]
                    if seqAllele in noiseAlleles:
                        flagnoise = True
                        #print "REF %s, count %d, freq %f, gene %s, pos %d, noisealleles %s" %(seqAllele, seq.count, seq.freq, seq.vgene, pos, ','.join(noiseAlleles))
                        break
        
        #CHECK J GENE
        elif seq.jgene in gene2pos2noiseAlleles:
            pos2noiseAlleles = gene2pos2noiseAlleles[seq.jgene]
            #Check to see if the sequence contains 'alternative allele' which is likely to be noise
            for pos, (ref, alt) in self.jpos2snp.iteritems():
                if pos not in pos2noiseAlleles:
                    continue
                noiseAlleles = pos2noiseAlleles[ pos ]
                if alt in noiseAlleles:
                    flagnoise = True
                    break
            #Check to see if the sequence contains the reference allele, but which is likely to be noise
            if not seq.jpos2snp:
                for pos, noiseAlleles in pos2noiseAlleles.iteritems():
                    seqAllele = seq.nuc[ pos - seq.jrefstart ] 
                    if seqAllele in noiseAlleles:
                        flagnoise = True
                        break
         
        if flagnoise:
            f.write( "%d\t%f\n" %(count, seq.freq) )
        else:
            newseqs.append(seq)
    f.close()
    return newseqs

def readGene2ref(file):
    #TRBV1 141645479 + chr7
    gene2offset = {}

    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        items = line.split('\t')
        if len(items) != 4:
            raise ValueError("Wrong file format. File %s. Expected format:<genename>\\t<start>\\t<strand>\\t<refname>\n" %file)
        gene2offset[items[0]] = (int(items[1]), items[2], items[3])
    f.close()
    return gene2offset

def getSnps(seq, refseq, refstart):
    pos2snp = {} #key = position, val = [ref, alt]
    #check for SNPs:
    for i, refNuc in enumerate(refseq):
        nuc = seq[i]
        if nuc != refNuc:
            pos2snp[i + refstart] = (refNuc, nuc)
    return pos2snp
    
def getVsnps(seq, genes, gene2pos):
    '''Compare the V-gene portion to the 5' end of the CDR3 region of the clone sequence with the corresponding reference V-gene sequence to look for SNPs.
    '''
    vfrag = seq.nuc[:seq.vindex].upper()
    vseq = genes[seq.vgene].upper()
    
    endPos = gene2pos[seq.vgene] - 3 #exclusive
    startPos = endPos - len(vfrag) #inclusive, base 0
    if startPos < 0:
        sys.stderr.write("clone goes upstream beyond the start of V gene\n")
        #vfrag = vfrag[abs(startPos):]
        #startPos = 0

    refvfrag = vseq[startPos:endPos]
    seq.vrefstart = startPos
    return getSnps(vfrag, refvfrag, startPos)
    
def getJsnps(seq, genes, gene2pos):
    pos2snp = {}
    jfrag = seq.nuc[seq.jindex:].upper()
    jseq = genes[seq.jgene].upper()
    
    #Find where the jfrag is in the jseq:
    matches = re.finditer(jfrag, jseq)
    starts = [ m.start() for m in matches ]
    okstarts = []

    pos = gene2pos[seq.jgene] + 3
    for s in starts:
        e = s + len(jfrag)
        if s <= 15 and e >= 10 and e < pos: #the primer should be somewhere to the left (5'end) of the J gene and should cover somewhere after the 10 leftmost nucleotides
            extraseq = jseq[e :pos]
            cdr3len = len(seq.nuc) - seq.vindex + len(extraseq)
            if cdr3len%3 == 0: #inframe
                #Translate to cdr3:
                newNuc = seq.nuc + extraseq
                newaa = iseqlib.nt2aa( newNuc[seq.vindex:] )
                if newaa[ len(newaa) -1 ] == 'F' and '*' not in newaa:
                    if seq.aa not in newaa:
                        sys.stderr.write('Warning: the new infered aa doesnot contain current aa\n')
                    okstarts.append(s)
    if len(okstarts) > 1:
        sys.stderr.write("Looking for SNPs if any on the J fragment. However, mapped to multiple Js. Sequence: %s, %s, %s\n" %(seq.id, seq.nuc, seq.aa))
        sys.exit(1)

    startPos = starts[0]
    endPos = min( [ startPos + len(jfrag), gene2pos[seq.jgene] + 3 ] )
    if endPos > len(jseq):
        sys.stderr.write("Clone goes beyond downstream of J gene\n")
    refjfrag = jseq[startPos:endPos]
    seq.jrefstart = startPos
    return getSnps(jfrag, refjfrag, startPos)

#========== DONE CALLING SNPs ===========

def main():
    parser = iseqlib.initOptions()
    parser.add_option('-i', '--indir', dest='indir')
    parser.add_option('-o', '--outdir', dest='outdir')
    parser.add_option('-r', '--raw', dest='raw', default=False, action='store_true', help='If specified, use raw counts instead of normalized counts. Default=%default')
    parser.add_option('-l', '--long', dest='long', default=False, action='store_true', help='If specified, print inframe (nt/aa) sequences instead of just the CDR3 parts. Default=%default')
    parser.add_option('-n', '--nucleotide', dest='nuc', default=False, action='store_true', help='If specified, print nucleotide sequences instead of amino acid sequences. Default=amino acid')
    parser.add_option('-p', '--productive', dest='productive', default=False, action="store_true", help="If specified, only prints out productive sequences. Default is to print all sequences")
    parser.add_option('-b', '--nonProductive', dest='nonProductive', default=False, action="store_true", help="If specified, only prints out non-productive sequences. Default is to print all.")
    parser.add_option('-v', '--uniqv', dest='uniqv', default=False, action="store_true", help="If specified, filter out sequences that mapped to multiple V genes")
    parser.add_option('-j', '--uniqj', dest='uniqj', default=False, action="store_true", help="If specified, filter out sequences that mapped to multiple J genes")
    parser.add_option('-s', '--snp', dest='snp', default=False, action="store_true", help="If specified, get SNPs")
    parser.add_option('--filterBySnpNoise', dest='filterBySnpNoise', default=False, action="store_true", help="If specified, print out fasta files that exclude errorneous reads. Errorneous reads are those that contain alternative alleles different from the reference genes with low frequency, and therefore more likely to be sequencing noise.")
    parser.add_option('--noFasta', dest='noFasta', default=False, action="store_true", help="If specified, don't convert to fasta files")
    parser.add_option('--jgenes', dest='jgenes', default='/hive/users/nknguyen/immuno/imgt/TRB/jgenes.txt', help='Required argument if would like to extend the 3\' end of sequences to get the full CDR3')
    parser.add_option('--jgeneToPos', dest='jgene2pos', default='/hive/users/nknguyen/immuno/imgt/TRB/jmappedF.txt', help='Required argument if would like to extend the 3\' end of sequences to get the full CDR3')
    parser.add_option('--vgenes', dest='vgenes', default= '/hive/users/nknguyen/immuno/imgt/TRB/vgenes.txt', help='Required argument if would like to check for SNPs in the V region')
    parser.add_option('--vgeneToPos', dest='vgene2pos', default='/hive/users/nknguyen/immuno/imgt/TRB/vmapped312.txt', help='Required argument if would like to check for SNPs in the V region')
    parser.add_option('--mapToRef', dest='mapToRef', default='/hive/users/nknguyen/immuno/imgt/TRB/genesToHg18.txt', help='If specified, map the SNPs back to the reference')
    parser.add_option('--noMapToRef', dest='noMapToRef', action='store_true', default=False)
    parser.add_option('--minAlleleScore', dest='minAlleleScore', type='int', default=100, help='Minimum reads required for an allele to be called. Default=%default')
    parser.add_option('--minAlleleFreq', dest='minAlleleFreq', type='float', default=1.0, help='Minimum percentage required for an allele to be called. Default=%default')

    options, args = parser.parse_args()
    if (options.jgenes and not options.jgene2pos) or (options.jgene2pos and not options.jgenes):
        parser.error('jgenes and jgene2pos are required if would like to extend the 3 prime end to get the full CDR3 sequences.\n')
    if options.jgenes and options.nonProductive:
        parser.error('Cannot extend the CDR3 sequences of non-productive sequences.\n')
    if options.nonProductive and not options.nuc:
        parser.error('No amino acid options is available for non-productive sequences. Please choose the nucleotide options.\n')
    if options.jgenes:
        options.jgenes = readGeneFile( options.jgenes )
    if options.jgene2pos:
        options.jgene2pos = readGene2pos( options.jgene2pos )
    if options.vgenes:
        options.vgenes = readGeneFile( options.vgenes )
    if options.vgene2pos:
        options.vgene2pos = readGene2pos( options.vgene2pos )
    
    if options.noMapToRef:
        options.mapToRef = None
    if options.mapToRef:
        options.gene2offset = readGene2ref(options.mapToRef)

    indir = options.indir
    outdir = options.outdir
    useNorm = not options.raw

    for file in os.listdir(indir):
        sample = file.split('.')[0]
        if len(file.split('.')) < 2 or file.split('.')[1] != 'tsv':
            continue
        seqs, hasNorm = readFile( os.path.join(indir,file) )
        if not hasNorm and useNorm:
            useNorm = False
            sys.stderr.write("Some sample/sequence does not have normalized frequency or normalized count. Using raw count instead\n")
        if not options.noFasta:
            outfile = os.path.join(outdir,"%s.fa" % sample)
            printSeqs(seqs, outfile, useNorm, options)
        
        if options.snp:
            snpoutfile = os.path.join(outdir, "%s.snp.txt" %sample)
            gene2pos2snp, gene2count, gene2pos2noise = findSnps(seqs, useNorm, options)
            #print gene2pos2noise
            printSnps(snpoutfile, gene2pos2snp, options)
            
            #DEBUG:
            #for seq in seqs:
            #    if seq.vpos2snp and len(seq.vpos2snp) > 0:
            #        print seq.vpos2snp
            #END DEBUG


            if options.filterBySnpNoise:
                noisefile = os.path.join(outdir, "%s-noise.txt" %sample)
                newseqs = filterNoiseReads(seqs, gene2pos2noise, noisefile, options)
                outfile = os.path.join(outdir, "%s-filtered.fa" %sample)
                printSeqs(newseqs, outfile, useNorm, options)



if __name__ == '__main__':
    main()














        

