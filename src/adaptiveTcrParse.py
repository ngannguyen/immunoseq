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
        self.dgene = items[12]
        self.jgene = items[13]
        self.jties = items[14]
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
        
        header = ">%s;%s|%s|%s;size=%d" %(samplename, str(i), vstr, jstr, seq.normCount)
        if not useNorm:
            header = ">%s;%s|%s|%s;size=%d" %(samplename, str(i), vstr, jstr, seq.count)
        f.write("%s\n" %header)
        if options.nuc:
            f.write("%s\n" %seq.nuc)
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
    return 1

def main():
    parser = iseqlib.initOptions()
    parser.add_option('-i', '--indir', dest='indir')
    parser.add_option('-o', '--outdir', dest='outdir')
    parser.add_option('-r', '--raw', dest='raw', default=False, action='store_true', help='If specified, use raw counts instead of normalized counts. Default=%default')
    parser.add_option('-n', '--nucleotide', dest='nuc', default=False, action='store_true', help='If specified, print nucleotide sequences instead of amino acid sequences. Default=amino acid')
    parser.add_option('-p', '--productive', dest='productive', default=False, action="store_true", help="If specified, only prints out productive sequences. Default is to print all sequences")
    parser.add_option('-b', '--nonProductive', dest='nonProductive', default=False, action="store_true", help="If specified, only prints out non-productive sequences. Default is to print all.")
    parser.add_option('-v', '--uniqv', dest='uniqv', default=False, action="store_true", help="If specified, filter out sequences that mapped to multiple V genes")
    parser.add_option('-j', '--uniqj', dest='uniqj', default=False, action="store_true", help="If specified, filter out sequences that mapped to multiple J genes")
    parser.add_option('--jgenes', dest='jgenes', help='Required argument if would like to extend the 3\' end of sequences to get the full CDR3')
    parser.add_option('--jgeneToPos', dest='jgene2pos', help='Required argument if would like to extend the 3\' end of sequences to get the full CDR3')

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
        outfile = os.path.join(outdir,"%s.fa" % sample)
        printSeqs(seqs, outfile, useNorm, options)
        #printSeqs(seqs, outfile, False)


if __name__ == '__main__':
    main()














        

