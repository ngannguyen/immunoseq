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
        
        self.vs = self.vgene
        if self.vties != '':
            self.vs = ','.join( self.vties.split(', ') )

        self.dgene = items[12]
        self.jgene = items[13]
        self.jties = items[14]
        
        self.js = self.jgene
        if self.jties != '':
            self.js = ','.join( self.jties.split(', ') )

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
        
        self.cdr3nuc = self.nuc[self.vindex: self.jindex] 

def getNt2seqs(vs, js, nts, nt2genes2seqs):
    nt2seqs = {}
    for nt in nts:
        if nt in nt2seqs:
            continue

        if nt not in nt2genes2seqs:
            #raise ValueError('Unmatch nucleotide sequence %s\n' %nt)
            continue

        genes2seqs = nt2genes2seqs[nt]
        found = False
        for genes, seqs in genes2seqs.iteritems():
            if checkGenes(vs, js, genes[0], genes[1]): 
                nt2seqs[nt] = seqs
                found = True
                break
        if not found:
            #print nt, vs, js
            #print nt2genes2seqs[nt].keys()
            #raise ValueError('Unmatch nucleotide sequence %s\n' %nt)
            continue
    #for nt, seqs in nt2seqs.iteritems():
    #    print nt
    #    print [seq.cdr3nuc for seq in seqs]

    return nt2seqs

def isJoint(i, nts, nt2seqs):
    cat2count = {'v': 0, 'n2':0, 'd': 0, 'n1': 0, 'j':0} 
    joints = 0
    nojoints = 0
    print "i = %d" %i
    #print [nt[i] for nt in nts]
    #print nts
    #print [ (seq.n2index - seq.vindex, seq.dindex - seq.vindex, seq.n1index - seq.vindex, seq.jindex - seq.vindex)  for seq in nt2seqs.values() ]

    for nt in nts:
        if nt not in nt2seqs:
            continue

        seqs = nt2seqs[nt]
        isjoint = False
        for seq in seqs:
            pos = i + seq.vindex
            print "<%d-%d> <%d-%d>" %(seq.n2index, seq.dindex, seq.n1index, seq.jindex)
            print pos

            if (seq.n2index <= pos and pos < seq.dindex) or (seq.n1index <= pos and pos < seq.jindex):
                isjoint = True
                #break
            if pos < seq.n2index:
                cat2count['v'] += 1
            elif pos < seq.dindex:
                cat2count['n2'] += 1
            elif pos < seq.n1index:
                cat2count['d'] += 1
            elif pos < seq.jindex:
                cat2count['n1'] += 1
            else:
                cat2count['j'] += 1

        if isjoint:
            joints += 1
        else:
            nojoints += 1
    
    if joints > 0 and joints != len(nts):
        print "Inconsistent sequences\n"
        #print nts
        #print nt2seqs
        print joints
    
    #if joints > nojoints:
    #    return True
    #else:
    #    return False
    return joints, nojoints, cat2count

def getSeqStats(nts, nt2seqs):
    joints = 0
    nojoints = 0
    cat2count = {'v': 0, 'n2':0, 'd': 0, 'n1': 0, 'j':0} 
    
    l = len(nts[0])
    for nt in nts:
        if len(nt) != l:
            raise ValueError("Nucleotide sequences do not have the same length\n")

    for i in xrange(l): #each position, check to see if there is any mismatch
        letters = [nt[i] for nt in nts]
        mismatch = False
        for le in letters:
            if le != letters[0]:
                mismatch = True
                break
        if mismatch:
            jt, njt, tempcat2count = isJoint(i, nts, nt2seqs)
            for k, c in tempcat2count.iteritems():
                cat2count[k] += c

            if jt > njt:
                joints += 1
            elif jt < njt:
                nojoints += 1
    return joints, nojoints, cat2count

def getStats( header2nts, nt2genes2seqs, outfile ):
    #outfile = os.path.join(outdir, "aa2nts_stats.txt")
    
    f = open(outfile, 'w')
    f.write("Seq\tRatio\tJoints\tNot Joints\tTotal\tTotalLength\tv\tn2\td\tn1\tj\n")
    for header, nts in header2nts.iteritems():
        items = header.split('|')
        aa = items[0]
        vs = items[1]
        js = items[2]
        nt2seqs = getNt2seqs(vs, js, nts, nt2genes2seqs)
        joints, nojoints, cat2count = getSeqStats( nts, nt2seqs )
        if joints == 0 and nojoints == 0:
            continue
        ratio = 100.0
        if nojoints > 0:
            ratio = float(joints)/nojoints
        f.write("%s\t%f\t%d\t%d\t%d\t%d\t" %(header, ratio, joints, nojoints, joints + nojoints, len(aa)*3) )
        f.write("%d\t%d\t%d\t%d\t%d\n" %(cat2count['v'], cat2count['n2'], cat2count['d'], cat2count['n1'], cat2count['j']) )
    f.close()

def checkGenes(v1str, j1str, vstr, jstr):
    #Return true if v1s has at least 1 V in vs, and j1s has at least 1 J in js
    v1s = v1str.split(',')
    j1s = j1str.split(',')
    vs = vstr.split(',')
    js = jstr.split(',')
    hasJ = False
    hasV = False
    for v in v1s:
        if v in vs:
            hasV = True
            break
    for j in j1s:
        if j in js:
            hasJ = True
            break
    return hasJ and hasV

def unionGenes(genes1str, genes2str):
    genes1 = genes1str.split(',')
    genes2 = genes2str.split(',')
    genes = genes2
    for g in genes1:
        if g not in genes2:
            genes.append(g)
    return ','.join(genes)

def addSeq(seqs, seq):
    if seq.cdr3nuc not in seqs:
        seqs[ seq.cdr3nuc ] = { (seq.vs, seq.js): [seq] }
    else:
        genes2seqs = seqs[ seq.cdr3nuc ]
        nomatch = True
        for genes, seqlist in genes2seqs.iteritems():
            if checkGenes( seq.vs, seq.js, genes[0], genes[1] ):
                newvs = unionGenes( seq.vs, genes[0] )
                newjs = unionGenes( seq.js, genes[1] )
                del genes2seqs[genes]
                seqlist.append(seq)
                genes2seqs[ (newvs, newjs) ] = seqlist
                nomatch = False
                break
        if nomatch:
            genes2seqs[ (seq.vs, seq.js) ] = [seq]
    return
                
def readTsvFile(file, seqs, options):
    f = sys.stdin
    #hasNorm = True
    if file != '-' and os.path.exists(file):
        f = open(file, 'r')
    #seqs = {} #key = nucleotide sequence, val = Seq()
    for line in f:
        items = line.strip().split('\t')
        if len(items) < 27 or items[0] =='sequenceID':
            continue
        seq = Sequence(line)
        #if seq.normFreq < 0 or seq.normCount < 0:
        #    hasNorm = False
        
        if seq.status != 'Productive':
            continue
           
        fillInSeq(seq, options.jgenes, options.jgene2pos)
        addSeq( seqs, seq )
    
    f.close()
    #return seqs, hasNorm

def readTsvFiles(indir, options):
    seqs = {} #key = nucleotide sequence, val = Seq()
    for file in os.listdir(indir):
        if os.path.isdir( os.path.join(indir, file) ):
            continue
        readTsvFile( os.path.join(indir, file), seqs, options )
    return seqs

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
    return 1

def readInfile( infile ):
    header2nts = {} #key = aa, vals = [nts]
    f = open(infile, 'r')
    header = ''
    nts = []
    for line in f:
        line = line.strip()
        if line[0] == '>':
            if header != '':
                header2nts[header] = nts
            header = line.lstrip('>')
            nts = []
        else:
            nts.append(line)
    if header != '':
        header2nts[header] = nts
    f.close()
    return header2nts

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

def main():
    parser = iseqlib.initOptions()
    parser.add_option('-i', '--infile', dest='infile', help='amino acid to nucleotide sequences file. Format: ')
    parser.add_option('-t', '--tsvdir', dest='tsvdir', help='Directory of tsv files')
    parser.add_option('-o', '--outfile', dest='outfile', help='Output file')
    parser.add_option('--jgenes', dest='jgenes', default='/hive/users/nknguyen/immuno/imgt/TRB/jgenes.txt', help='Required argument if would like to extend the 3\' end of sequences to get the full CDR3')
    parser.add_option('--jgeneToPos', dest='jgene2pos', default='/hive/users/nknguyen/immuno/imgt/TRB/jmappedF.txt', help='Required argument if would like to extend the 3\' end of sequences to get the full CDR3')
    
    options, args = parser.parse_args()
    if options.jgenes:
        options.jgenes = readGeneFile( options.jgenes )
    if options.jgene2pos:
        options.jgene2pos = readGene2pos( options.jgene2pos )

    header2nts = readInfile( options.infile )
    nt2genes2seqs = readTsvFiles( options.tsvdir, options )
    getStats( header2nts, nt2genes2seqs, options.outfile )

if __name__ == '__main__':
    main()














        

