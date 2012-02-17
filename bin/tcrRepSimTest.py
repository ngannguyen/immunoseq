#!/usr/bin/env python

"""
nknguyen soe ucsc edu
Feb 09 2012
Unit Test for tcrRep
"""

import tcrRepSim
import unittest

#class TestCheckOptions( unittest.TestCase ):
#    def test_input_types(self):
#        '''checkOptions raises errors if input type is wrong'''

class TestGetOverlap( unittest.TestCase ):
    aaUsageFile = "../data/tcrRepSimTest/aaUsage.txt"
    codonUsageFile = "../data/tcrRepSimTest/codonUsage.txt"
    aaUsage = tcrRepSim.readAaUsage(aaUsageFile)
    codonUsage = tcrRepSim.readCodonUsage(codonUsageFile)
    seqs = tcrRepSim.getUniqSeqRep(10, aaUsage[9], codonUsage)
    tcrRepSim.getRep( seqs, 50, None )
    samseqs = tcrRepSim.sampling( seqs, 20 )
    
    def test_getOverlap_self(self): 
        '''test_getOverlap_self: the overlap between the same set of sequences must be 100%'''
        reads1, reads2, clones1, clones2, stats1, stats2 = tcrRepSim.getOverlap(self.samseqs, self.samseqs, [0, 0.01, 10, 20, 100])
        
        self.assertTrue( clones1[0] == len(self.samseqs) ) #number of clones/ uniq sequences 
        for i, c1 in enumerate( clones1 ): #each cutoff
            #every stats of sample 2 is exactly the same with sample 1
            self.assertTrue( c1 == clones2[i] ) 
            self.assertTrue( reads1[i] == reads2[i] )
            for k in stats1:
                self.assertTrue( stats1[k][i] == stats2[k][i] )
            
            #The overlaping stats is equal the the original stats:
            self.assertTrue( stats1["oclones"][i] == c1 )
            self.assertTrue( stats1["oreads"][i] == reads1[i] )

            #Number of reads must be >= number of clones
            self.assertTrue( reads1[i] >= clones1[i] )
            self.assertTrue( stats1["oreads"][i] >= stats1["oclones"][i] )

        #at 100% cutoff, numClones must be 0 if total numClones > 1
        if clones1[0] > 1:
            self.assertTrue( clones1[4] == 0 )

        #Number of clones decreases as the cutoff increases
        for i in xrange( 4 ):
            self.assertTrue( clones1[i] >= clones1[i + 1] )

    def test_getOverlap_zeroOverlap(self):
        '''Two repertoires with zero overlap should return zero overlap'''
        samseqs1 = self.samseqs
        
        #create another repertoire with totally different length
        seqs2 = tcrRepSim.getUniqSeqRep(10, self.aaUsage[11], self.codonUsage)
        tcrRepSim.getRep( seqs2, 50, None )
        samseqs2 = tcrRepSim.sampling( seqs2, 20 )
        reads1, reads2, clones1, clones2, stats1, stats2 = tcrRepSim.getOverlap(samseqs1, samseqs2, [0, 0.01, 10, 20, 100])
        for i in xrange(5):
            self.assertTrue( stats1["oclones"][i] == 0 )
            self.assertTrue( stats2["oclones"][i] == 0 )
            self.assertTrue( stats1["oreads"][i] == 0 )
            self.assertTrue( stats2["oreads"][i] == 0 )

class TestSampling( unittest.TestCase ):
    aaUsageFile = "../data/tcrRepSimTest/aaUsage.txt"
    codonUsageFile = "../data/tcrRepSimTest/codonUsage.txt"
    def test_sampling(self): 
        '''sampling should return the desired number of sequences '''
        aaUsage = tcrRepSim.readAaUsage(self.aaUsageFile)
        codonUsage = tcrRepSim.readCodonUsage(self.codonUsageFile)
        seqs = tcrRepSim.getUniqSeqRep(10, aaUsage[9], codonUsage)
        tcrRepSim.getRep( seqs, 50, None )

        samseqs = tcrRepSim.sampling( seqs, 20 )
        total = sum([s.count for s in samseqs])

        self.assertTrue( len(samseqs) <= 10 )
        self.assertTrue( total == 20 )

class TestGetRep( unittest.TestCase ):
    aaUsageFile = "../data/tcrRepSimTest/aaUsage.txt"
    codonUsageFile = "../data/tcrRepSimTest/codonUsage.txt"
    def test_getRep(self): 
        '''getRep should return the desired number of sequences '''
        aaUsage = tcrRepSim.readAaUsage(self.aaUsageFile)
        codonUsage = tcrRepSim.readCodonUsage(self.codonUsageFile)
        seqs = tcrRepSim.getUniqSeqRep(10, aaUsage[9], codonUsage)
        tcrRepSim.getRep( seqs, 50, None )
        total = sum([s.count for s in seqs])
        self.assertTrue( total == 50 )

class TestGetUniqSeqRep( unittest.TestCase ):
    aaUsageFile = "../data/tcrRepSimTest/aaUsage.txt"
    codonUsageFile = "../data/tcrRepSimTest/codonUsage.txt"
    def test_getUniqSeqRep(self): 
        '''getUniqSeqRep should return the desired number of uniq sequences '''
        aaUsage = tcrRepSim.readAaUsage(self.aaUsageFile)
        codonUsage = tcrRepSim.readCodonUsage(self.codonUsageFile)
        seqs = tcrRepSim.getUniqSeqRep(100, aaUsage[9], codonUsage)
        self.assertTrue( len(seqs) == 100)
        visited = []
        for s in seqs:
            self.assertTrue( s.nuc not in visited )
            visited .append(s.nuc)

class TestUtils( unittest.TestCase ):
    def test_getAaSeq(self):
        '''getAaSeq '''
        aaUsage = ['AAAB', 'CD']
        a2ef = {'AC': 0.75*0.5, 'AD':0.75*0.5, 'BC':0.25*0.5, 'BD':0.25*0.5}
        a2f = {}
        for i in xrange(10000):
            a = tcrRepSim.getAaSeq(aaUsage)
            if a in a2f:
                a2f[a] += 1
            else:
                a2f[a] = 1
        for a in a2f:
            self.assertTrue( a in a2ef )
            self.assertTrue( abs( float(a2f[a])/10000 - a2ef[a] ) < 0.1 )

    def test_aa2nuc(self):
        '''aa2nuc '''
        aa2codons = {'A':['abc','abc','abc','bcd'], 'B':['abb', 'bbb']}
        testaa = 'AB'
        nuc2expFreq = {}
        for c1 in aa2codons['A']:
            for c2 in aa2codons['B']:
                nuc = c1 + c2
                if nuc in nuc2expFreq:
                    nuc2expFreq[nuc] += 1
                else:
                    nuc2expFreq[nuc] = 1
        total = sum( nuc2expFreq.values() )
        for nuc in nuc2expFreq:
            nuc2expFreq[nuc] = float(nuc2expFreq[nuc])/total

        n2f = {}
        for i in xrange(10000):
            nucseq = tcrRepSim.aa2nuc(testaa, aa2codons)
            if nucseq not in n2f:
                n2f[nucseq] = 1
            else:
                n2f[nucseq] += 1
        for n in n2f:
            n2f[n] = float(n2f[n])/10000
            self.assertTrue( n in nuc2expFreq )
            self.assertTrue( abs(n2f[n] - nuc2expFreq[n]) < 0.1 )

class TestReadCodonUsage( unittest.TestCase ):
    goodEx = "../data/tcrRepSimTest/codonUsage.txt"

    def test_readCodonUsage_goodInput(self):
        '''readCodonUsage: usage of codon should add up to 1 for each a.a'''
        aa2codons = tcrRepSim.readCodonUsage(self.goodEx)
        self.assertTrue( len( aa2codons.keys() ) == 21 )
        for aa, codons in aa2codons.iteritems():
            self.assertTrue( abs(len(codons) - 100) <= 1 )
            for c in codons:
                for letter in c:
                    self.assertTrue( letter in "ACTG" )
        codonsF = aa2codons['F']
        c2c = {}
        for c in codonsF:
            if c in c2c:
                c2c[c] += 1
            else:
                c2c[c] = 1
        self.assertTrue( c2c['TTT'] == 45 )
        self.assertTrue( c2c['TTC'] == 55 )

class TestReadAaUsage( unittest.TestCase ):
    goodEx = "../data/tcrRepSimTest/aaUsage-good.txt"
    wrongFormatEx = "../data/tcrRepSimTest/aaUsage-wrongFormat.txt"
    incompleteEx = "../data/tcrRepSimTest/aaUsage-incomplete.txt"

    def test_readAaUsage_goodInput(self):
        '''readAaUsage should return correct known usage'''
        len2aaUsage = tcrRepSim.readAaUsage( self.goodEx )
        lens = [2, 4]
        knownUsage = {2:[{'C':9999.0/10000, 'A':1.0/10000}, {'F':3.0/3}], 4:[ {'C':100.0/100}, {'A':79.0/100, 'R': 1.0/100, 'S':20.0/100}, {'A':500.0/1000, 'C':200.0/1000, 'G':300.0/1000},{'F':1}] }
        
        for l in lens:
            self.assertTrue( l in len2aaUsage )
            self.assertTrue( len( len2aaUsage[l] ) == l )
            aaUsage = len2aaUsage[l]
            known = knownUsage[l]
            for i in xrange(l):
                s = aaUsage[i]
                freqs = known[i]
                a2c = {}
                for a in s:
                    if a in a2c:
                        a2c[a] += 1
                    else:
                        a2c[a] = 1
                for a in a2c:
                    a2c[a] = float( a2c[a] )/len(s)
                    self.assertTrue( a2c[a] == freqs[a] )

    def test_readAaUsage_badInput(self):
        '''readAaUsage should raise errors if the input was bad'''
        self.assertRaises( tcrRepSim.FileFormatError, tcrRepSim.readAaUsage, self.wrongFormatEx )
        self.assertRaises( tcrRepSim.FileFormatError, tcrRepSim.readAaUsage, self.incompleteEx )

class TestUniAaUsage( unittest.TestCase ):
    def test_uniAaUsage(self):
        '''getUniAaUsage should return dictionary of aaUsage for each input length'''
        lens = ['a', 1, 3]
        self.assertRaises( ValueError, tcrRepSim.getUniAaUsage, lens )
        
        aas = "ARNDCEQGHILKFPSTWYV"
        lens = [3, 8, 1]
        len2aaUsage = tcrRepSim.getUniAaUsage(lens)
        for l in len2aaUsage:
            self.assertTrue( len( len2aaUsage[l] ) == l )
            for s in len2aaUsage[l]:
                self.assertTrue( s == aas )

class TestReadCdr3LenDist( unittest.TestCase ):
    goodEx = "../data/tcrRepSimTest/cdr3lenDist.txt"
    good2Ex = "../data/tcrRepSimTest/cdr3lenDist-2.txt"
    wrongFormatEx = "../data/tcrRepSimTest/cdr3lenDist-wrongFormat.txt"
    wrongFormat2Ex = "../data/tcrRepSimTest/cdr3lenDist-wrongFormat2.txt"
    
    def test_readCdr3LenDist_goodInput(self):
        '''Good input should return a list of positive floats whose sum = 1'''
        len2freq = tcrRepSim.readCdr3LenDist(self.goodEx)
        for l, f in len2freq.iteritems():
            self.assertTrue( isinstance(l, int) and isinstance(f, float) )
            self.assertTrue( f >= 0 and f <= 1 )
        self.assertTrue( abs( 1 - sum(len2freq.values()) ) < 0.0001 )

        len2freq = tcrRepSim.readCdr3LenDist(self.good2Ex)
        for l, f in len2freq.iteritems():
            self.assertTrue( isinstance(l, int) and isinstance(f, float) )
            self.assertTrue( f >= 0 and f <= 1 )
        self.assertTrue( abs( 1 - sum(len2freq.values()) ) < 0.001 )

    def test_readCdr3LenDist_wrongFormat(self):
        self.assertRaises( ValueError, tcrRepSim.readCdr3LenDist, self.wrongFormatEx)
        self.assertRaises( tcrRepSim.FileFormatError, tcrRepSim.readCdr3LenDist, self.wrongFormat2Ex )

class TestTopFreqs( unittest.TestCase ):
    goodEx = "../data/tcrRepSimTest/topFreqs-good.txt"
    outOfRangeEx = "../data/tcrRepSimTest/topFreqs-bad-outOfRange.txt"
    wrongFormatEx = "../data/tcrRepSimTest/topFreqs-bad-wrongFormat.txt"
    negEx = "../data/tcrRepSimTest/topFreqs-bad-neg.txt"

    def test_topFreqs_goodInput(self):
        '''Good input should return a list of positive floats whose sum < 1 - (clones-numFreqs)/seqs'''
        freqs = tcrRepSim.getTopFreqs(self.goodEx, 100, 210)
        for f in freqs:
            self.assertTrue( isinstance(f, float) )
            self.assertTrue( f > 0 )
        self.assertTrue( sum(freqs) < 1 - float( (100-len(freqs)) )/210 )

    def test_topFreqs_negInput(self):
        '''Input with negative values should raise ValueError exception'''
        self.assertRaises( ValueError, tcrRepSim.getTopFreqs, self.negEx, 100, 210)

    def test_topFreqs_outOfRangeInput(self):
        '''Input with out of range values should raise exception'''
        self.assertRaises( tcrRepSim.TopFreqError, tcrRepSim.getTopFreqs, self.outOfRangeEx, 100, 210)

    def test_topFreqs_wrongFormatInput(self):
         '''Wrong format input should raise exception'''
         self.assertRaises( ValueError, tcrRepSim.getTopFreqs, self.wrongFormatEx, 100, 210)

class TestGetUniCdr3LenDist( unittest.TestCase ):
    def test_uniCdr3LenDist(self): 
        '''getUniCdr3LenDist should return a dictionary of len to frequencies. All freqs should be equal and add up to 1. Len goes from 9 - 23'''
        len2freq = tcrRepSim.getUniCdr3LenDist()
        self.assertTrue( len(len2freq.keys()) == 23 - 9 + 1 )
        expectedFreq = 1.0/(23 - 9 + 1)
        for l, f in len2freq.iteritems():
            self.assertTrue( l >= 9 and l <= 23 )
            self.assertTrue( f == expectedFreq )
        self.assertTrue( abs(1 - sum(len2freq.values())) < 0.0001)

if __name__ == '__main__':
    unittest.main()
