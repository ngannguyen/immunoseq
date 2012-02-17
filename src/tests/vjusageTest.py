#!/usr/bin/env python

import immunoseq.src.vjusage as vju
import unittest

class TestVjusage( unittest.TestCase ):
    
    def test_vjusage(self):
        indir = "/hive/users/nknguyen/reconGit/immunoseq/data/vjusageTest/seqs"
        samples = vju.readFiles(indir) 
        for s in samples:
            s.getVJusage()
        
        #V usage:
        s2v2c = { 'sample1':{'TRBV1-1':35, 'TRBV1-2':15, 'TRBV2':20, 'TRBV3':20, 'TRBV4':10, 'TRBV5':0}, 'sample2':{'TRBV1-1':95, 'TRBV1-2':15, 'TRBV2':20, 'TRBV3':20, 'TRBV4':0, 'TRBV5':50} }
        s2v2uc = { 'sample1':{'TRBV1-1':1, 'TRBV1-2':0, 'TRBV2':0, 'TRBV3':0, 'TRBV4':1, 'TRBV5':0}, 'sample2':{'TRBV1-1':2, 'TRBV1-2':0, 'TRBV2':0, 'TRBV3':0, 'TRBV4':0, 'TRBV5':1} }
        knownVs = 'TRBV1-1,TRBV1-2,TRBV2,TRBV3,TRBV4,TRBV5'.split(',')
        
        vgenes = vju.getUnionGeneList(samples, 'v')
        for v in knownVs:
            self.assertTrue( v in vgenes )
        for s in samples:
            if s.name not in ['sample1', 'sample2']:
                print s.name
                continue
            v2c = s.usage['v']
            for v in v2c:
                count = int( v2c[v][0] )
                uniq = int( v2c[v][1] )
                knowncount = s2v2c[s.name][v]
                knownuniq = s2v2uc[s.name][v]
                
                self.assertTrue( count == knowncount )
                self.assertTrue( uniq == knownuniq )
        
        #VJ usage:
        vjs = vju.getUnionGeneList(samples, 'vj')
        vj2c = samples[0].usage['vj']
        count = int( vj2c['TRBV1-1|TRBJ1-1'][0] )
        self.assertTrue( count == 35 )
        
        uniqcount = int( vj2c['TRBV1-1|TRBJ1-1'][1] )
        self.assertTrue( uniqcount == 1 )
        
        count = int( vj2c['TRBV2|TRBJ2'][0] )
        self.assertTrue( count == 10 )
        
        uniqcount = int( vj2c['TRBV2|TRBJ2'][1] )
        self.assertTrue( uniqcount == 0 )
        
        

if __name__ == '__main__':
    unittest.main()
