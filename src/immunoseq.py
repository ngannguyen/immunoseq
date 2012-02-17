#!/usr/bin/env python

#nknguyen at soe ucsc edu
#Feb 14 2012
#Immuno-seq pipeline

'''
Immuno seq pipepline
Inputs: 
    1/ Directory contains fasta files of productive CDR3 sequences (1 file per sample/patient)  in the following format:
        >sampleName;userDefinedCloneID|Vgenes|Jgenes;size=#
        Sequence
    Vgenes and Jgenes are comma separated, size = clone size (number of reads)
    Example: 
        >SBC1;Clusters6438|TRBV7-2|TRBJ2-3;size=1
        CASSFTQGETDTQYF

    2/ File that maps sample to patient, with format: <SampleName> <PatientName>
    Example: SBC2 B

    3/ 

Outputs/Functions:
    1/ V, J, VJ usage
    2/ Clone Sizes
    3/ CDR3 length distribution
    4/ Pairwise-overlapping
    5/ 

'''

import os, sys, re, time
import xml.etree.ElementTree as ET

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import getTempDirectory
from sonLib.bioio import setLogLevel

import immunoseqLib

############## Adaptive Biotechnologies File Processing ##################
class AdaptiveBioPreProcess( Target ):
    '''Process Adaptive Biotechnologies output files (.tsv) into productive CDR3 fasta files
    '''
    def __init__(self, indir, outdir, collapse):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.outdir = outdir
        self.collapse = collapse

    def run(self):
        files = os.listdir(self.indir)
        for f in files:
            self.addChildTarget( AdaptiveBioTsvToFa(self.indir, self.outdir, self.collapse, filename) )

class AdaptiveBioTsvToFa( Target ):
    def __init__(self, indir, outdir, collapse, tsvfile):
        self.indir = indir
        self.outdir = outdir
        self.collapse = collapse
        self.tsvfile = tsvfile

    def run(self):
        samplename = self.tsvfile.rstrip(".tsv")
        localTempDir = self.getLocalTempDir()
        system("cp %s %s" %( os.path.join(self.indir, self.tsvfile), localTempDir ))
        system("adaptiveTcrParse.py %s %s" %(localTempDir, localTempDir))
        localfa = os.path.join( localTempDir, "%s.fa" %samplename )
        if not self.collapse:
            system("mv %s %s" %( localfa, self.outdir ))
        else:
            collapseFa = os.path.join( localTempDir, "%sCol.fa" %samplename )
            system("faCollapse.py %s %s" %(localfa, collapseFa))
            system("mv %s %s" %(collapseFa, self.outdir))

############## Main pipeline ###############
class Setup( Target ):
    """Setting up jobs
    """
    def __init__(self, options):
        Target.__init__(self, time=0.00025)
        self.options = options
    
    def run(self):
        #V, J, VJ usage:
        usageOutdir = os.path.join(self.options.outdir, "vjusage")
        system("mkdir -p %s" %usageOutdir)
        self.addChildTarget( VJusage(self.options.indir, usageOutdir) )
        
        #Clone sizes
        clonesizeOutdir = os.path.join( self.options.outdir, "clonesize")
        system("mkdir -p %s" %(clonesizeOutdir))
        self.addChildTarget( Clonesize(self.options.indir, clonesizeOutdir) )

        #CDR3 length distribution
        lendistOutdir = os.path.join( self.options.outdir, "lendist")
        system("mkdir -p %s" %(lendistOutdir))
        self.addChildTarget( Lendist(self.options.indir, lendistOutdir) )

        #Pairwise overlapping
        overlapOutdir = os.path.join( self.options.outdir, "overlap" )
        system("mkdir -p %s" %overlapOutdir)
        self.addChildTarget( Overlap(self.options.indir, overlapOutdir) )

class VJusage( Target ):
    '''Get V, J, VJ usage
    '''
    def __init__(self, indir, outdir):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.outdir = outdir

    def run(self):
        system("vjusage.py -i %s -o %s" %(self.indir, self.outdir))
        vjdir = os.path.join(self.outdir, "vj")
        self.addChildTarget( VJusagePlot(vjdir, self.outdir) )

class VJusagePlot( Target ):
    '''Plot V, J, VJ usage
    '''
    def __init__(self, indir, outdir):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.outdir = outdir

    def run(self):
        system("usagePlot.py -i %s -o %s" %(self.indir, self.outdir))

class Clonesize( Target ):
    def __init__(self, indir, outdir):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.outdir = outdir

    def run(self):
        system()


def checkOptions( parser, options, args ):
    if not options.indir:
        raise InputOptionError("Input directory is required, none was given.\n")
    if not os.path.isdir( options.indir ):
        raise InputOptionError("Input directory (%s) is not a directory\n" %options.indir)
    if not os.path.exists( options.outdir ):
        system( "mkdir %s\n" %(options.outdir) )

def addOptions( parser ):
    parser.add_option("-s", "--sequences", dest = 'seqdir', help="Required argument. Input directory that contains fasta files of productive CDR3 sequences (each file represents a sample)\n")
    parser.add_option("-o", "--outdir", dest = 'outdir', default='immunoseqOut', help="Required argument. Output directory. Default=%default\n")

def main():
    parser = immunoseqLib.initOptions()
    Stack.addJobTreeOptions( parser )

    addOptions( parser )

    options, args = parser.parse_args()
    checkOptions( parser, options, args )

    i = Stack( Setup(options) ).startJobTree( options )
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" %i)

if __name__ == "__main__":
    from immunoseq.src.immunoseq import *
    main()




