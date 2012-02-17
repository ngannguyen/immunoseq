#!/usr/bin/env python

#nknguyen at soe ucsc edu
#Aug 3 2011
#Immuno-seq pipeline

'''
Immuno seq pipepline
Input: reads fastq files
Functions:

'''

import os, sys, re, time
from optparse import OptionParser
import xml.etree.ElementTree as ET

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import getTempDirectory
from sonLib.bioio import setLogLevel

from immunoseq.lib import *

class Setup( Target ):
    """Setting up jobs
    """
    def __init__(self, options):
        Target.__init__(self, time=0.00025)
        self.options = options
    
    def run(self):
        setLogLevel("DEBUG")
        infiles = os.listdir(self.options.indir)
        for file in infiles:
            if re.search('\.fastq$', file) or re.search('\.fq$', file):
                sample = os.path.basename(file).split('.')[0]
                self.addChildTarget( Sample(sample, self.options) )


class Sample( Target ):
    """Run different filtering for the sample
    """
    def __init__(self, name, options):
        self.name = name
        self.options = options
    
    def run(self):
         




def checkOptions( parser, options, args ):
    if not options.indir:
        parser.error("Input directory is required, none was given.\n")
    if not options.outdir:
        parser.error("Output directory is required, none was given.\n")
    if not os.path.exists( options.indir ):
        parser.error("Input directory %s does not exist.\n" %(options.indir) )
    if not os.path.exists( options.outdir ):
        system( "mkdir %s\n" %(options.outdir) )

def addOptions( parser ):
    parser.add_option("-i", "--indir", dest = 'indir', help="Required argument. Input directory that contains input fastq files (each represents a sample)\n")
    parser.add_option("-o", "--outdir", dest = 'outdir', help="Required argument. Output directory\n")
    parser.add_option("-v", "--vref", dest='vref', help="Fasta file contains V reference sequences\n")
    parser.add_option("-j", "--jref", dest='jref', help="Fasta file contains J reference sequences\n")
    parser.add_option("-l", "--lenghCutoff", dest'lcutoff', default=50, type="int", help="Read length cutoff. This is the first filtering step. Reads with length less than this cutoff are filtered out. Default = 50.\n")
    parser.add_option("-q", "--qualCutoff", dest='qcutoff', default=30, type="int", help="Quality cutoff. This is the second filtering step. Reads with quality less than this cutoff are filtered out. Default = 30.\n")
    parser.add_option("-p", "--qualPCutoff", dest='qpcutoff', default=95, type="int", help="Percentage of read length with quality >= qualCutoff for the reads to be considered passed quality filtering. Default = 95, meaning a read of length 100 bases must have 95 bases or more with quality >= qualCutoff.\n")
    parser.add_option("-c", "--coverageCutoff", dest='cov', default=1, type='int', help="Coverage cutoff. Reads with less than this many coverage are filtered out. Default = 1, meaning no coverage filtering")
    parser.add_option("-d", "--dcutoff", dest'dcutoff', default=100, type='int', help="Unique sequences with corresponding counts (number of reads) are sorted in descending order. Sequences whose reads account for this 'dcutoff' percent of total reads are kept. The rest of the reads are filtered out. Default = 100, meaning no filtering.")

def main():
    parser = OptionParser()
    Stack.addJobTreeOptions( parser )

    addOptions( parser )

    options, args = parser.parse_args()
    checkOptions( parser, options, args )

    i = Stack( Setup(options) ).startJobTree( options )
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" %i)

if __name__ == "__main__":
    from immunoseq import *
    main()




