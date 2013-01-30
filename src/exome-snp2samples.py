#!/usr/bin/env python

#Thu Aug  2 17:14:42 PDT 2012
#nknguyen soe ucsc edu
#
#Input: directory of mutation files, one file/sample
#       group2samples (Format:<group> <sample1,sample2,...>), eg: group=patients/controls
#       
#Output: returns SNPs that are present in at least N number of patients & at most M number of controls
#

import os, sys, re, time
import immunoseq.lib.immunoseqLib as iseqlib


def addOptions(parser):
    parser.add_option('-')

def main():
    parser = iseqlib.initOptions3()
    addOptions(parser)
    options, args = parser.parse_args()
    iseqlib.checkOptions3(parser, args, options)


if __name__ == '__main__':
    main()

