#!/usr/bin/env python

"""
Jan 13 2012
nknguyen@soe.ucsc.edu
Extract (productive a.a) sequences that present in both adaptiveRepertoire & iRepertoire with >= cutoffCount and >= minPercentageOfTotalSampleReads
Input: Input directory that contains fasta files of the samples with header format: >id;size=#
        sample2patient.txt file with format: SampleName\tPatient
       Output directory
Output: fasta files of sequences of each patient that present in at least 2 samples with >= minCutoff and minPercentage
"""

import os, re, sys
from optparse import OptionParser
from sonLib.bioio import system

def addseq(aa2headers, seq, header):
    if seq != '' and header != '':
        if seq in aa2headers and header not in aa2headers[seq]:
            aa2headers[seq].append(header)
        else:
            aa2headers[seq] = [header]

def readFaFile(file):
    f = open(file, 'r')
    aa2headers = {}
    header = ''
    seq = ''
    for line in f:
        line = line.strip()
        if line[0] == '>':
            addseq(aa2headers, seq, header)
            seq = ''
            header = line.lstrip('>')
        else:
            seq = line
    addseq(aa2headers, seq, header)
    f.close()
    return aa2headers

def getTotalCount( aa2headers ):
    total = 0
    for aa in aa2headers:
        total += getCount( aa2headers[aa] )
    return total

def getCount( headers ):
    count = 0
    for h in headers:
        items = h.split('size=')
        count += int(items[1])
    return count

def readSamples(indir, sample2patient, minCount):
    patient2aa2samples = {} 
    aa2samples = {}
    sample2aaNum = {}
    for file in os.listdir(indir):
        fileItems = file.split('.')
        if len(fileItems) == 1:
            continue
        path = os.path.join(indir, file)
        if os.path.isdir(path):
            continue
        sample = file.split('.')[0]
        patient = sample2patient[sample]
        aa2headers = readFaFile(path)
        totalCount = getTotalCount( aa2headers )
        if totalCount == 0:
            sys.stderr.write('File %s has no sequences\n' %file)
            continue
        #sample2aaNum[ sample ] = len( aa2headers.keys() )
        sample2aaNum[ sample ] = [0,0] #(aaNumbers, totalCount)

        for aa, headers in aa2headers.iteritems():
            count = getCount(headers)
            #percent = 100.0*count/totalCount
            if count >= minCount:
                sample2aaNum[ sample ][0] += 1
                sample2aaNum[ sample ][1] += count
                if patient not in patient2aa2samples:
                    patient2aa2samples[patient] = {aa:[(sample, count)]}
                elif aa not in patient2aa2samples[patient]:
                    patient2aa2samples[patient][aa] = [(sample, count)]
                else:
                    patient2aa2samples[patient][aa].append( (sample, count) )
                if aa not in aa2samples:
                    aa2samples[aa] = [ (sample, count) ]
                else:
                    aa2samples[aa].append( (sample, count) )
    
    return patient2aa2samples, aa2samples, sample2aaNum

def filter(patient2aa2samples, patient2samples, outdir, required, minCount, minPercentage, sample2aaNum, mode):
    numPassRequired = 1
    if mode == 1:
        numPassRequired = required #If mode 1, only take sequences that are present with frequencies >= minPercentage in at least "required" number of samples. 

    for patient in patient2aa2samples:
        f = open( os.path.join(outdir, "%s.fa" % patient), 'w')
        numSamples = len(patient2samples[patient])
        if numSamples < required:
            required = numSamples
        aa2samples = patient2aa2samples[patient]
        for aa in aa2samples:
            samples = aa2samples[aa]
            numPass = 0
            numPresent = 0
            for s in samples:
                count = s[1]
                percent = 100.0*count/sample2aaNum[s[0]][1]
                if count >= minCount and percent >= minPercentage:
                    numPass += 1
                if count > 0:
                    numPresent += 1

            if numPresent >= required and numPass >= numPassRequired: #The clone is present in at least 'required' number of samples, and has >= minCount, minPercentage in at least one sample
                header = "%s" %( '|'.join( [ "%s-%d" %(s[0], s[1]) for s in samples ] ) )
                total = sum([ s[1] for s in samples])
                f.write(">%s;size=%d\n" %(header, total))
                f.write("%s\n" %aa)
        f.close()

def readSample2patient(file):
    f = open(file, 'r')
    sample2patient = {}
    patient2samples = {}
    for line in f:
        items = line.strip().split()
        if len(items) < 2:
            continue
        s = items[0]
        p = items[1]
        sample2patient[ s ] = p
        if p not in patient2samples:
            patient2samples[p] = [s]
        else:
            patient2samples[p].append(s)
    f.close()
    return sample2patient, patient2samples

def checkOptions(parser, args, options):
    if options.indir == '':
        parser.error("Input directory was not specified\n")
    if options.outdir == '':
        parser.error("Output directory was not specified\n")
    if not os.path.exists(options.sample2patient):
        parser.error("sample2patient File %s does not exists\n" %options.sample2patient)
    if options.mode != 1 and options.mode !=2:
        parser.error("Invalid mode %d. Please select either 1 or 2\n" %options.mode)
    options.minPercentage = [float(p) for p in options.minPercentage.strip(',').split(',')]

def initOptions(parser):
    parser.add_option('-i', '--indir', dest='indir', help="Required argument. Input directory")
    parser.add_option('-o', '--outdir', dest='outdir', help="Required argument. Output directory")
    parser.add_option('-r', '--minRequiredSamples', dest='minSamples', default=1, type='int', help='Minimun number of samples for each patient the clone must be present in to pass filtering. Default = 1')
    parser.add_option('-c', '--minCount', dest='minCount', default=1, type='int', help='Minimun clone count to be included. Default=%default')
    parser.add_option('-p', '--minPercentage', dest='minPercentage', default="0.01,0.05", type='string', help='Comma separated list of percentage cutoffs. Clone must have at least this percentage of total reads to be included. Default=%deafult')
    parser.add_option('-s', '--sample2patient', dest='sample2patient', help="File that maps each sample to corresponding patient. Format:[Sample Patient]")
    parser.add_option('-m', '--mode', dest='mode', type='int', default=2, help="Mode 1: for each sample, only look at clones that pass the minCount and minPercentage. Mode 2: clones must pass minCount and minPercentage in only at least one sample. Default=%default")

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage = usage)
    initOptions(parser)
    options, args = parser.parse_args()
    checkOptions(parser, args, options)

    sample2patient, patient2samples = readSample2patient( options.sample2patient )
    sys.stderr.write("read sample2patient info\n")

    patient2aa2samples, aa2samples, sample2aaNum = readSamples(options.indir, sample2patient, options.minCount)
    sys.stderr.write("read sample sequences\n")
    
    for minPercent in options.minPercentage:
        outdir = os.path.join(options.outdir, "%.4f" %minPercent )
        system("mkdir -p %s" %outdir)
        filter(patient2aa2samples, patient2samples, outdir, options.minSamples, options.minCount, minPercent, sample2aaNum, options.mode)
        sys.stderr.write("Filtered minPercentage %f\n" %minPercent)
    
if __name__ == '__main__':
    main()

