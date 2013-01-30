#!/usr/bin/env python2.6

"""
Oct 5 2012
cp from parseLiteratureBlast.py
Parser Blast output file into a table of clones that have hits matches certain keywords
Filtering criteria includes:
    1/ Match score (E.g minPositives = 90%)
    2/ Keywords (E.g: autologous keywords)
    3/ Min sequence length (to avoid short sequences which are easier to be false negatives)
    4/ Min number of samples 

nknguyen at soe ucsc edu
Sep 1
"""

import os, sys, re
from Bio.Blast import NCBIXML
import immunoseq.lib.immunoseqLib as iseqlib

def reformatTitle(title):
    title = title.replace("_", " ")
    title = title.replace(";", " ")
    return title

def checkKeywords(title, keywords):
    t = title.lower()
    for keyword in keywords:
        if re.search(keyword, t):
            return True
    return False

def getCloneInfo(clonestr):
    #>adapt11D,adapt15D,adapt16D;111621|TRBV19|TRBJ2-1,TRBJ2-7|TRBD1-1;5,45504,1427;size=45504
    items = clonestr.lstrip('>').split(';')
    samples = items[0].split(',')
    size = int( items[-1].lstrip('size=') )
    sizes = items[-2].split(',')
    sample2size = {}
    for i, s in enumerate(samples):
        sample2size[s] = int(sizes[i])
    genes = items[1].split('|')
    vs = genes[1].split(',')
    js = genes[2].split(',')
    
    vs = ','.join([v.lstrip("TRBV") for v in vs])
    js = ','.join([j.lstrip("TRBJ") for j in js])
    if len(genes) > 3:
        ds = genes[3].split(',')
        ds = ','.join([d.lstrip("TRBD") for d in ds])
    return vs, js, ds, sample2size

def parsePaperInfo(paperStr):
    #gnl|BL ORD ID|14836 AJ224294|Silins S.L.    Submitted (12-FEB-1998) to the EMBL/GenBank/DDBJ databases. Silins S.L., Queensland Institute of Medical Research, The Bancroft Centre, 300 Herston Road, Brisbane, AUSTRALIA 4029   Silins S.L., Cross S.M., Krauer K.G., Moss D.J., Schmidt C.W., Misko I.S.  "A functional link for major TCR expansions in healthy adults caused by persistent EBV-infection"  J. Clin. Invest. 102(8):1551-1558(1998).|
    items = paperStr.split('"') 
    return items[1]

def checkNumSamples(title):
    items = title.split(';')
    samples = items[0].lstrip('>').split(',')
    sample2patient={'SBC8': 'B', 'asBD':'B', 'asBR':'B', 
                    'SBC7': '1', 'as1D':'1', 'as1R':'1',
                    'adapt11D': '11', 'as11R': '11'}
    patients = []
    for s in samples:
        if s not in sample2patient:
            patient = s
        else:
            patient = sample2patient[s]
        if patient not in patients:
            patients.append(patient)
    return len(patients)

def readNcbiXml(infile, minPositives, minNumsams, minLen):
    #Read in blast-output file:
    rh = open(infile)
    records = NCBIXML.parse( rh)

    clone2hits = {} #key = cloneName, val = [ (list of hit papers, identity) ]
    clones = []

    for record in records:
        clone = record.query
        numsams = checkNumSamples(clone)
        if numsams < minNumsams:
            continue
        if record.query_length < minLen:
            continue

        clones.append(clone)
        for aln in record.alignments:
            for hit in aln.hsps:
                if float(hit.positives)/len(hit.query) < minPositives:
                #if float(hit.identities)/len(hit.query) < minPositives:
                    continue
                if clone in clone2hits:
                    clone2hits[ clone ].append( (reformatTitle(aln.title), hit.identities, hit.query, hit.match, hit.sbjct) )
                else:
                    clone2hits[ clone ] = [ (reformatTitle(aln.title), hit.identities, hit.query, hit.match, hit.sbjct) ]
    return clones, clone2hits

def getPc(count, total):
    if total == 0:
        return 0
    return 100.0*count/total

def getDefaultKeywords():
    autoimmuneKeywords=['arthritis', 'ankylosing', 'spondy', 'autoreactive', 'autoantigen', 'reactive arthritis', 'rheumatoid arthritis', 'multiple sclerosis', 'self', 'cross-reactive', 'mimicry', 'synovial', 'crohn', 'psoriasis', 'inflammatory bowel disease', 'ibd', 'ulcerative colitis', 'uveitis']
    b27Keywords=['b27']
    pathogenKeywords=['chlamydia', 'salmonella', 'yersinia', 'shigella', 'campylobacter', 'vipr1', 'ebv', 'epstein-barr', 'lmp2']
    group2keywords = {'autoimmune': autoimmuneKeywords, 'b27': b27Keywords, 'pathogen': pathogenKeywords}
    return group2keywords

def printTab(clones, clone2hits, group2keywords, options, outdir): 
    outfile = os.path.join(outdir, 'hits.txt')
    fh = open(outfile, 'w')
    fh.write("#MinPositives: %f; MinNumberOfSamples: %d; MinLen: %d\n" %(options.minPos, options.minNumSamples, options.minLen))
    fh.write("#Keywords:\n")
    for g, k in group2keywords.iteritems():
        fh.write("#\t%s:\t%s\n" %(g, ','.join(k)))
    
    cutoff = 1
    numAuto = 0 #number of clones with at least one hit passed cutoff and matches one of the autoimmuneKeywords
    numB27 = 0
    numPathogen = 0

    for clone in clones:
        #vs, ds, js, sample2size = getCloneInfo(clone)
        items = clone.split(';')
        id = '.'.join( [items[0], items[1]] )
        
        if clone in clone2hits:
            hits = clone2hits[ clone ]
            
            matchAuto = False
            matchB27 = False
            matchPathogen = False
            for i, hit in enumerate(hits):
                #Check to see if any keywords matched:
                if not matchAuto and checkKeywords(hit[0], group2keywords['autoimmune']):
                    matchAuto = True
                if not matchB27 and checkKeywords(hit[0], group2keywords['b27']):
                    matchB27 = True
                if not matchPathogen and checkKeywords(hit[0], group2keywords['pathogen']):
                    matchPathogen = True
                
            if matchAuto or matchB27 or matchPathogen:
                fh.write("\n>%s\n" %clone)
                for i, hit in enumerate(hits):
                    hasKeyword = False
                    for g, keywords in group2keywords.iteritems():
                        if checkKeywords(hit[0], keywords):
                            hasKeyword = True
                            break
                    if hasKeyword:
                        fh.write("\t%d/ %s\n" %(i, hit[0]))
                        fh.write("\t\t%s\n" %hit[2])
                        fh.write("\t\t%s\n" %hit[3])
                        fh.write("\t\t%s\n" %hit[4])
                
            if matchAuto:
                numAuto += 1
            elif matchB27:
                numB27 += 1
            elif matchPathogen:
                numPathogen += 1

    total = len(clones)
    numhits = len(clone2hits)
    fh.write("\n### Summary ###\n")
    fh.write("Total\tNumHits\t% hits/total\tnumAuto\t% auto/total\t% auto/hits\tnumB27\t% b27/total\t% b27/hits\tnumPathogen\t% pathogen/total\t% pathogen/hits\n")
    fh.write("%d\t%d\t%f\t%d\t%f\t%f\t%d\t%f\t%f\t%d\t%f\t%f\n" %(total, numhits, getPc(numhits, total), numAuto, getPc(numAuto, total), getPc(numAuto, numhits), numB27, getPc(numB27, total), getPc(numB27, numhits), numPathogen, getPc(numPathogen, total), getPc(numPathogen, numhits)) )
    fh.close()

####### LATEX TABLE ########
def myTabHeader(f):
    #f.write("\\begin{sidewaystable}\n")
    f.write("\\begin{table}\n")
    f.write("\\centering\n")
    f.write("\\scalebox{0.4}{%\n")
    f.write("\\begin{tabular}{c|c|c|c|c|c|c|c}\n")
    f.write(" \\multicolumn{3}{c|}{Clones} & \\multicolumn{2}{c|}{Samples} & \\multicolumn{3}{c}{Hits} \\\\\n")
    #f.write("\\cline{2-%d}\n" %( len(colnames)*2 + 1 ))
    f.write("\\hline\n")
    f.write("V & CDR3 & J & Name & Size & CDR3 & Alignment & Paper \\\\\n")
    f.write("\\hline\n")
    
def tab(f, clones, clone2hits, group2keywords, options):
    for clone in clones:
        vs, js, ds, sample2size = getCloneInfo(clone)
        if clone in clone2hits:
            hits = clone2hits[clone]
            hitsWithKeyword = [] #list of hits that have at least 1 keyword
            for hit in hits:
                for g, kw in group2keywords.iteritems():
                    if g == 'b27' or g == 'pathogen':
                        continue
                    if checkKeywords(hit[0], kw):
                        hitsWithKeyword.append(hit)
                        break
            if len(hitsWithKeyword) == 0: #no hit with keyword
                continue
            
            seq = hits[0][2]
            samples = sorted( [iseqlib.properName(s) for s in sample2size.keys()] )
            numrow = max( [len(samples), len(hitsWithKeyword)] )
            
            f.write("\\multirow{%d}{*}{%s} & \\multirow{%d}{*}{%s} & \\multirow{%d}{*}{%s} & " %(numrow, vs, numrow, seq, numrow, js) ) #Write V, CDR3, J
            #f.write("%s & %d & %s & %s & %s \\\\\n " %( samples[0], sample2size[iseqlib.properName2name(samples[0])], hitsWithKeyword[0][4], hitsWithKeyword[0][3], hitsWithKeyword[0][0] )) #first row
            f.write("%s & %d & %s & %s & %s \\\\\n " %( samples[0], sample2size[iseqlib.properName2name(samples[0])], hitsWithKeyword[0][4], hitsWithKeyword[0][3], parsePaperInfo(hitsWithKeyword[0][0]) ))
            for i in xrange(1, numrow):
                f.write("\\cline{4-8}\n")
                f.write(" & & & ")
                if i < len(samples):
                    s = samples[i]
                    f.write(" %s & %d &" %(s, sample2size[iseqlib.properName2name(s)]) )
                else:
                    f.write(" & & ")
                if i < len(hitsWithKeyword):
                    h = hitsWithKeyword[i]
                    #f.write( "%s & %s & %s \\\\\n" %(h[4], h[3], h[0]) )
                    f.write( "%s & %s & %s \\\\\n" %(h[4], h[3], parsePaperInfo(h[0])) )
                else:
                    f.write(" & & \\\\\n")
            f.write("\\hline\n")

def printTexTab(clones, clone2hits, group2keywords, options, outdir):
    outfile = os.path.join(outdir, 'hits.tex')
    f = open(outfile, 'w')
    iseqlib.writeDocumentStart(f)
    myTabHeader(f)
    tab(f, clones, clone2hits, group2keywords, options)
    label = ''
    captionStr = ''
    iseqlib.tableCloser(f, captionStr, label)
    #iseqlib.sidewaystableCloser(f, captionStr, label)
    iseqlib.writeDocumentEnd(f)
    f.close()

def addOptions(parser):
    parser.add_option('-i', '--infile', dest='infile', help='Input xml file')
    parser.add_option('-o', '--outdir', dest='outdir', help='Output directory')
    parser.add_option('-p', '--positive', dest='minPos', type='float', default=0.9, help='Minimum portion of positive matches. Default=%default')
    parser.add_option('-k', '--keywords', dest='keywords', default=None, help='Only hits matching at least one keyword are reported')
    parser.add_option('-l', '--len', dest='minLen', type='int', default=10, help='Minimum sequence length to be included in the output. Default=%default')
    parser.add_option('-s', '--samples', dest='minNumSamples', type='int', default=1, help='Minimum number of samples containing the sequence. Default=%default')

def main():
    parser = iseqlib.initOptions()
    addOptions(parser)
    options, args = parser.parse_args()
    group2keywords = {} #key = keywordGroup, val = list of keywords
    if options.keywords:
        if options.keywords == '-':
            group2keywords = getDefaultKeywords()
        else:
            group2keywords, kw2group = iseqlib.readGroup2samples(options.keywords)

    clones, clone2hits = readNcbiXml(options.infile, options.minPos, options.minNumSamples, options.minLen)
    #printTab(clones, clone2hits, group2keywords, options, options.outdir)
    printTexTab(clones, clone2hits, group2keywords, options, options.outdir)

if __name__ == '__main__':
    main()


