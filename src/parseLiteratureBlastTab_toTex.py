#!/usr/bin/env python2.6

"""
Oct 5 2012
cp from parseLiteratureBlast_toTex.py
Parser Max's Blast output tab-delimited file into a table of clones that have hits matches certain keywords
Filtering criteria includes:
    1/ Match score (E.g minPositives = 90%)
    2/ Keywords (E.g: autologous keywords)
    3/ Min sequence length (to avoid short sequences which are easier to be false negatives)
    4/ Min number of samples 

nknguyen at soe ucsc edu
Sep 1
"""

import os, sys, re
import immunoseq.lib.immunoseqLib as iseqlib

class Record():
    def __init__(self, line):
    #length  identities      positives       eVal    score   tFile   tId     tLen    qFile   qId     qLen    tAln    match   qAln    authors journal year    title   abstract        fulltextUrl
        items = line.strip().split('\t')
        self.length = int(items[0])
        self.identities = int(items[1])
        self.positives = int(items[2])
        self.db = items[5].split("/")[-2]
        self.qId = items[9]
        self.qLen = int(items[10])
        self.tAln = items[11]
        self.match = items[12]
        self.qAln = items[13]
        #self.authors = items[14]
        self.title = items[17]
        self.abstract = items[18]

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
    if len(items) == 4:
        sizes = items[-2].split(',')
    else:
        sizes = [1 for s in samples]
    sample2size = {}
    
    for i, s in enumerate(samples):
        sample2size[s] = int(sizes[i])
    #print sample2size
    genes = items[1].split(',,')[0].split('|')
    vs = genes[1].split(',')
    js = genes[2].split(',')
    
    vs = ','.join([v.lstrip("TRBV") for v in vs])
    js = ','.join([j.lstrip("TRBJ") for j in js])
    ds = ''
    if len(genes) > 3:
        ds = genes[3].split(',')
        ds = ','.join([d.lstrip("TRBD") for d in ds])
    return vs, js, ds, sample2size

def parsePaperInfo(paperStr):
    #gnl|BL ORD ID|14836 AJ224294|Silins S.L.    Submitted (12-FEB-1998) to the EMBL/GenBank/DDBJ databases. Silins S.L., Queensland Institute of Medical Research, The Bancroft Centre, 300 Herston Road, Brisbane, AUSTRALIA 4029   Silins S.L., Cross S.M., Krauer K.G., Moss D.J., Schmidt C.W., Misko I.S.  "A functional link for major TCR expansions in healthy adults caused by persistent EBV-infection"  J. Clin. Invest. 102(8):1551-1558(1998).|
    #items = paperStr.split('"') 
    diseases = ["Rheumatic Heart Disease", "Autoimmune", "Ankylosing Spondylitis", "Rheumatoid Arthritis", "Reactive Arthritis", "Multiple Sclerosis", "Psoriatic Arthritis", "Spondyloarthropathy", 'Lupus', "Diabetes", "Vitiligo", "haematopoietic stem cell transplantation"]
    disease2short = {"Rheumatic Heart Disease": "RHD", "Ankylosing Spondylitis":"AS", "Rheumatoid Arthritis":"RA", "Reactive Arthritis":"ReA", "Multiple Sclerosis":"MS", "Psoriatic Arthritis":"PA", "Spondyloarthropathy":"AS, SpA", 'Lupus':'SLE', 'Vitiligo': 'V', "haematopoietic stem cell transplantation": "AS"} 
    #title = items[1].lower()
    title = paperStr.lower()
    matchDiseases = []
    for d in diseases:
        if d.lower() in title:
            if d in disease2short:
                matchDiseases.append( disease2short[d] )
            else:
                matchDiseases.append(d)
    
    #return items[1]
    if len(matchDiseases) > 0:
        return ', '.join( matchDiseases )
    else:
        return paperStr
        #return items[1]

def checkNumSamples(title):
    items = title.split(';')
    samples = items[0].lstrip('>').split(',')
    sample2patient={'SBC8': 'B', 'asBD':'B', 'asBR':'B', 'adaptBDdraw2':'B', 'asBDdraw2':'B', 
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

def checkNumPatients(title, group2sample2host, minPatientCount, minControlCount, maxPatientCount, maxControlCount):
    items = title.split(';')
    samples = items[0].lstrip('>').split(',')
    
    if len(items) == 3:
        sizes = [1 for s in samples]
    else:
        sizes = [int(size) for size in items[-2].split(',')]
    #controls = ["B", "adaptBD", "asBD", "20", "as20D", "adapt20D", "adaptBDdraw2", 'asBDdraw2']
    #sample2patient={'SBC8': 'B', 'asBD':'B', 'asBR':'B', 'adaptBDdraw2':'B', 'asBDdraw2':'B',
    #                'SBC7': '1', 'as1D':'1', 'as1R':'1',
    #                'adapt11D': '11', 'as11R': '11'}
    
    sample2patient = {}
    sample2control = {}
    if 'patient' in group2sample2host:
        sample2patient = group2sample2host['patient']
    if 'control' in group2sample2host:
        sample2control = group2sample2host['control']
    
    patients = []
    controls = []
    numPatientOutofrange = 0
    numControlOutofrange = 0
    for i, s in enumerate(samples):
        size = sizes[i]
        if s in sample2patient and s not in patients:
            if size >= minPatientCount and size <= maxPatientCount:
                patients.append(s)
            else:
                numPatientOutofrange += 1
        elif s in sample2control and s not in controls:
            if size >= minControlCount and size <= maxControlCount:
                controls.append(s)
            else:
                numControlOutofrange += 1

    return len(patients), len(controls), numPatientOutofrange, numControlOutofrange

def readInfile(options, group2sample2host): 
    #infile, minPositives, minNumsams, minLen, minNumpatients, minNumcontrols, minPatientCount, minControlCount, group2sample2host):
    #clones, clone2hits = readNcbiXml(options.infile, options.minPos, options.minNumSamples, options.minLen, options.minNumPatients, options.minNumControls, options.minPatientCount, options.minControlCount, group2sample2host)
    #Read in blast-output file:
    f = open(options.infile, 'r')

    clone2hits = {} #key = cloneName, val = [ (list of hit papers, identity) ]
    clones = []

    f.readline()
    for line in f:
        #length  identities      positives       eVal    score   tFile   tId     tLen    qFile   qId     qLen    tAln    match   qAln    authors journal year    title   abstract        fulltextUrl
        record = Record(line)
        clone = record.qId
        numsams = checkNumSamples(clone)
        if numsams < options.minNumSamples:
            continue
        numpatients, numcontrols, numPoutofrange, numCoutofrange = checkNumPatients(clone, group2sample2host, options.minPatientCount, options.minControlCount, options.maxPatientCount, options.maxControlCount)
        if numpatients < options.minNumPatients or numcontrols < options.minNumControls or numPoutofrange > options.maxPatientOutofrange or numCoutofrange > options.maxControlOutofrange:
            #print numpatients
            #print numcontrols
            continue
        if record.qLen < options.minLen or record.length < record.qLen:
            continue

        if clone not in clones:
            clones.append(clone)
        #if float( record.positives )/record.qLen < options.minPos or float(record.identities)/record.qLen < options.minIden:
        if float( record.positives )/record.length < options.minPos or float(record.identities)/record.length < options.minIden:
            continue
        if clone not in clone2hits:
            clone2hits[ clone ] = [ record ]
        else:
            clone2hits[ clone ].append( record )
    f.close()
    return clones, clone2hits

def getPc(count, total):
    if total == 0:
        return 0
    return 100.0*count/total

def getDefaultKeywords():
    #autoimmuneKeywords= ['haematopoietic stem cell transplantation']
    autoimmuneKeywords=['arthritis', 'ankylosing', 'spondy', 'autoreactive', 'autoantigen', 'reactive arthritis', 'rheumatoid arthritis', 'multiple sclerosis', 'self', 'cross-reactive', 'mimicry', 'synovial', 'crohn', 'psoriasis', 'inflammatory bowel disease', 'ibd', 'ulcerative colitis', 'uveitis']
    #autoimmuneKeywords=['arthritis', 'ankylosing', 'spondy', 'autoreactive', 'autoantigen', 'reactive arthritis', 'rheumatoid arthritis', 'multiple sclerosis', 'self', 'cross-reactive', 'mimicry', 'synovial', 'crohn', 'psoriasis', 'inflammatory bowel disease', 'ibd', 'ulcerative colitis', 'uveitis', 'haematopoietic stem cell transplantation']
    b27Keywords=['b27']
    pathogenKeywords=['chlamydia', 'salmonella', 'yersinia', 'shigella', 'campylobacter', 'vipr1', 'ebv', 'epstein-barr', 'lmp2']
    miceKeywords=['mice', 'murine', 'mouse']
    group2keywords = {'autoimmune': autoimmuneKeywords, 'b27': b27Keywords, 'pathogen': pathogenKeywords, 'mice':miceKeywords}
    return group2keywords

def printTab(clones, clone2hits, group2keywords, options, outbasename): 
    outfile = "%s.txt" % outbasename
    fh = open(outfile, 'w')
    fh.write("#MinPositives: %f; MinIdentities: %f, MinNumberOfSamples: %d; MinLen: %d\n" %(options.minPos, options.minIden, options.minNumSamples, options.minLen))
    fh.write("#Keywords:\n")
    for g, k in group2keywords.iteritems():
        fh.write("#\t%s:\t%s\n" %(g, ','.join(k)))
    
    cutoff = 1
    numAuto = 0 #number of clones with at least one hit passed cutoff and matches one of the autoimmuneKeywords
    numB27 = 0
    numPathogen = 0

    db2title2count = {'aai':{}, 'imgt':{}, 'pmc':{}, 'elsevier':{}}
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
                if checkKeywords(hit.title, group2keywords['mice']) or (options.abstract and checkKeywords(hit.abstract, group2keywords['mice'])): #mice studies, ignore
                    continue
                #Check to see if any keywords matched:
                if not matchAuto and ( checkKeywords(hit.title, group2keywords['autoimmune']) or (options.abstract and checkKeywords(hit.abstract, group2keywords['autoimmune'])) ):
                    matchAuto = True
                if not matchB27 and ( checkKeywords(hit.title, group2keywords['b27']) or (options.abstract and checkKeywords(hit.abstract, group2keywords['b27'])) ):
                    matchB27 = True
                if not matchPathogen and ( checkKeywords(hit.title, group2keywords['pathogen']) or (options.abstract and checkKeywords(hit.abstract, group2keywords['pathogen'])) ):
                    matchPathogen = True
                
            if matchAuto or matchB27 or matchPathogen:
                fh.write("\n>%s\n" %clone)
                for i, hit in enumerate(hits):
                    hasKeyword = False
                    if checkKeywords(hit.title, group2keywords['mice']) or (options.abstract and checkKeywords(hit.abstract, group2keywords['mice'])): #mice studies, ignore
                        continue
                    for g, keywords in group2keywords.iteritems():
                        if checkKeywords(hit.title, keywords) or (options.abstract and checkKeywords(hit.abstract, keywords)):
                            hasKeyword = True
                            break
                    if hasKeyword:
                        fh.write("\t%d/ %s; db=%s\n" %(i, hit.title, hit.db))
                        fh.write("\t\t%s\n" %hit.qAln)
                        fh.write("\t\t%s\n" %hit.match)
                        fh.write("\t\t%s\n" %hit.tAln)
                        if hit.title not in db2title2count[hit.db]:
                            db2title2count[hit.db][hit.title] = 1
                        else:
                            db2title2count[hit.db][hit.title] += 1

                    #PRINT LENGTH OF THE CLONE:
                    #if i == 0 and matchAuto:
                    #    sys.stdout.write("%d\n" % hit.qLen)
                
            if matchAuto:
                numAuto += 1
            elif matchB27:
                numB27 += 1
            elif matchPathogen:
                numPathogen += 1

    for db, t2c in db2title2count.iteritems():
        print "%s: %d" %(db, len(t2c.keys()))

    total = len(clones)
    numhits = len(clone2hits)
    fh.write("\n### Summary ###\n")
    fh.write("Total\tNumHits\t% hits/total\tnumAuto\t% auto/total\t% auto/hits\tnumB27\t% b27/total\t% b27/hits\tnumPathogen\t% pathogen/total\t% pathogen/hits\n")
    fh.write("%d\t%d\t%f\t%d\t%f\t%f\t%d\t%f\t%f\t%d\t%f\t%f\n" %(total, numhits, getPc(numhits, total), numAuto, getPc(numAuto, total), getPc(numAuto, numhits), numB27, getPc(numB27, total), getPc(numB27, numhits), numPathogen, getPc(numPathogen, total), getPc(numPathogen, numhits)) )
    fh.close()

####### LATEX TABLE ########
def myTabHeader(f, samples):
    f.write("\\begin{sidewaystable}\n")
    #f.write("\\begin{table}\n")
    f.write("\\centering\n")
    f.write("\\scalebox{0.9}{%\n")
    #f.write("\\begin{tabular}{c|c|c|%s|c|c|c}\n" %( "|".join(["c" for s in samples]) ) )
    f.write("\\begin{tabular}{l|l|l|%s|l|l|l}\n" %( "|".join(["l" for s in samples]) ) )
    #f.write(" \\multicolumn{3}{c|}{Clones} & \\multicolumn{%d}{c|}{Samples} & \\multicolumn{3}{c}{Hits} \\\\\n" %(len(samples)) )
    f.write(" \\multicolumn{3}{c|}{\\textbf{Clones}} & \\multicolumn{%d}{c|}{\\textbf{Samples}} & \\multicolumn{3}{c}{\\textbf{Hits}} \\\\\n" %(len(samples)) )
    #f.write("\\cline{2-%d}\n" %( len(colnames)*2 + 1 ))
    f.write("\\hline\n")
    #f.write("V & CDR3 & J & %s & CDR3 & Alignment & Disease \\\\\n" %(" & ".join(samples)))
    f.write("\\textbf{V} & \\textbf{CDR3} & \\textbf{J} & \\textbf{%s} & \\textbf{CDR3} & \\textbf{Alignment} & \\textbf{Disease} \\\\\n" %("} & \\textbf{".join(samples)))
    f.write("\\hline\n")
    
def tab(f, clones, clone2hits, group2keywords, options, samples):
    for clone in clones:
        vs, js, ds, sample2size = getCloneInfo(clone)
        if clone in clone2hits:
            hits = clone2hits[clone]
            hitsWithKeyword = [] #list of hits that have at least 1 keyword
            for hit in hits:
                if checkKeywords(hit.title, group2keywords['mice']) or (options.abstract and checkKeywords(hit.abstract, group2keywords['mice'])): #mice studies, ignore
                    continue
                for g, kw in group2keywords.iteritems():
                    if g == 'b27' or g == 'pathogen':
                        continue
                    if checkKeywords(hit.title, kw) or (options.abstract and checkKeywords(hit.abstract, kw)):
                        hitsWithKeyword.append(hit)
                        break
            if len(hitsWithKeyword) == 0: #no hit with keyword
                continue
            
            seq = hits[0].qAln
            numrow = len(hitsWithKeyword)
            
            #First line
            f.write("\\multirow{%d}{*}{%s} & \\multirow{%d}{*}{%s} & \\multirow{%d}{*}{%s} & " %(numrow, vs, numrow, seq, numrow, js) ) #Write V, CDR3, J
            for s in samples:
                name = iseqlib.properName2name(s)
                if name in sample2size:
                    count = sample2size[name]
                    f.write("\\multirow{%d}{*}{%d} & " % (numrow, count))
                else:
                    f.write("\\multirow{%d}{*}{} & " % (numrow))
            f.write("%s & %s & %s \\\\\n " %( hitsWithKeyword[0].tAln, hitsWithKeyword[0].match, parsePaperInfo(hitsWithKeyword[0].title) ))
            
            #Other hits:
            for i in xrange(1, numrow):
                f.write("\\cline{%d-%d}\n" %(3 + len(samples) + 1, 3 + len(samples) + 3))
                f.write(" &"*( 3 + len(samples) ) )
                h = hitsWithKeyword[i]
                f.write( "%s & %s & %s \\\\\n" %(h.tAln, h.match, parsePaperInfo(h.title)) )
            f.write("\\hline\n")

def printTexTab(clones, clone2hits, group2keywords, options, outbasename):
    outfile = "%s.tex" %outbasename
    f = open(outfile, 'w')
    iseqlib.writeDocumentStart(f)
    samples = ['AS1', 'AS2', 'AS3', 'AS4', 'AS5', 'H1', 'H2']
    myTabHeader(f, samples)
    tab(f, clones, clone2hits, group2keywords, options, samples)
    label = ''
    captionStr = ''
    #iseqlib.tableCloser(f, captionStr, label)
    iseqlib.sidewaystableCloser(f, captionStr, label)
    iseqlib.writeDocumentEnd(f)
    f.close()

####### LATEX TABLE FORMAT 0 ===============
def myTabHeader0(f):
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
    
def tab0(f, clones, clone2hits, group2keywords, options):
    for clone in clones:
        vs, js, ds, sample2size = getCloneInfo(clone)
        if clone in clone2hits:
            hits = clone2hits[clone]
            hitsWithKeyword = [] #list of hits that have at least 1 keyword
            for hit in hits:
                if checkKeywords(hit.title, group2keywords['mice']) or (options.abstract and checkKeywords(hit.abstract, group2keywords['mice'])): #mice studies, ignore
                    continue
                for g, kw in group2keywords.iteritems():
                    if g == 'b27' or g == 'pathogen':
                        continue
                    if checkKeywords(hit.title, kw) or (options.abstract and checkKeywords(hit.abstract, kw)):
                        hitsWithKeyword.append(hit)
                        break
            if len(hitsWithKeyword) == 0: #no hit with keyword
                continue
            
            seq = hits[0].qAln
            samples = sorted( [iseqlib.properName(s) for s in sample2size.keys()] )
            numrow = max( [len(samples), len(hitsWithKeyword)] )
            
            f.write("\\multirow{%d}{*}{%s} & \\multirow{%d}{*}{%s} & \\multirow{%d}{*}{%s} & " %(numrow, vs, numrow, seq, numrow, js) ) #Write V, CDR3, J
            #f.write("%s & %d & %s & %s & %s \\\\\n " %( samples[0], sample2size[iseqlib.properName2name(samples[0])], hitsWithKeyword[0][4], hitsWithKeyword[0][3], hitsWithKeyword[0][0] )) #first row
            f.write("%s & %d & %s & %s & %s \\\\\n " %( samples[0], sample2size[iseqlib.properName2name(samples[0])], hitsWithKeyword[0].tAln, hitsWithKeyword[0].match, parsePaperInfo(hitsWithKeyword[0].title) ))
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
                    f.write( "%s & %s & %s \\\\\n" %(h.tAln, h.match, parsePaperInfo(h.title)) )
                else:
                    f.write(" & & \\\\\n")
            f.write("\\hline\n")

def printTexTab0(clones, clone2hits, group2keywords, options, outbasename):
    outfile = "%s.tex" %outbasename
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

def readSample2host(file):
    f = open(file, 'r')
    group2sample2host = {}
    for line in f:
        items = line.strip().split()
        if len(items) < 3:
            continue
        group = items[0]
        sample = items[1]
        host = items[2]
        if group not in group2sample2host:
            group2sample2host[group] = {sample:host}
        else:
            group2sample2host[group][sample] = host
    f.close()
    return group2sample2host

def addOptions(parser):
    parser.add_option('-i', '--infile', dest='infile', help='Input xml file')
    parser.add_option('-o', '--outdir', dest='outdir', help='Output directory')
    parser.add_option('-b', '--basename', dest='basename', default='hits', help='Output files basename. Default=%default')
    parser.add_option('-p', '--positive', dest='minPos', type='float', default=0.9, help='Minimum portion of positive matches. Default=%default')
    parser.add_option('-I', '--identity', dest='minIden', type='float', default=0.9, help='Minimum portion of positive matches. Default=%default')
    parser.add_option('-k', '--keywords', dest='keywords', default=None, help='Only hits matching at least one keyword are reported')
    parser.add_option('-l', '--len', dest='minLen', type='int', default=10, help='Minimum sequence length to be included in the output. Default=%default')
    parser.add_option('-s', '--samples', dest='minNumSamples', type='int', default=1, help='Minimum number of samples containing the sequence. Default=%default')
    parser.add_option('--patients', dest='minNumPatients', type='int', default=0, help='Minimum number of patients containing the sequence. Default=%default')
    parser.add_option('--controls', dest='minNumControls', type='int', default=0, help='Minimum number of controls containing the sequence. Default=%default')
    parser.add_option('--minPatientCount', dest='minPatientCount', type='int', default=1, help='Minimum size a clone must have in a patient sample to be considered as "present" in that sample. Default=%default')
    parser.add_option('--minControlCount', dest='minControlCount', type='int', default=1, help='Minimum size a clone must have in a control sample to be considered as "present" in that sample. Default=%default')
    parser.add_option('--maxPatientCount', dest='maxPatientCount', type='int', default=10000000, help='Maximun size a clone must have in a patient sample to be considered as "present" in that sample. Default=%default')
    parser.add_option('--maxControlCount', dest='maxControlCount', type='int', default=10000000, help='Maximun size a clone must have in a control sample to be considered as "present" in that sample. Default=%default')
    parser.add_option('--maxPatientOutofrange', dest='maxPatientOutofrange', type='int', default=100, help='Max number of patients with outofrange counts allowed. Default=%default')
    parser.add_option('--maxControlOutofrange', dest='maxControlOutofrange', type='int', default=100, help='Max number of controls with outofrange counts allowed. Default=%default')
    parser.add_option('--sample2host', dest='sample2host', help='Optional. File contains mapping between samples and host. Format:<Group> <sample> <host>. Ex: control asBD B ')
    parser.add_option('-a', '--abstract', dest='abstract', action='store_true', default=False, help='If specified, will include the abstract in the search of matching keywords')

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

    group2sample2host = {}
    if options.sample2host:
        group2sample2host = readSample2host(options.sample2host)
    #clones, clone2hits = readNcbiXml(options.infile, options.minPos, options.minNumSamples, options.minLen, options.minNumPatients, options.minNumControls, options.minPatientCount, options.minControlCount, group2sample2host)
    #clones, clone2hits = readNcbiXml(options, group2sample2host)
    clones, clone2hits = readInfile(options, group2sample2host)
    outbasename = os.path.join(options.outdir, options.basename)
    printTab(clones, clone2hits, group2keywords, options, outbasename)
    printTexTab(clones, clone2hits, group2keywords, options, outbasename)

if __name__ == '__main__':
    main()


