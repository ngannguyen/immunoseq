binPath = ../bin
libPath = ../lib

progs = $(notdir $(wildcard parseLiteratureBlastTab_toTex.py clusterBlastXml.py clusterType_blastclust.py clusterType.py filterBlastXml.py blastSummaryTab.py sampleInfoTab.py geneCorrTab.py geneusedTab.py splitFaByVJ.py parseLiteratureBlast_toTex.py parseLiteratureBlast.py overlapWithSampling.py diversityPlot.py diversityPlot.R repertoireSurvey.py cloneFilter.py similarity.py rnaseq-genesfpkm.py pairwiseOverlapTab.py pairwiseScatter.py nt2aa.py overlapNt2aa.py overlap2.py lendist.py aausage.py clonesizeDist2.py immunoSeq.py sampling.py faCollapse.py cloneSizeTab.py overlap.py getOverlapSequences.py cloneSizeDist.py adaptiveTcrParse.py tcrRepSim.py tcrRepSimTest.py getOverlapTab.py vjusage.py usagePlot.py))

.PHONY: %.all

targets = ${progs:%=${binPath}/%} 

all: ${targets}

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

clean :
	rm -f ${progs:%=${binPath}/%}
