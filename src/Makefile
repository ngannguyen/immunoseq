binPath = ../bin
libPath = ../lib

progs = $(notdir $(wildcard cloneSizeDist.py adaptiveTcrParse.py tcrRepSim.py tcrRepSimTest.py getOverlapTab.py vjusage.py usagePlot.py))

.PHONY: %.all

targets = ${progs:%=${binPath}/%} 

all: ${targets}

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

clean :
	rm -f ${progs:%=${binPath}/%}
