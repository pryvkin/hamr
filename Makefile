SAMTOOLS_DIR=~/src/samtools

all: rnapileup filter_pileup rnapileup2mismatchbed

rnapileup: rnapileup.cpp
	g++ -O2 -L $(SAMTOOLS_DIR) -I $(SAMTOOLS_DIR) $? -o $@ -lbam -lz

filter_pileup: filter_pileup.cpp
	g++ -O2 -L $(SAMTOOLS_DIR) -I $(SAMTOOLS_DIR) $? -o $@ -lbam -lz

rnapileup2mismatchbed: rnapileup2mismatchbed.cpp
	g++ -O2 rnapileup2mismatchbed.cpp $? -o $@

clean:
	rm rnapileup filter_pileup rnapileup2mismatchbed


