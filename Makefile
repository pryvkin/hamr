SAMTOOLS_DIR=samtools

#DEBUG=
DEBUG=-g

CXX = g++
CXXFLAGS = -O2 -Wall $(DEBUG) -I $(SAMTOOLS_DIR)
LFLAGS = -L $(SAMTOOLS_DIR) -lbam -lz -lpthread

PROG = hamr_cmd
SRCS = main.cpp rnapileup.cpp filter_pileup.cpp rnapileup2mismatchbed.cpp util.cpp
HDRS = hamr.h
OBJS = $(SRCS:cpp=o)

all: $(PROG)

$(PROG): $(OBJS) $(HDRS)
	$(CXX) $(CXXFLAGS) $(LFLAGS) $(OBJS) -o $@ $(LFLAGS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(PROG)


