# $Id: makefile,v 1.5 2002/08/14 13:46:18 clark Exp $

BINDIR = bin

COMPILER = g++
CDEBUG = -g
#CXXFLAGS = $(CDEBUG)  -Wall -O3 -DNDEBUG
CXXFLAGS = -g -Wall -O3

OBJECTS = matrix.o clusters.o simplecorpus.o hmm.o twohmm.o
BINARIES = cluster_neyessen cluster_baseline cluster_neyessenmorph cluster_perplexity cluster_merge  \
	cluster_traintwohmm cluster_tagtwohmm 
%.o:	%.cpp
	$(COMPILER) -c $(CXXFLAGS) $< -o $@



all: $(BINARIES)

clean:
	rm *.o

$(BINDIR):
	mkdir -p $(BINDIR)

cluster_neyessen: neyessen.o $(OBJECTS) $(BINDIR)
	$(COMPILER) -O3 -o $(BINDIR)/cluster_neyessen neyessen.o $(OBJECTS) -lpthread

cluster_baseline: baseline.o $(OBJECTS) $(BINDIR)
	$(COMPILER) -O3 -o $(BINDIR)/cluster_baseline baseline.o $(OBJECTS) -lpthread

cluster_neyessenmorph: neyessenmorph.o  $(OBJECTS) $(BINDIR)
	$(COMPILER) -O3 -o $(BINDIR)/$@ $< $(OBJECTS) -lpthread

cluster_perplexity: perplexity.o $(OBJECTS) $(BINDIR)
	$(COMPILER) -O3 -o $(BINDIR)/cluster_perplexity perplexity.o $(OBJECTS) -lpthread

cluster_merge: merge.o $(OBJECTS) $(BINDIR)
	$(COMPILER) -O3 -o $(BINDIR)/cluster_merge merge.o $(OBJECTS) -lpthread

cluster_traintwohmm: traintwohmm.o $(OBJECTS) $(BINDIR)
	$(COMPILER) -O3 -o $(BINDIR)/cluster_traintwohmm traintwohmm.o $(OBJECTS) -lpthread

cluster_tagtwohmm: tagtwohmm.o $(OBJECTS) $(BINDIR)
	$(COMPILER) -O3 -o $(BINDIR)/cluster_tagtwohmm tagtwohmm.o $(OBJECTS) -lpthread
