CC=g++
CFLAGS=-O3 -c -Wall
DFLAGS=-g -Wall
LDFLAGS= -lz -lm

#KSWSOURCE=ksw.c 
ALGNSOURCES=Aln.cpp Del_local_aln.cpp mk_hash.cpp LCtrl_option.cpp 
INDSOURCES=Index.cpp mk_hash.cpp readfl.cpp Ctrl_option.cpp

INDOBJECTS=$(INDSOURCES:.cpp=.o)
ALGNOBJECTS=$(ALGNSOURCES:.cpp=.o)

INDEXER=indexer
ALIGNER=aligner

all: $(INDSOURCES) $(ALGNSOURCES) $(ALIGNER) $(INDEXER)
	
$(ALIGNER): $(ALGNOBJECTS)
	$(CC)  $(ALGNOBJECTS) -o $@ $(LDFLAGS)

$(INDEXER): $(INDOBJECTS)
	$(CC)  $(INDOBJECTS)  -o $@ $(LDFLAGS)


.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
.c.o:
	$(CC) $(CFLAGS) $< -o $@
debug:
	$(CC)   $(DFLAGS)  $(ALGNSOURCES)   $(LDFLAGS)

clean:
	rm -rf *.o $(ALIGNER) $(INDEXER) a.out

 
