UNAME_S := $(shell uname -s)

CPPFLAGS := $(CPPFLAGS) -DNDEBUG -pthread

CXX = g++
CXXFLAGS := $(CXXFLAGS) -O3 -fopenmp -ffast-math -msse -mfpmath=sse

LDFLAGS := $(LDFLAGS) -O3 -fopenmp -pthread -lpthread
ifeq ($(UNAME_S),Linux)
    LDFLAGS += -static
endif

HDRS = \
  alpha.h \
  bipal.h \
  cmds.h \
  countsort.h \
  fastaseqsource.h \
  fastq.h \
  fastqseqsource.h \
  fileseqsource.h \
  filetype.h \
  gobuff.h \
  linereader.h \
  lockobj.h \
  lockobjs.h \
  myopts.h \
  myutils.h \
  obj.h \
  objmgr.h \
  objtype.h \
  objtypes.h \
  omplock.h \
  seqdb.h \
  seqinfo.h \
  seqsource.h \
  sort.h \
  translator.h \
  xlat.h \

OBJS = \
  o/alpha.o \
  o/bicodons.o \
  o/bipal.o \
  o/dsdn.o \
  o/endlessprob.o \
  o/evolve_neutral.o \
  o/fastaseqsource.o \
  o/fastq.o \
  o/fastqseqsource.o \
  o/fileseqsource.o \
  o/filetype.o \
  o/freq3.o \
  o/freq3afa.o \
  o/getorfs_serratus.o \
  o/kmer.o \
  o/mutable.o \
  o/phasim.o \
  o/rotate_motif.o \
  o/translate.o \
  o/codes.o \
  o/linereader.o \
  o/lockobj.o \
  o/objmgr.o \
  o/reorf_main.o \
  o/myutils.o \
  o/seqdb.o \
  o/seqdbfromfasta.o \
  o/seqinfo.o \
  o/seqsource.o \
  o/translatede.o \
  o/translator.o \
  o/wobblemers.o \
  o/wobblesim.o \

.PHONY: clean

o/reorf : o/ $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $@
	strip -d o/reorf

o/ :
	mkdir -p o/

o/%.o : %.cpp $(HDRS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	-rm -rf o/
