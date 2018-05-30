LIBS = -lm
CC = g++
WARNINGS = -Wall
CFLAGS = -g -std=c++11 -march=native -O3 -flto -fstrict-aliasing $(WARNINGS)
LFLAGS = -g -march=native -O3 -flto

SRCDIR = src
OBJDIR = obj
INCDIR = .

################################# for bloom filter code #################################
CXX=g++
# directory of sdsl library.so
SDSL_PATH=$(HOME)/lib
# directory of TCLAP head file
TCLAP_PATH=./include
# also requires head file of roaring, sdsl, jellyfish, tclap in ~/include
CXXFLAGS= -g -std=c++11 -pthread -Wall -O3 
CXXFLAGS+=-I$(HOME)/include -I$(TCLAP_PATH)

LD_LIB=-L $(SDSL_PATH)
LD_FLAG=-lsdsl -ldivsufsort -ldivsufsort64 

################################ end bloom filter def ##################################

.PHONY: all clean
all: gbf_lite vargeno_lite

OBJECTS = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(wildcard $(SRCDIR)/*.c))
HEADERS = $(wildcard $(INCDIR)/*.h)

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(HEADERS)
	$(CXX) $(CFLAGS) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	$(CXX) -c -o $@ $^ $(CXXFLAGS) -I$(INCDIR)

.PRECIOUS: $(TARGET) $(OBJECTS)

gbf_lite: $(OBJDIR)/gbf.o $(OBJDIR)/generate_bf.o $(OBJDIR)/allsome_util.o $(OBJDIR)/util.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_LIB) $(LD_FLAG)

kmerc: $(OBJDIR)/kmerc.o $(OBJDIR)/kmer_counting.o $(OBJDIR)/generate_bf.o $(OBJDIR)/allsome_util.o $(OBJDIR)/util.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_LIB) $(LD_FLAG)

vargeno_lite: $(OBJDIR)/qv.o $(OBJDIR)/generate_bf.o $(OBJDIR)/allsome_util.o $(OBJDIR)/dict_filt.o $(OBJECTS) 
	$(CC) $(LFLAGS) $(CXXFLAGS) $(OBJDIR)/qv.o $(OBJDIR)/generate_bf.o $(OBJDIR)/allsome_util.o $(OBJDIR)/dict_filt.o $(OBJECTS) $(LIBS) -o $@ $(LD_LIB) $(LD_FLAG)

clean:
	-rm -f $(OBJDIR)/*.o
	-rm -f gbf_lite
	-rm -f vargeno_lite
