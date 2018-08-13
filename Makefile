LIBS = -lm
CC = g++
WARNINGS = -Wall -Wextra
CFLAGS = -std=c++11 -march=native -O3 -flto -fstrict-aliasing $(WARNINGS)
LFLAGS = -march=native -O3 -flto

SRCDIR = src
OBJDIR = obj
INCDIR = .

################################# for bloom filter code #################################
CXX=g++
# directory of sdsl library.so
SDSL_PATH=$(PREFIX)/lib
# also requires head file of roaring, sdsl, jellyfish, tclap in ~/include
CXXFLAGS= -std=c++11 -Wall -O3 -I$(PREFIX)/include

LD_LIB=-L $(SDSL_PATH)
LD_FLAG=-lsdsl -ldivsufsort -ldivsufsort64 

################################ end bloom filter def ##################################

.PHONY: all clean
all: gbf vargeno

OBJECTS = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(wildcard $(SRCDIR)/*.c))
HEADERS = $(wildcard $(INCDIR)/*.h)

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(HEADERS)
	$(CXX) $(CFLAGS) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	$(CXX) -c -o $@ $^ $(CXXFLAGS) -I$(INCDIR)

.PRECIOUS: $(TARGET) $(OBJECTS)

gbf: $(OBJDIR)/gbf.o $(OBJDIR)/generate_bf.o $(OBJDIR)/allsome_util.o $(OBJDIR)/util.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_LIB) $(LD_FLAG)

kmerc: $(OBJDIR)/kmerc.o $(OBJDIR)/kmer_counting.o $(OBJDIR)/generate_bf.o $(OBJDIR)/allsome_util.o $(OBJDIR)/util.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_LIB) $(LD_FLAG)

vargeno: $(OBJDIR)/qv.o $(OBJDIR)/generate_bf.o $(OBJDIR)/allsome_util.o $(OBJECTS) 
	$(CC) $(LFLAGS) $(CXXFLAGS) $(OBJDIR)/qv.o $(OBJDIR)/generate_bf.o $(OBJDIR)/allsome_util.o $(OBJECTS) $(LIBS) -o $@ $(LD_LIB) $(LD_FLAG)

clean:
	-rm -f $(OBJDIR)/*.o
	-rm -f gbf
	-rm -f vargeno
