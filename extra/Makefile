CC=mpic++.mpich
CFLAGS=-Wall -O3 -std=c++11
DCFLAGS=-g -DDEBUG
SGPPDIR = $(HOME)/Tools/combi
GENEDIR = $(HOME)/Work/HiWi_Fault_Tolerance/Code/repos/gene_python_interface/gene_20140826
OPTDIR = $(HOME)/Work/HiWi_Fault_Tolerance/Code/repos/manager_test_ft/opt_combi_technique
LDIR=-L$(SGPPDIR)/lib/sgpp -L$(OPTDIR)/lib
LIBS=-lsgppcombigrid -lboost_serialization -lsgppbase -llp_opt
INC=-I$(SGPPDIR)/src/sgpp -I$(GENEDIR)/src -I$(OPTDIR)/src

SGPP_LIB=$(SGPPDIR)/lib/sgpp/libsgppcombigrid.so
OPT_LIB =$(OPTDIR)/lib/liblp_opt.so

all: manager 

debug: manager_debug

timing.o: timing.c timing.h
	$(CC) $(CFLAGS) -c timing.c

manager: manager.cpp $(SGPP_LIB) $(OPT_LIB)
	$(CC) $(CFLAGS) $(INC) $(LDIR) -o $@ $< $(LIBS)
	
errorCalc: errorCalc.cpp $(SGPP_LIB) timing.o
	$(CC) $(CFLAGS) $(INC) $(LDIR) -o $@ $< timing.o $(LIBS) -fopenmp 
	
convertFGtoCP: convertFGtoCP.cpp $(SGPP_LIB) timing.o
	$(CC) $(CFLAGS) $(INC) $(LDIR) -o $@ $< timing.o $(LIBS) 

evCalc: evCalc.cpp $(SGPP_LIB) timing.o
	$(CC) $(CFLAGS) $(INC) $(LDIR) -o $@ $< timing.o $(LIBS) -fopenmp

clean:
	rm -rf manager testLoadModel testGeneInstance
	rm -rf *.o

