# includes
incdir = ./include
algodir = ./algorithms

# IBM XL compiler
BG_XL = /bgsys/drivers/ppcfloor/comm/xl
BG_INC = -I$(BG_PATH)/include
BG_LIB = -L$(BG_PATH)/lib
BG_ALG = -I$(BG_PATH)/algorithms

# compilers/flags
compiler = g++ -O3 
pcompiler =mpic++ -O3 -std=c++0x
flags = -I$(incdir) -I$(algodir) -I$(algodir)/topology
pflags = $(flags) -include mpi.h

# RPI CCI AMOS compilers/flags
qcompiler = mpic++ -O5  -qarch=qp -qtune=qp -qflag=w -qstrict -qreport
#qcompiler = mpic++ -O5 -qarch=qp -qtune=qp -qflag=w -qstrict -qprefetch=aggressive -qsimd=auto -qhot=fastmath -qinline=level=10
#qflags = $(CFLAGS) $(BG_INC) $(BG_LIB) $(LDFLAGS) $(flags)

# ONLY uncomment the following if <module load experimental/zlib> FAILS.
qflags = $(BG_INC) $(BG_LIB) $(BG_ALG) $(flags) -I/bgsys/apps/CCNI/zlib/zlib-1.2.7/include -L/bgsys/apps/CCNI/zlib/zlib-1.2.7/lib 


# dependencies
core = $(incdir)/MMSP.utility.h \
       $(incdir)/MMSP.grid.h \
       $(incdir)/MMSP.sparse.h


# the program
gg.out: gg.cpp
	$(compiler) $(flags) $< -o $@ -lz

parallel: main.cpp gg.cpp tessellate.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@ -lz -lpthread	

bgqmc: main.cpp gg.cpp tessellate.hpp $(core)
	$(qcompiler) $(qflags)   -DBGQ -DSILENT -DRAW $< -o q_MC.out
	
wrongendian: wrongendian.cpp 
	$(compiler) $< -o $@.out -lz -pthread

clean:
	rm -rf solidification.out parallel q_MC.out
