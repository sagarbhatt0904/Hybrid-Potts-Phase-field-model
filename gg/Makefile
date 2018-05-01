# Makefile
# GNU makefile for dendritic solidification code
# Questions/comments to trevor.keller@gmail.com (Trevor Keller)

# includes
incdir = ./include
incdir = ./algorithms


# compilers/flags
compiler = g++ -O3
pcompiler = mpic++ -O3
flags = -I $(incdir)
pflags = $(flags) -include mpi.h

# the program
gg.out: gg.cpp
	$(compiler) $(flags) $< -o $@ -lz

parallel: gg.cpp
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -rf solidification.out parallel
