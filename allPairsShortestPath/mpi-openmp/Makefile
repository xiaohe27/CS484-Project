SHELL=/bin/bash

include gppDefs

# required files
GLOBAL_REQUIRE = graph.h
LOCAL_REQUIRE =
OBJS = APSP.o

CC=mpicxx
CCFLAGS= -D ENABLE_DEBUG -std=c++11 -g


all: clean main.exe

main.exe:
	$(CC) -DITERS=1 -DNUM_THREADS=2 $(CCFLAGS) main.cpp -o main.exe

clean:
	rm -f *.exe $(TARGETS)

seq: clean seq.exe

seq.exe:
	$(CC) -DITERS=5 $(CCFLAGS) sequential.cpp -o seq.exe


seqMPI: clean seqMPI.exe

seqMPI.exe:
	$(CC) -DITERS=5 $(CCFLAGS) seq_mpi.cpp -o seq_mpi.exe
