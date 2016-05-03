SHELL=/bin/bash

include gppDefs

# required files
GLOBAL_REQUIRE = graph.h
LOCAL_REQUIRE =
OBJS = APSP.o

CC=g++
CCFLAGS= -D ENABLE_DEBUG -pthread -std=c++11 -g


all: clean main.exe

main.exe: *.cpp
	$(CC) $(CCFLAGS) *.cpp -o main.exe

clean:
	rm -f *.exe $(TARGETS)
