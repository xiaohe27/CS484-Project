SHELL=/bin/bash

include gppDefs

# required files
GLOBAL_REQUIRE = graph.h
LOCAL_REQUIRE =
OBJS = APSP.o

CC=g++
CCFLAGS= -lrt -D ENABLE_DEBUG -std=c++11 -g


all: main.exe

main.exe:
	$(CC) $(CCFLAGS) *.cpp -o main.exe

clean:
	rm -f *.exe $(TARGETS)
