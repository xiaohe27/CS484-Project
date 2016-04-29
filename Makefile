clean: 
	-rm *.o *~ 

test: 
	cd allPairsShortestPath/sequential; make all; ./main.exe

all: clean test
