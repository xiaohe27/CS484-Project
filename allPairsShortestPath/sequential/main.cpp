//
// Created by xiaohe on 4/3/16.
//

//#include "parseCommandLine.h"
#include "APSP.h"

//using namespace std;

int main (int argc, char** argv) {
//    commandLine P(argc,argv,"-o <outFile>");
//    char* iFile = P.getArgument(0);
//    char* oFile = P.getOptionValue("-o");

    char * iFile = argv[0];
    graph<intT> G = benchIO::readGraphFromFile<intT>(iFile);

//    allPairsShortestPath(G);
}

