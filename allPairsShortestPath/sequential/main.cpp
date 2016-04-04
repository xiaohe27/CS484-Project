//
// Created by xiaohe on 4/3/16.
//
#include "main.h"

//#include "parseCommandLine.h"

using namespace std;

char* allPairsShortestPath(wghEdgeArray<intT> Gr) {
    intT n = Gr.n;
    wghEdge<intT> * edgeList = Gr.E;

    cout << "n is " << n << ", and m is " << Gr.m << "." << endl;
}

int main (int argc, char** argv) {
//    commandLine P(argc,argv,"-o <outFile>");
//    char* iFile = P.getArgument(0);
//    char* oFile = P.getOptionValue("-o");

    char * iFile = argv[1];
    wghEdgeArray<intT> G = benchIO::readWghEdgeArrayFromFile<intT>(iFile);

    allPairsShortestPath(G);
}

