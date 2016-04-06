//
// Created by xiaohe on 4/3/16.
//
#include "main.h"

#include "parseCommandLine.h"

using namespace std;

char* allPairsShortestPath(wghEdgeArray<intT> Gr) {
    intT n = Gr.n;
    wghEdge<intT> * edgeList = Gr.E;

    cout << "number of nodes is " << n << ", and number of edges is " << Gr.m << "." << endl;

    wghEdge<intT> curEdge;
    for (int i = 0; i < Gr.m; ++i) {
        curEdge = edgeList[i];
        cout << "edge " << i << " connects node " << curEdge.u << " and node " << curEdge.v
                << ", with weight " << curEdge.weight << endl;
        
    }
}

int main (int argc, char** argv) {
    commandLine P(argc,argv,"-o <outFile>");
    char* iFile = P.getArgument(0);
    char* oFile = P.getOptionValue("-o");

    wghEdgeArray<intT> G = benchIO::readWghEdgeArrayFromFile<intT>(iFile);

    allPairsShortestPath(G);
}

