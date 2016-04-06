//
// Created by xiaohe on 4/3/16.
//
#include "main.h"
//#include <limits>
#include "parseCommandLine.h"

using namespace std;

char *allPairsShortestPath(wghEdgeArray<intT> Gr) {
    intT n = Gr.n;
    wghEdge<intT> *edgeList = Gr.E;

    cout << "number of nodes is " << n << ", and number of edges is " << Gr.m << "." << endl;

    wghEdge<intT> curEdge;

    //the initial shortest path table is a n*n matrix
    //so we need to analyze the original weighted multi-graph
    //to obtain the initial configuration.
    double **weightTable = new double *[n];

    for (int i = 0; i < n; ++i) {
        weightTable[i] = new double[n];

        for (int j = 0; j < n; ++j) {
            weightTable[i][j] = numeric_limits<double>::max();
        }
    }

    for (int i = 0; i < Gr.m; ++i) {
        curEdge = edgeList[i];
        int u = curEdge.u;
        int v = curEdge.v;
        double weight = curEdge.weight;

        cout << "edge " << i << " connects node " << u << " and node " << v
        << ", with weight " << weight << endl;

        if (weightTable[u][v] > weight)
            weightTable[u][v] = weight;
    }

    //print the init config
    cout << "Init config is " << endl;
    printMatrix(weightTable, n);

    //TODO: the sequential algorithm
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                double viaK = weightTable[i][k] + weightTable[k][j];
                if (weightTable[i][j] > viaK)
                    weightTable[i][j] = viaK;
            }
        }
    }

    cout << "The all pairs shortest path table is:" << endl;
    printMatrix(weightTable, n);//the final config.

    //clean
    for (int k = 0; k < n; ++k) {
        delete[] weightTable[k];
    }

    delete[] weightTable;
}

int main(int argc, char **argv) {
    commandLine P(argc, argv, "-o <outFile>");
    char *iFile = P.getArgument(0);
    char *oFile = P.getOptionValue("-o");

    wghEdgeArray<intT> G = benchIO::readWghEdgeArrayFromFile<intT>(iFile);

    allPairsShortestPath(G);
}

