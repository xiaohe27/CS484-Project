//
// Created by xiaohe on 4/3/16.
//
#include "apsp.h"
#include "graphIO.h"
using namespace std;


inline void clean(double **weightTable, int n) {
    //clean
    for (int k = 0; k < n; ++k) {
        delete[] weightTable[k];
    }

    delete[] weightTable;
}

wghEdgeArray<intT> allPairsShortestPath(wghEdgeArray<intT> Gr) {
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

    wghEdge<intT> *finalEdgeList = new wghEdge<intT>[n * n];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            finalEdgeList[i * n + j] = wghEdge<intT>(i, j, weightTable[i][j]);
        }
    }

    clean(weightTable, n);

    return wghEdgeArray<intT>(finalEdgeList, n, n * n);
}

void write_allPairsShortestPath_2_file(char *iFile, char *oFile) {
    wghEdgeArray<intT> G = benchIO::readWghEdgeArrayFromFile<intT>(iFile);

    wghEdgeArray<intT> finalG = allPairsShortestPath(G);


    if (oFile != NULL) {
        cout << "output file is " << oFile << endl;

        //the final graph can be used to check correctness of other implementations
        writeWghEdgeArrayToFile(finalG, oFile);

        finalG.del();
    } else {
        cout << "output file is null" << endl;
    }
}


//int main(){}