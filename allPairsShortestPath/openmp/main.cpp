//
// Created by xiaohe on 4/3/16.
//
#include "main.h"
//#include <limits>
#include "parseCommandLine.h"
#include <omp.h>

using namespace std;

char *allPairsShortestPath(wghEdgeArray<intT> Gr, int num_threads) {
    intT n = Gr.n;
    wghEdge<intT> *edgeList = Gr.E;

    // cout << "number of nodes is " << n << ", and number of edges is " << Gr.m << "." << endl;

    wghEdge<intT> curEdge;

    //the initial shortest path table is a n*n matrix
    //so we need to analyze the original weighted multi-graph
    //to obtain the initial configuration.
    double **weightTable = new double *[n];
    double **weightTable_omp = new double *[n];

    for (int i = 0; i < n; ++i) {
        weightTable[i] = new double[n];
        weightTable_omp[i] = new double[n];

        for (int j = 0; j < n; ++j) {
            weightTable[i][j] = numeric_limits<double>::max();
            weightTable_omp[i][j] = numeric_limits<double>::max();
        }
    }

    for (int i = 0; i < Gr.m; ++i) {
        curEdge = edgeList[i];
        int u = curEdge.u;
        int v = curEdge.v;
        double weight = curEdge.weight;

        // cout << "edge " << i << " connects node " << u << " and node " << v
        // << ", with weight " << weight << endl;

        if (weightTable[u][v] > weight){
            weightTable[u][v] = weight;
            weightTable_omp[u][v] = weight;
        }
    }

    //print the init config
    // cout << "Init config is " << endl;
    // printMatrix(weightTable, n);

    //TODO: the sequential algorithm
    double t0 = omp_get_wtime();

    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                double viaK = weightTable[i][k] + weightTable[k][j];
                if (weightTable[i][j] > viaK)
                    weightTable[i][j] = viaK;
            }
        }
    }

    double t1 = omp_get_wtime();
    printf("seq: Time = %.5lf\n", t1 - t0);


    // parallel
    omp_set_num_threads(num_threads);
    t0 = omp_get_wtime();

    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
        #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                double viaK = weightTable_omp[i][k] + weightTable_omp[k][j];
                if (weightTable_omp[i][j] > viaK)
                    weightTable_omp[i][j] = viaK;
            }
        }
    }

    t1 = omp_get_wtime();
    printf("par: Time = %.5lf\n", t1 - t0);

    // cout << "The all pairs shortest path table is:" << endl;
    // printMatrix(weightTable, n);//the final config.
    // printf("\n");
    // printMatrix(weightTable_omp, n);//the final config.

    bool break_flag = false;
    for (int i = 0; i < n && !break_flag; ++i) {
        for (int j = 0; j < n && !break_flag; ++j) {
            if(weightTable[i][j] - weightTable_omp[i][j] != 0){
                printf("error!\n");
                // break;
                break_flag = true;
            }
        }

    }

    //clean
    for (int k = 0; k < n; ++k) {
        delete[] weightTable[k];
        delete[] weightTable_omp[k];
    }

    delete[] weightTable;
    delete[] weightTable_omp;
}

int main(int argc, char **argv) {
    // commandLine P(argc, argv, "-o <outFile>");
    // char *iFile = P.getArgument(0);
    // char *oFile = P.getOptionValue("-o");

    char *iFile = argv[1];
    int num_threads = atoi(argv[2]);

    wghEdgeArray<intT> G = benchIO::readWghEdgeArrayFromFile<intT>(iFile);

    allPairsShortestPath(G, num_threads);
}

