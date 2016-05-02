//
// Created by xiaohe on 5/1/16.
//
#include "main.h"
//#include <limits>
#include "parseCommandLine.h"

#define BILLION 1000000000L
#define ind(i, j) (n * i) + j

using namespace std;

int main(int argc, char **argv) {

    int n;

    double diff;
    double seqTime = 0.0;

    double *weightTable;

    commandLine P(argc, argv, "-o <outFile>");
    char *iFile = P.getArgument(0);
    char *oFile = P.getOptionValue("-o");

    wghEdgeArray<intT> Gr = benchIO::readWghEdgeArrayFromFile<intT>(iFile);

// double **weightTable = allPairsShortestPath(G);

// initialization for table
// intT n = Gr.n;
    n = Gr.n;
    wghEdge<intT> *edgeList = Gr.E;
    wghEdge<intT> curEdge;
    weightTable = new double[n * n];
    for (int i = 0; i < n; ++i) {
// weightTable[i] = new double[n];
// weightTable_mpi[i] = new double[n];
        for (int j = 0; j < n; ++j) {
            weightTable[ind(i, j)] = numeric_limits<double>::max();
        }
    }
    for (int i = 0; i < Gr.m; ++i) {
        curEdge = edgeList[i];
        int u = curEdge.u;
        int v = curEdge.v;
        double weight = curEdge.weight;

// cout << "edge " << i << " connects node " << u << " and node " << v
// << ", with weight " << weight << endl;

        if (weightTable[ind(u, v)] > weight) {
            weightTable[ind(u, v)] = weight;
        }
    }


// if (rank == 0) {
//     printf("original table \n", rank);
//     for(int i = 0; i < n; i++){
//         printf("row %d: ", i);
//         for(int j = 0; j < n; j++)
//             if(weightTable_mpi[ind(i,j)] > 10)
//                 printf(" 0.00");
//             else
//                 printf(" %1.2f",weightTable_mpi[ind(i,j)]);
//         printf("\n");
//     }
//     printf("\n");
// }

    struct timespec start_ser, end_ser;


// sequential version
    clock_gettime(CLOCK_MONOTONIC, &start_ser);

    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                double viaK = weightTable[ind(i, k)] + weightTable[ind(k, j)];
                if (weightTable[ind(i, j)] > viaK)
                    weightTable[ind(i, j)] = viaK;
            }
        }
    }


// if (rank == 0) {
//     printf("original table \n", rank);
//     for(int i = 0; i < n; i++){
//         printf("row %d: ", i);
//         for(int j = 0; j < n; j++)
//             if(weightTable[ind(i,j)] > 10)
//                 printf(" 0.00");
//             else
//                 printf(" %1.2f",weightTable[ind(i,j)]);
//         printf("\n");
//     }
//     printf("\n");
// }

    clock_gettime(CLOCK_MONOTONIC, &end_ser);

    diff = (double) ((double) BILLION * (end_ser.tv_sec - start_ser.tv_sec) + end_ser.tv_nsec - start_ser.tv_nsec) /
           1000000;

    seqTime = diff;
    printf("Time taken for sequential version = %f milliseconds\n", (double) diff);


}
