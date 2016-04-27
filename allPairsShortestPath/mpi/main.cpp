//
// Created by xiaohe on 4/3/16.
//
#include "main.h"
//#include <limits>
#include "parseCommandLine.h"
#include <mpi.h>

#define BILLION 1000000000L
#define ind(i, j) (n * i) + j

using namespace std;

int main(int argc, char **argv) {

    int n, *sendcounts, *displs, start;

    double *local, *rowK_space, *rowK;

    int has_rowK, row_per_pro;

    double diff;
    double seqTime, parallelTime;

    double *weightTable;
    double *weightTable_mpi;

    MPI_Init(&argc, &argv);

    int rank, num_pro;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_pro);

    if (rank == 0) {

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
        weightTable_mpi = new double[n * n];
        for (int i = 0; i < n; ++i) {
            // weightTable[i] = new double[n];
            // weightTable_mpi[i] = new double[n];
            for (int j = 0; j < n; ++j) {
                weightTable[ind(i, j)] = numeric_limits<double>::max();
                weightTable_mpi[ind(i, j)] = numeric_limits<double>::max();
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
                weightTable_mpi[ind(u, v)] = weight;
            }
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
    //                 printf(" %1.2f",i,weightTable_mpi[ind(i,j)]);
    //         printf("\n");
    //     }
    //     printf("\n");
    // }

    struct timespec start_ser, end_ser;


    // sequential version
    clock_gettime(CLOCK_MONOTONIC, &start_ser);

    if (rank == 0) {
        for (int k = 0; k < n; ++k) {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < n; ++i) {
                    double viaK = weightTable[ind(i, k)] + weightTable[ind(k, j)];
                    if (weightTable[ind(i, j)] > viaK)
                        weightTable[ind(i, j)] = viaK;
                }
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &end_ser);

    diff = (double) ((double) BILLION * (end_ser.tv_sec - start_ser.tv_sec) + end_ser.tv_nsec - start_ser.tv_nsec) /
           1000000;

    if (rank == 0) {
        seqTime = diff;
        printf("Time taken for sequential version = %f milliseconds\n", (double) diff);
    }

    // parallel version
    clock_gettime(CLOCK_MONOTONIC, &start_ser);

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);

    // printf("here rank %d has n %d\n", rank, n);

    row_per_pro = n / num_pro;
    start = row_per_pro * rank;
    // local  = (double*) malloc(sizeof(double) * row_per_pro * n);
    // rowK_space = (double *) malloc(sizeof(double) * n);

    local = new double[row_per_pro * n];
    rowK_space = new double[n];


    if (rank == 0) {
        // sendcounts = (int*) malloc(sizeof(int) * num_pro);
        // displs = (int*) malloc(sizeof(int) * num_pro);

        sendcounts = new int[num_pro];
        displs = new int[num_pro];
        for (int i = 0; i < num_pro; ++i) {
            sendcounts[i] = row_per_pro * n;
            displs[i] = i * row_per_pro * n;
        }

        MPI_Scatterv(weightTable_mpi, sendcounts, displs, MPI_DOUBLE, local, row_per_pro * n, MPI_DOUBLE, 0,
                     MPI_COMM_WORLD);
    }
    else
        MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, local, row_per_pro * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    // for(int r = 0; r < num_pro; r++) {
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     if (r == rank) {
    //         printf("rank: %d\n", rank);
    //         for(int i = 0; i < row_per_pro; i++){
    //             printf("row %d: ", i);
    //             for(int j = 0; j < n; j++)
    //                 if(local[ind(i,j)] > 10)
    //                     printf(" 0.00");
    //                 else
    //                     printf(" %1.2f",i,local[ind(i,j)]);
    //             printf("\n");
    //         }
    //     }
    // }




    // do algo
    for (int k = 0; k < n; ++k) {
        has_rowK = k / row_per_pro;
        if (has_rowK >= num_pro)
            has_rowK = num_pro - 1;

        if (has_rowK == rank)
            rowK = &local[(k - start) * n];
        else
            rowK = rowK_space;

        MPI_Bcast(rowK, n, MPI_DOUBLE, has_rowK, MPI_COMM_WORLD);

        // MPI_Barrier(MPI_COMM_WORLD);

        for (int i = 0; i < row_per_pro; ++i)
            for (int j = 0; j < n; ++j)
                if (local[ind(i, k)] + rowK[j] < local[ind(i, j)])
                    local[ind(i, j)] = local[ind(i, k)] + rowK[j];
    }
    free(rowK_space);

    MPI_Barrier(MPI_COMM_WORLD);


    if (rank == 0)
        MPI_Gatherv(local, row_per_pro * n, MPI_DOUBLE, weightTable_mpi, sendcounts, displs, MPI_DOUBLE, 0,
                    MPI_COMM_WORLD);
    else
        MPI_Gatherv(local, row_per_pro * n, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(local);

    clock_gettime(CLOCK_MONOTONIC, &end_ser);

    diff = (double) ((double) BILLION * (end_ser.tv_sec - start_ser.tv_sec) + end_ser.tv_nsec - start_ser.tv_nsec) /
           1000000;

    if (rank == 0) {
        parallelTime = diff;
        printf("Time taken for parallel   version = %f milliseconds\n", (double) diff);

        double iso_efficiency = seqTime / (num_pro * parallelTime);
        printf("Iso efficiency is %f.", iso_efficiency);
    }

    if (rank == 0) {
        bool break_flag = false;
        for (int i = 0; i < n && !break_flag; ++i) {
            for (int j = 0; j < n && !break_flag; ++j) {
                if (weightTable[ind(i, j)] - weightTable_mpi[ind(i, j)] != 0) {
                    printf("error!\n");
                    // break;
                    break_flag = true;
                }
            }
        }
    }

    MPI_Finalize();
}

