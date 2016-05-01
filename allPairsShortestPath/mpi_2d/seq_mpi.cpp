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

    int n, *sendcounts, *displs;

    double *local, *rowK, *colK, *out_order_buf;

    int sub_matrix_size, row_per_pro;

    double diff;
    double seqTime, parallelTime = 0.0;

    double *weightTable;
    double *weightTable_mpi;

    MPI_Status status;
    MPI_Comm comm_row;
    MPI_Comm comm_col;

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
    //                 printf(" %1.2f",weightTable_mpi[ind(i,j)]);
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

    if (rank == 0) {
        seqTime = diff;
        printf("Time taken for sequential version = %f milliseconds\n", (double) diff);
    }

    for (int x = 0; x < ITERS; ++x) {

        // parallel version
        clock_gettime(CLOCK_MONOTONIC, &start_ser);

        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // printf("here rank %d has n %d\n", rank, n);

        sub_matrix_size = n / (int) sqrt(num_pro);

        local = new double[sub_matrix_size * sub_matrix_size];

        if (rank == 0) {
            for (int i = 0; i < sub_matrix_size; i++)
                for (int j = 0; j < sub_matrix_size; j++)
                    local[j + i * sub_matrix_size] = weightTable_mpi[i * n + j];

            // send to others
            double *tmp_m = new double[sub_matrix_size * sub_matrix_size];
            for (int r = 1; r < num_pro; r++) {
                int sub_matrix_X = r / (int) sqrt(num_pro);
                sub_matrix_X *= sub_matrix_size;
                int sub_matrix_Y = r % (int) sqrt(num_pro);
                sub_matrix_Y *= sub_matrix_size;

                // printf("%d %d\n", sub_matrix_X, sub_matrix_Y);

                for (int i = 0; i < sub_matrix_size; i++) {
                    for (int j = 0; j < sub_matrix_size; j++) {
                        tmp_m[j + i * sub_matrix_size] = weightTable_mpi[(i + sub_matrix_X) * n + (j + sub_matrix_Y)];
                    }
                }

                MPI_Send(tmp_m, sub_matrix_size * sub_matrix_size, MPI_DOUBLE, r, 0, MPI_COMM_WORLD);
            }
            free(tmp_m);
        }
        else {
            MPI_Recv(local, sub_matrix_size * sub_matrix_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // build communication
        int sub_matrix_X = rank / (int) sqrt(num_pro);
        int sub_matrix_Y = rank % (int) sqrt(num_pro);
        MPI_Comm_split(MPI_COMM_WORLD, sub_matrix_X, sub_matrix_Y, &comm_row);
        MPI_Comm_split(MPI_COMM_WORLD, sub_matrix_Y, sub_matrix_X, &comm_col);

        rowK = new double[sub_matrix_size];
        colK = new double[sub_matrix_size];

        for (int k = 0; k < n; k++) {

            int k_rank = k / sub_matrix_size;

            if (sub_matrix_X == k_rank) {
                int i = k - sub_matrix_X * sub_matrix_size;
                for (int j = 0; j < sub_matrix_size; j++) {
                    rowK[j] = local[i * sub_matrix_size + j];
                }
            }

            if (sub_matrix_Y == k_rank) {
                int j = k - sub_matrix_Y * sub_matrix_size;
                for (int i = 0; i < sub_matrix_size; i++) {
                    colK[i] = local[i * sub_matrix_size + j];
                }
            }

            MPI_Bcast(rowK, sub_matrix_size, MPI_DOUBLE, k_rank, comm_col);
            MPI_Bcast(colK, sub_matrix_size, MPI_DOUBLE, k_rank, comm_row);

            MPI_Barrier(MPI_COMM_WORLD);

            // build new local
            for (int i = 0; i < sub_matrix_size; i++) {
                for (int j = 0; j < sub_matrix_size; j++) {
                    if (rowK[j] + colK[i] < local[i * sub_matrix_size + j])
                        local[i * sub_matrix_size + j] = rowK[j] + colK[i];
                }
            }

        }

        // for(int r = 0; r < num_pro; r++) {
        //     MPI_Barrier(MPI_COMM_WORLD);
        //     if (r == rank) {
        //         printf("rank: %d\n", rank);
        //         for(int i = 0; i < sub_matrix_size; i++){
        //             printf("row %d: ", i);
        //             for(int j = 0; j < sub_matrix_size; j++)
        //                 if(local[i*sub_matrix_size + j] > 10)
        //                     printf(" 0.00");
        //                 else
        //                     printf(" %1.2f",local[i*sub_matrix_size + j]);
        //             printf("\n");
        //         }
        //     }
        // }

        // gather
        if (rank == 0) {
            out_order_buf = new double[n * n];
            sendcounts = new int[num_pro];
            displs = new int[num_pro];
            for (int i = 0; i < num_pro; ++i) {
                sendcounts[i] = sub_matrix_size * sub_matrix_size;
                displs[i] = i * sendcounts[i];
            }
            MPI_Gatherv(local, sub_matrix_size * sub_matrix_size, MPI_DOUBLE, out_order_buf, sendcounts, displs,
                        MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
        }
        else {
            MPI_Gatherv(local, sub_matrix_size * sub_matrix_size, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0,
                        MPI_COMM_WORLD);
        }

        // reorder buf
        if (rank == 0) {
            int i = 0, j = 0, k = 0, offset = 0;
            for (int x = 0; x < n * n; x++) {
                weightTable_mpi[x] = out_order_buf[k * n * sub_matrix_size + j * sub_matrix_size * sub_matrix_size + i +
                                                   offset];
                if (++i % sub_matrix_size == 0) {
                    i = 0;
                    if (++j % (int) sqrt(num_pro) == 0) {
                        j = 0;
                        offset += sub_matrix_size;
                        if (offset == sub_matrix_size * sub_matrix_size) {
                            offset = 0;
                            k++;
                        }
                    }
                }
            }
            free(out_order_buf);
        }

        free(local);

        clock_gettime(CLOCK_MONOTONIC, &end_ser);

        diff = (double) ((double) BILLION * (end_ser.tv_sec - start_ser.tv_sec) + end_ser.tv_nsec - start_ser.tv_nsec) /
               1000000;

        if (rank == 0) {
            parallelTime += diff;
//            printf("Time taken for parallel   version = %f milliseconds\n", (double) diff);
        }
    }

    if (rank == 0) {
	parallelTime /= (double)ITERS;
        printf("Avg parallel time is %f\n", parallelTime);

        double iso_efficiency = seqTime / (num_pro * parallelTime);
        printf("Iso efficiency is %f\n", iso_efficiency);
    }

    if (rank == 0) {
        bool break_flag = false;
        for (int i = 0; i < n && !break_flag; ++i) {
            for (int j = 0; j < n && !break_flag; ++j) {
                if (weightTable[ind(i, j)] - weightTable_mpi[i * n + j] != 0) {
                    printf("error!\n");
                    // break;
                    break_flag = true;
                }
            }
        }
    }

    MPI_Finalize();
}

