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

    double *local, *row_m, *col_m, *out_order_buf, *res_m;

    int sub_matrix_size, row_per_pro;

    double diff;

    double *weightTable;
    double *weightTable_mpi;

    MPI_Status status;
    MPI_Comm comm_row;
    MPI_Comm comm_col;

    MPI_Init(&argc, &argv);

    int rank, num_pro;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_pro);

    if(rank == 0){
    
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
        weightTable = new double [n*n];
        weightTable_mpi = new double [n*n];
        for (int i = 0; i < n; ++i) {
            // weightTable[i] = new double[n];
            // weightTable_mpi[i] = new double[n];
            for (int j = 0; j < n; ++j){
                weightTable[ind(i,j)] = numeric_limits<double>::max();
                weightTable_mpi[ind(i,j)] = numeric_limits<double>::max();
            }
        }
        for (int i = 0; i < Gr.m; ++i) {
            curEdge = edgeList[i];
            int u = curEdge.u;
            int v = curEdge.v;
            double weight = curEdge.weight;

            // cout << "edge " << i << " connects node " << u << " and node " << v
            // << ", with weight " << weight << endl;

            if (weightTable[ind(u,v)] > weight){
                weightTable[ind(u,v)] = weight;
                weightTable_mpi[ind(u,v)] = weight;
            }
        }
    }

    struct timespec start_ser, end_ser;


    // sequential version
    clock_gettime(CLOCK_MONOTONIC,&start_ser);

    if(rank == 0){
        for (int k = 0; k < n; ++k) {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < n; ++i) {
                    double viaK = weightTable[ind(i,k)] + weightTable[ind(k,j)];
                    if (weightTable[ind(i,j)] > viaK)
                        weightTable[ind(i,j)] = viaK;
                }
            }
        }
    }


    clock_gettime(CLOCK_MONOTONIC,&end_ser);

    diff = (double)((double)BILLION*(end_ser.tv_sec-start_ser.tv_sec)+end_ser.tv_nsec-start_ser.tv_nsec)/1000000;

    if(rank==0)
    printf("Time taken for sequential version = %f milliseconds\n",(double)diff);


    // parallel version
    clock_gettime(CLOCK_MONOTONIC,&start_ser);

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // printf("here rank %d has n %d\n", rank, n);

    int sqrt_q = (int)sqrt(num_pro);

    sub_matrix_size = n / sqrt_q;

    local = new double[sub_matrix_size * sub_matrix_size];

    if(rank==0){
        for(int i = 0; i < sub_matrix_size; i++)
            for(int j = 0; j < sub_matrix_size; j++)
                local[j + i * sub_matrix_size] = weightTable_mpi[i*n + j];

        // send to others
        double *tmp_m = new double[sub_matrix_size * sub_matrix_size];
        for(int r = 1; r < num_pro; r++){
            int sub_matrix_X = r / sqrt_q; sub_matrix_X *= sub_matrix_size;
            int sub_matrix_Y = r % sqrt_q; sub_matrix_Y *= sub_matrix_size;

            // printf("%d %d\n", sub_matrix_X, sub_matrix_Y);

            for(int i = 0; i < sub_matrix_size; i++){
                for(int j = 0; j < sub_matrix_size; j++){
                    tmp_m[j + i * sub_matrix_size] = weightTable_mpi[(i + sub_matrix_X) * n + (j + sub_matrix_Y)];
                }
            }

            MPI_Send(tmp_m, sub_matrix_size * sub_matrix_size, MPI_DOUBLE, r, 0, MPI_COMM_WORLD);
        }
        free(tmp_m);
    }
    else{
        MPI_Recv(local, sub_matrix_size * sub_matrix_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // build communication
    int sub_matrix_X = rank / sqrt_q;
    int sub_matrix_Y = rank % sqrt_q;


    // http://stackoverflow.com/questions/11372012/mpi-several-broadcast-at-the-same-time
    int dims[2] = {sqrt_q, sqrt_q};
    int periods[2] = {1, 1};
    MPI_Comm comm_cart;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_cart);

    int row_dims[2] = {0, 1};
    MPI_Cart_sub(comm_cart, row_dims, &comm_row);

    int col_dims[2] = {1, 0};
    MPI_Cart_sub(comm_cart, col_dims, &comm_col);

    // MPI_Comm_split(MPI_COMM_WORLD, sub_matrix_X, sub_matrix_Y, &comm_row);
    // MPI_Comm_split(MPI_COMM_WORLD, sub_matrix_Y, sub_matrix_X, &comm_col);

    row_m = new double[sub_matrix_size*sub_matrix_size];
    col_m = new double[sub_matrix_size*sub_matrix_size];
    res_m = new double[sub_matrix_size*sub_matrix_size];

    memcpy(res_m, local, sub_matrix_size * sub_matrix_size * sizeof(double));

    int src_col = (sub_matrix_X + 1) % sqrt_q;
    int dst_col = (sub_matrix_X + sqrt_q - 1) % sqrt_q;


    for(int d = 1; d < n; d <<= 1){
    // for(int d = 0; d < 3; d++){
        memcpy(col_m, local, sub_matrix_size * sub_matrix_size * sizeof(double));
        
        for(int step = 0; step < sqrt_q; step++){
            int bcast = (sub_matrix_X + step) % sqrt_q;
            if(sub_matrix_Y == bcast){
                memcpy(row_m, local, sub_matrix_size * sub_matrix_size * sizeof(double));
                MPI_Bcast(local, sub_matrix_size * sub_matrix_size, MPI_DOUBLE, bcast, comm_row);
            }
            else{
                // int sour = sub_matrix_X * sqrt_q + bcast;
                MPI_Bcast(row_m, sub_matrix_size * sub_matrix_size, MPI_DOUBLE, bcast, comm_row);
            }

            for (int k = 0; k < sub_matrix_size; ++k) {
                for (int j = 0; j < sub_matrix_size; ++j) {
                    for (int i = 0; i < sub_matrix_size; ++i) {
                        double viaK = row_m[i*sub_matrix_size + k] + col_m[k * sub_matrix_size + j];
                        if (res_m[i*sub_matrix_size + j] > viaK)
                            res_m[i*sub_matrix_size + j] = viaK;
                    }
                }
            }


            MPI_Sendrecv_replace(col_m, sub_matrix_size * sub_matrix_size, MPI_DOUBLE, dst_col, 0, src_col, 0, comm_col, MPI_STATUS_IGNORE);
        
        }

        MPI_Barrier(MPI_COMM_WORLD);
        memcpy(local, res_m, sub_matrix_size * sub_matrix_size * sizeof(double));
    }

    // gather
    if(rank == 0){
        out_order_buf = new double[n * n];
        sendcounts = new int[num_pro];
        displs = new int[num_pro];
        for (int i = 0; i <num_pro; ++i) {
            sendcounts[i] = sub_matrix_size * sub_matrix_size;
            displs[i] = i * sendcounts[i];
        }
        MPI_Gatherv(local, sub_matrix_size * sub_matrix_size, MPI_DOUBLE, out_order_buf, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    else{
        MPI_Gatherv(local, sub_matrix_size * sub_matrix_size, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // reorder buf
    if(rank == 0){
        int i = 0, j = 0, k = 0, offset = 0;
        for(int x = 0; x < n*n; x++){
            weightTable_mpi[x] = out_order_buf[k*n*sub_matrix_size + j*sub_matrix_size*sub_matrix_size + i + offset];
            if(++i % sub_matrix_size == 0){
                i = 0;
                if(++j % sqrt_q == 0){
                    j = 0;
                    offset += sub_matrix_size;
                    if(offset == sub_matrix_size * sub_matrix_size){
                        offset = 0;
                        k++;
                    }
                }
            }
        }
        free(out_order_buf);
    }

    free(local);

    clock_gettime(CLOCK_MONOTONIC,&end_ser);

    diff = (double)((double)BILLION*(end_ser.tv_sec-start_ser.tv_sec)+end_ser.tv_nsec-start_ser.tv_nsec)/1000000;

    if(rank==0)
    printf("Time taken for parallel   version = %f milliseconds\n",(double)diff);

    if(rank==0){
        bool break_flag = false;
        for (int i = 0; i < n && !break_flag; ++i) {
            for (int j = 0; j < n && !break_flag; ++j) {
                if(abs(weightTable[ind(i,j)] - weightTable_mpi[i*n + j]) > 0.00001){
                    printf("error at %d, %d \n", i, j);
                    // break;
                    break_flag = true;
                }
            }
        }
    }

    MPI_Finalize();
}

