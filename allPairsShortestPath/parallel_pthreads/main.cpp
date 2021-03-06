#include "main.h"
#include "gettime.h"
#include <pthread.h>
#include "parseCommandLine.h"

using namespace std;

int num_threads = 16;

struct data {
 double ** weightTable;
 int size;
 int index;
 int j;
 int k;
};

void * computeWeightTable(void * data)
{
    struct data * paralleldata = (struct data *)data;

    double ** weightTable = paralleldata->weightTable;
    int j = paralleldata->j;
    int k = paralleldata->k;
    int index = paralleldata->index;
    int size = paralleldata->size;

    //printf("j = %d, k = %d\n", j, k);

    for (int i = 0; i < size; ++i) {
        double viaK = weightTable[i + index][k] + weightTable[k][j];
        if (weightTable[i + index][j] > viaK)
            weightTable[i + index][j] = viaK;
    }

    free(data);
    return NULL;
}

void allPairsShortestPath(wghEdgeArray<intT> Gr) {
    intT n = Gr.n;
    wghEdge<intT> *edgeList = Gr.E;

    cout << "number of nodes is " << n << ", and number of edges is " << Gr.m << "." << endl;

    wghEdge<intT> curEdge;

    //the initial shortest path table is a n*n matrix
    //so we need to analyze the original weighted multi-graph
    //to obtain the initial configuration.
    double **weightTable = new double *[n];
    double **weightTable_pthreads = new double *[n];

    printf("Number of threads = %d\n", num_threads);

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

//        cout << "edge " << i << " connects node " << u << " and node " << v
//        << ", with weight " << weight << endl;

        if (weightTable[u][v] > weight)
            weightTable[u][v] = weight;
    }
    for (int i = 0; i < n; ++i)
    {
        weightTable_pthreads[i] = new double[n];
        for (int j = 0 ; j < n; j++)
        {

            weightTable_pthreads[i][j] = weightTable[i][j];
        }
    }
   
    /*

    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                double viaK = weightTable[i][k] + weightTable[k][j];
                if (weightTable[i][j] > viaK)
                    weightTable[i][j] = viaK;
            }
        }
    }

    */

    pthread_t * threads = (pthread_t *)malloc(sizeof(pthread_t) * num_threads);
    int remain = 0;
    int split = 0;
    int count = 0;

    split = (int)(n / num_threads);

    _tm.start();

    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            count = 0;
            remain = n % num_threads;

            for (int i = 0; i < num_threads; ++i) {

                struct data * data = (struct data *)malloc(sizeof( struct data));
                data->weightTable = weightTable_pthreads;
                data->size = split;
                data->j = j;
                data->k = k;

                if(remain-- > 0)
                    data->size++;

                data->index = count;
                count += data->size;
                pthread_create(&(threads[i]), NULL, computeWeightTable, data);
            }

            for(int i = 0; i < num_threads; i++)
                pthread_join(threads[i], NULL);

        }
    }
    printf("Time = %f\n", _tm.stop());
    free(threads);

    /*

    bool break_flag = false;
    for (int i = 0; i < n && !break_flag; ++i) {
        for (int j = 0; j < n && !break_flag; ++j) {
            if(weightTable[i][j] - weightTable_pthreads[i][j] != 0){
                printf("%f %f\n", weightTable[i][j] , weightTable_pthreads[i][j]);
                break_flag = true;
            }
        }
    }

    */

    for (int k = 0; k < n; ++k) {
        delete[] weightTable[k];
    }

    delete[] weightTable;
}

int main(int argc, char **argv) {
    commandLine P(argc, argv, "-n <numThreads> -o <outFile>");
    char *iFile = P.getArgument(0);
    char *oFile = P.getOptionValue("-o");

    num_threads = P.getOptionIntValue("-n", 16);

    wghEdgeArray<intT> G = benchIO::readWghEdgeArrayFromFile<intT>(iFile);

    allPairsShortestPath(G);
}

