//
// Created by xiaohe on 4/3/16.
//
#include "graph.h"
#include "graphIO.h"
#include "parallel.h"

#ifndef CS484_PROJECT_APSP_H
#define CS484_PROJECT_APSP_H

//#define ind(x,y)

using namespace std;

extern "C" {
 char* allPairsShortestPath(graph<intT> Gr);
}

inline void printMatrix(double **matrix, intT N) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << matrix[i][j] << ", ";
        }

        cout << "\n" << endl;
    }
}
#endif //CS484_PROJECT_APSP_H
