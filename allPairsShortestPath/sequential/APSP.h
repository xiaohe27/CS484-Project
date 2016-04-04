//
// Created by xiaohe on 4/3/16.
//
#include "graph.h"
#include "graphIO.h"
#include "parallel.h"

#ifndef CS484_PROJECT_APSP_H
#define CS484_PROJECT_APSP_H

using namespace std;

extern "C" {
 char* allPairsShortestPath(graph<intT> Gr);
}
#endif //CS484_PROJECT_APSP_H
