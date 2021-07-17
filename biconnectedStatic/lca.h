#ifndef LCA_H
#define LCA_H
#include "unweightedGraph.h"
#include <omp.h>
#include <vector>
#include <unordered_map>
#include <set>
#include <utility>

bool isTreeEdge(int *parent, int u, int v);

int findLCAandUpdate(unweightedGraph &G,int *parent, int *level, int u, int v, bool* partOfFundamental, bool* unfinishedTraversal, 
					dynamicGraph &lcaGraph, omp_lock_t* lLcaGraph, vector<pair<int,int>> &edgeList, bool* isBase);

#endif