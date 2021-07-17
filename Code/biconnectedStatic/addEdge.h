#ifndef ADDEDGE_H
#define ADDEDGE_H
#include "unweightedGraph.h"
#include <omp.h>
#include "dynamicGraph.h"
#include<vector>
#include<set>


void addNonTreeEdge(int u,int v, unweightedGraph &G, int* parent, int* level, dynamicGraph &lcaGraph, bool* isLCA, bool* isBase,
					bool* partOfFundamental, bool* unfinishedTraversal, omp_lock_t* lLcaGraph, vector<pair<int,int>> &edgeList);


#endif 