#include"addEdge.h"
#include"lca.h"
#include<iostream>
#include<vector>
#include<set>
#include <sstream>
#include<omp.h>


void addNonTreeEdge(int u,int v, unweightedGraph &G, int* parent, int* level, dynamicGraph &lcaGraph, bool* isLCA, bool* isBase,
					bool* partOfFundamental, bool* unfinishedTraversal, omp_lock_t* lLcaGraph, vector<pair<int,int>> &edgeList)
{
	int LCAVertex = findLCAandUpdate(G, parent, level, u, v, partOfFundamental, unfinishedTraversal, lcaGraph, lLcaGraph, edgeList, isBase);
	isLCA[LCAVertex] += true;
}

