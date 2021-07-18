#ifndef ADDEDGE_H
#define ADDEDGE_H
#include "unweightedGraph.h"
#include <omp.h>
#include "dynamicGraph.h"
#include<vector>
#include<set>

void copyFromSetToVector(set<pair<int,int>> s, vector<pair<int,int>> &v);

void addNonTreeEdge(int u,int v, unweightedGraph &G, int* parent, int* isLCA, int* unfinishedTraversal, 
					vector<set<pair<int,int>> > &setOfTreeEdges, omp_lock_t* lLcaGraph,
					omp_lock_t* lSetOfTreeEdges, vector<int> &result, vector<int> &partiallyAffectedLCA, 
					vector<int> &isVisited, int &findLCACount, dynamicGraph &globalGraph, int* isBase);

void addNonTreeEdge(int u,int v, unweightedGraph &G, SpanningTree &T, Biconnected &B, LCAGraph &L,int tid);

void addTreeEdge(int u,int v, unweightedGraph &G, int* parent,int* isLCA, int* unfinishedTraversal,
					vector<set<pair<int,int>> > &setOfTreeEdges, omp_lock_t* lLcaGraph,
					omp_lock_t* lSetOfTreeEdges, vector<int> &result, int count, vector<int> &isSpecial,
					vector<int> &rank, int tid, dynamicGraph &globalGraph, int* isBase);


void addTreeEdge(int u,int v, int count,unweightedGraph &G, SpanningTree &T, Biconnected &B, LCAGraph &L,int tid);

#endif 