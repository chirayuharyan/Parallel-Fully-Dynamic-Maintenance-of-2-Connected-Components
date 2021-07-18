#ifndef DELETEEDGE_H
#define DELETEEDGE_H
#include "unweightedGraph.h"
#include "dynamicGraph.h"
#include<vector>
#include<set>
#include<omp.h>


void deleteNonTreeEdge(int u,int v, unweightedGraph &G, int* parent, int* isLCA, int* unfinishedTraversal,
						vector<set<pair<int,int>> > &setOfTreeEdges, omp_lock_t* lLcaGraph,
						omp_lock_t* lSetOfTreeEdges, vector<int> &result, vector<int> &partiallyAffectedLCA, 
						vector<int> &isVisited, int &findLCACount, dynamicGraph &globalGraph, int* isBase);

void deleteEdgeParallel(unweightedGraph &G, int* parent,int* isLCA, int* unfinishedTraversal, 
						vector<set<pair<int,int>> > &setOfTreeEdges, int count, int totalQueries,
						vector<query> &queries, omp_lock_t* lLcaGraph, omp_lock_t* lSetOfTreeEdges, 
						vector<vector<int>> &affectedLCA, vector<vector<int>> &partiallyAffectedLCA, 
						int &countTreeEdge, vector<vector<int>> &isSpecial,
						vector<vector<int>> &rank, vector<vector<int>> &isVisited, int* findLCACount,
						dynamicGraph &globalGraph, int* isBase);


void deleteNonTreeEdge(int u,int v, unweightedGraph &G, SpanningTree &T, Biconnected &B, LCAGraph &L,int tid);


void deleteEdgeParallel(unweightedGraph &G, int count, int totalQueries, vector<query> &queries,
						  SpanningTree &T, Biconnected &B, LCAGraph &L);

#endif