#ifndef LCA_H
#define LCA_H
#include "unweightedGraph.h"
#include <omp.h>
#include "dynamicGraph.h"
#include <vector>
#include <unordered_map>
#include <set>

/*void LCA(unweightedGraph &G, int *parent, vector<dynamicGraph> &lcaGraph, vector<int> &isLCA, int* isBridge, int* unfinishedTraversal, vector<set<pair<int,int>> > &setOfTreeEdges,
			omp_lock_t* lParent, omp_lock_t* lLcaGraph, omp_lock_t* lIsLCA, omp_lock_t* lIsBridge, omp_lock_t* lUnfinishedTraversal, omp_lock_t* lSetOfTreeEdges);*/

bool isTreeEdge(int *parent, int u, int v);

int findLCA(unweightedGraph &G,int *parent, int u, int v,vector<int> &isVisited, int &count);

int findLCA(int u, int v,unweightedGraph &G, SpanningTree &T, Biconnected &B, int tid);

void assignXY(int &x, int &y, int u, int v);

int updateTreeEdgeValues(unweightedGraph &G, int* parent, int lca, int u,int v, int* unfinishedTraversal,
						 vector<set<pair<int,int>> > &setOfTreeEdges, int* isLCA, vector<int> &result,
						 omp_lock_t* lSetOfTreeEdges);

int updateTreeEdgeValues(int lca, int u, int v, unweightedGraph &G, SpanningTree &T, Biconnected &B, LCAGraph &L);


int updateTreeEdgeValuesDelete(unweightedGraph &G, int* parent, int lca, int u,int v, int* unfinishedTraversal,
						 		vector<set<pair<int,int>> > &setOfTreeEdges, int* isLCA, vector<int> &result,
								omp_lock_t* lSetOfTreeEdges);

int updateTreeEdgeValuesDelete(int lca, int u, int v, unweightedGraph &G, SpanningTree &T, Biconnected &B, 
								LCAGraph &L);

void constructLCAGraph(int lca, unweightedGraph &G, int* parent, int u, int v, int* unfinishedTraversal,
						vector<set<pair<int,int>> > &setOfTreeEdges, int* isLCA, vector<int> &partiallyAffectedLCA,
						omp_lock_t* lLcaGraph, omp_lock_t* lSetOfTreeEdges, dynamicGraph &globalGraph,
						int* isBase);

bool constructLCAGraph(int lca, int u, int v, unweightedGraph &G, SpanningTree &T, Biconnected &B, LCAGraph &L);

void destructLCAGraph(int lca, unweightedGraph &G, int* parent, int u, int v, int* unfinishedTraversal,
					 	vector<set<pair<int,int>> > &setOfTreeEdges, int* isLCA, vector<int> &partiallyAffectedLCA,
					 	omp_lock_t* lLcaGraph, omp_lock_t* lSetOfTreeEdges, dynamicGraph &globalGraph,
						int* isBase);

bool destructLCAGraph(int lca, int u, int v, unweightedGraph &G, SpanningTree &T, Biconnected &B, LCAGraph &L);


#endif