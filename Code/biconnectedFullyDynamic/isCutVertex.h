#ifndef ISCUTVERTEX_H
#define ISCUTVERTEX_H
#include"dynamicGraph.h"
#include"unweightedGraph.h"


bool isCutVertex(dynamicGraph &G, int* unfinishedTraversal, int lca, int* parent,int &totalComponentsGlobal);

void isCutVertexAll(dynamicGraph &G, int* unfinishedTraversal, int* parent, int* isLCA, bool* isCutVertex, int totalVertices, int root, int* isBase,
					int* &localParent, bool* &isSafe, int* &level);

void isCutVertexAffectedOld(dynamicGraph &G, int* unfinishedTraversal, int* parent, int* isLCA, bool* isCutVertex, int totalVertices, int root, int* isBase,
					int* &localParent, bool* &isSafe, int* &level, vector<vector<int>> &children, vector<vector<int>> &affectedLCA,
					int* postProcessed, int batchNumber);

// void isCutVertexAffected(dynamicGraph &G, int* unfinishedTraversal, int* parent, int* isLCA, bool* isCutVertex,
//                         int totalVertices, int root, int* isBase, int* &localParent, bool* &isSafe, int* &level,
//                         vector<vector<int>> &children, int* affectedLCA, int* postProcessed,
//                         int batchNumber, int totalAffectedLCA, int** q, unweightedGraph &OG);

void isCutVertexAffected(unweightedGraph &OG, unweightedGraph &G, SpanningTree &T,
						 Biconnected &B, LCAGraph &L, int* postProcessed, int batchNumber);

#endif
 