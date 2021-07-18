#ifndef UNWEIGHTEDGRAPH_H
#define UNWEIGHTEDGRAPH_H

#include<fstream>
#include<vector>
#include<set>
#include<omp.h>
#include"dynamicGraph.h"
using namespace std;

typedef struct SpanningTree
{
	int *parent,*levels;
}SpanningTree;

typedef struct Biconnected
{
	int *unfinishedTraversal,*isBase,*isLCA;
	vector<set<pair<int,int>> > setOfTreeEdges;
	omp_lock_t* lSetOfTreeEdges;
	bool* cutVertex;
	bool* cutVertexAsBridge;
	int** isSpecial;
	int countTreeEdge;
	int** rank;
	int** findLCAisVisited;
	int* findLCACount;
}Biconnected;

typedef struct LCAGraph
{
	int *affectedLCA,headAffectedLCA;
	int *partiallyAffectedLCA,headPartiallyAffectedLCA;

	int** partiallyAffectedBaseVertices;
	int* headPartiallyAffectedBaseVertices;

	int** affectedBaseVertices;
	int* headAffectedBaseVertices;

	dynamicGraph globalGraph;
	omp_lock_t* lLcaGraph;
	int* localParent;
	bool* isSafe;
	int* level;
	int** myQueue;
}LCAGraph;


typedef struct edge
{
	int u,v; 
}edge;

typedef struct query
{
	int u,v;
}query;

class unweightedGraph
{
public:
	int totalVertices, totalEdges, *offset, *neighbour,root, maxDegree;
	vector<vector<edge> > adjacencyList;
	int *degree;

public:
	// Graph();
	unweightedGraph(ifstream &edgeList);
	void addEdge(int u, int v);
	void deleteEdge(int u, int v);
	void print();
	void printCSR();
	void addEdgeDegree(int u, int v);
	void deleteEdgeDegree(int u, int v);
};

void printEdge(edge e);
#define out_degree(G, n) (G.offset[n+1] - G.offset[n])
#define out_vertices(G, n) (&G.neighbour[G.offset[n]])

#endif