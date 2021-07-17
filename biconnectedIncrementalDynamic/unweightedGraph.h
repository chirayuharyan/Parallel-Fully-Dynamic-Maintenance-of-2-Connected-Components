#ifndef UNWEIGHTEDGRAPH_H
#define UNWEIGHTEDGRAPH_H

#include<fstream>
#include<vector>
#include<set>
using namespace std;

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
	int totalVertices, totalEdges, *offset, *neighbour, *U, *V, root,*degree;
	vector<vector<edge> > adjacencyList;
	set<pair<int,int>> edgeSet;

public:
	// Graph();
	unweightedGraph();
	unweightedGraph(ifstream &edgeList);
	void addEdge(int u, int v);
	void deleteEdge(int u, int v);
	void print();
	void printCSR();
	void addEdgeSet(int u, int v);
	void deleteEdgeSet(int u, int v);
	void recreate();
};

void printEdge(edge e);

#define out_degree(G, n) (G.offset[n+1] - G.offset[n])
#define out_vertices(G, n) (&G.neighbour[G.offset[n]])

#endif