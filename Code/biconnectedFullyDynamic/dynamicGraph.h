#ifndef DYNAMICGRAPH_H
#define DYNAMICGRAPH_H

#include<vector>
#include<map>
using namespace std;
 
typedef struct unweightedEdge
{
	int u,v,cost;
}unweightedEdge;

class dynamicGraph
{
public:
	int totalVertices, totalEdges;
	map<int,map<int,unweightedEdge> > adjacencyList;
	vector<map<int,int> > adjacencyList1;

public:
	dynamicGraph();
	void print();
	bool addEdge(int,int);
	void deleteEdge(int,int);
	bool addEdge1(int,int);
	bool deleteEdge1(int,int);
	void addVertex(int);
	void deleteVertex(int);
};

#endif