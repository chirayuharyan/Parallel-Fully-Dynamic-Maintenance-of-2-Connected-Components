#include"dynamicGraph.h"
#include<iostream>
#include<unordered_map>
#include<vector>
using namespace std; 

dynamicGraph::dynamicGraph()
{
	totalVertices = 0;
	totalEdges = 0;
}

void dynamicGraph::addVertex(int u)
{
	if(adjacencyList.find(u)==adjacencyList.end())
	{
		map<int,unweightedEdge> temp;
		adjacencyList[u] = temp;
		totalVertices++;
	}
}

void dynamicGraph::deleteVertex(int u)
{
	map<int,map<int,unweightedEdge>>::iterator it = adjacencyList.find(u);
	if(it != adjacencyList.end())
	{
		adjacencyList.erase(it);
		totalVertices--;
	}
}

bool dynamicGraph::addEdge(int u, int v) //adds u,v and v,u
{
	addVertex(u);
	addVertex(v);

	bool flag = false;
	if(adjacencyList[u].find(v) == adjacencyList[u].end())
	{
		unweightedEdge temp;
		temp.u = u;
		temp.v = v;
		temp.cost = 0;
		adjacencyList[u][v] = temp;
		flag = true;
		totalEdges++;
	}
	adjacencyList[u][v].cost++;


	if(adjacencyList[v].find(u) == adjacencyList[v].end())
	{
		unweightedEdge temp;
		temp.u = v;
		temp.v = u;
		temp.cost = 0;
		adjacencyList[v][u] = temp;
	}
	adjacencyList[v][u].cost++;

	return flag;
}

void dynamicGraph::deleteEdge(int u, int v)
{
	adjacencyList[u][v].cost--;
	adjacencyList[v][u].cost--;
	if(adjacencyList[u][v].cost == 0)
	{
		map<int,unweightedEdge>::iterator it;
		it = adjacencyList[u].find(v);
		adjacencyList[u].erase(it);
		totalEdges--;
	}

	if(adjacencyList[v][u].cost == 0)
	{
		map<int,unweightedEdge>::iterator it;
		it = adjacencyList[v].find(u);
		adjacencyList[v].erase(it);
	}

	if(adjacencyList[u].size() == 0)
	{
		deleteVertex(u);
	}
	if(adjacencyList[v].size() == 0)
	{
		deleteVertex(v);
	}
}

void printEdge(unweightedEdge e)
{
	cout<<e.u<<" "<<e.v<<" "<<e.cost<<endl;
}

void dynamicGraph::print()
{
	cout<<totalVertices<<" "<<totalEdges<<endl;
	map<int,map<int,unweightedEdge>>::iterator it;
	map<int,unweightedEdge>::iterator it2;
	for(it = adjacencyList.begin();it!=adjacencyList.end();it++)
	{
		for(it2 = it->second.begin(); it2 != it->second.end(); it2++)
		{
			printEdge(it2->second);
		}
	}
}

bool dynamicGraph::addEdge1(int u, int v)
{
	auto it = adjacencyList1[u].find(v);
	if(it == adjacencyList1[u].end())
	{
		adjacencyList1[u][v] = 1;
		return true;
	} 
	else
	{
		adjacencyList1[u][v]++;
		return false;
	}

	// adjacencyList1[u][v]++;
}


bool dynamicGraph::deleteEdge1(int u, int v)
{
	auto it = adjacencyList1[u].find(v);
	if(it == adjacencyList1[u].end())
	{
		// printf("Problem deleting base edge %d,%d\n",u,v);
		return false; 	
	} 
	if(it->second == 1)
	{
		adjacencyList1[u].erase(it);
		return true;
	}
	else
	{
		adjacencyList1[u][v]--;
		return false;
	}
}
 