#include<omp.h>
#include"lca.h"
#include "unweightedGraph.h"
#include <vector>
#include <unordered_map>
#include <iostream>
#include <set>


bool isTreeEdge(int *parent, int u, int v)
{
	if(parent[u] == v || parent[v] == u)
		return true;
	return false; 
}

void updateValues(int u, bool* partOfFundamental, bool* unfinishedTraversal)
{
	unfinishedTraversal[u]+= true;
	partOfFundamental[u]+= true;	
}

void updateValues(int u, bool* partOfFundamental)
{
	partOfFundamental[u]+=true;	
}

int findLCAandUpdate(unweightedGraph &G,int *parent, int *level, int u, int v, bool* partOfFundamental, bool* unfinishedTraversal, 
					dynamicGraph &lcaGraph, omp_lock_t* lLcaGraph, vector<pair<int,int>> &edgeList, bool* isBase)
{
	int currU = u;
	int currV = v;
	if(level[currU] < level[currV])
	{
		updateValues(currV,partOfFundamental,unfinishedTraversal);
		currV = parent[currV];
	}
	if(level[currV] < level[currU])
	{
		updateValues(currU,partOfFundamental,unfinishedTraversal);
		currU = parent[currU];
	}

	while(parent[currU] != parent[currV])
	{
		updateValues(currV,partOfFundamental,unfinishedTraversal);
		currV = parent[currV];

		updateValues(currU,partOfFundamental,unfinishedTraversal);
		currU = parent[currU];
	}

	updateValues(currU, partOfFundamental);
	updateValues(currV, partOfFundamental);
	
	/*omp_set_lock(&(lLcaGraph[currU]));
	lcaGraph.addEdge(currU,currV);
	omp_unset_lock(&(lLcaGraph[currU]));	

	omp_set_lock(&(lLcaGraph[currV]));
	lcaGraph.addEdge(currV,currU);
	omp_unset_lock(&(lLcaGraph[currV]));*/	

	isBase[currU] += true;
	isBase[currV] += true;
	// edgeList.push_back(make_pair(currU, currV));

	return parent[currU];
}


