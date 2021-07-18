#include<omp.h>
#include"lca.h"
#include "unweightedGraph.h"
#include "dynamicGraph.h"
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

int getRoot(int i, int* parent)
{
	int root = i;
	while(parent[root]!=root)
	{
		root = parent[root];
	}
	return root;
}

int findLCA(unweightedGraph &G,int *parent, int u, int v, vector<int> &isVisited, int &count)
{
	int countU = count;
	int countV = count+1;
	isVisited[u] = countU;
	isVisited[v] = countV;
	int currU = u;
	int currV = v;
	int lca = -1;
	bool turn = 0;
	
	while(1)
	{
		if(turn==0) //turn 0 means u will mark
		{
			currU = parent[currU];

			if(isVisited[currU] == countV)
			{
				lca = currU;
				break;
			}
			else
			{
				isVisited[currU] = countU;
				turn = 1;
			}
		}
		else
		{
			currV = parent[currV];
			if(isVisited[currV] == countU)
			{
				lca = currV;
				break;
			}
			else
			{
				isVisited[currV] = countV;
				turn = 0;
			}
		}
	} 
	
	count+=2;
	return lca;
}

int findLCA(int u, int v, unweightedGraph &G, SpanningTree &T, Biconnected &B, int tid)
{
	int count = B.findLCACount[tid];
	int countU = count;
	int countV = count+1;
	int *isVisited = B.findLCAisVisited[tid];
	isVisited[u] = countU;
	isVisited[v] = countV;
	int currU = u;
	int currV = v;
	int lca = -1;
	bool turn = 0;
	
	while(1)
	{
		if(turn==0) //turn 0 means u will mark
		{
			currU = T.parent[currU];

			if(isVisited[currU] == countV)
			{
				lca = currU;
				break;
			}
			else
			{
				isVisited[currU] = countU;
				turn = 1;
			}
		}
		else
		{
			currV = T.parent[currV];
			if(isVisited[currV] == countU)
			{
				lca = currV;
				break;
			}
			else
			{
				isVisited[currV] = countV;
				turn = 0;
			}
		}
	} 
	
	B.findLCACount[tid]+=2;
	return lca;
}

void assignXY(int &x, int &y, int u, int v)
{
	x = u;
	y = v;
	if(x>y)
	{
		x = v;
		y = u;
	}
}

int updateTreeEdgeValues(unweightedGraph &G, int* parent, int lca, int u,int v, int* unfinishedTraversal,
						 vector<set<pair<int,int>> > &setOfTreeEdges, int* isLCA, vector<int> &result,
						 omp_lock_t* lSetOfTreeEdges)
{
	int currU = u;
	int x,y;
	int u1, v1;
	int oldValue;
	assignXY(u1, v1, u ,v);
	
	pair<int,int> currNonTreeEdge = make_pair(u1,v1); 

	while(parent[currU] != lca)
	{
		assignXY(x,y,currU,parent[currU]);

		oldValue = __sync_fetch_and_add(&unfinishedTraversal[currU],1);

		omp_set_lock(&(lSetOfTreeEdges[currU]));
		setOfTreeEdges[currU].insert(currNonTreeEdge);
		omp_unset_lock(&(lSetOfTreeEdges[currU]));


		if(oldValue == 0 && isLCA[parent[currU]] > 0)
		{
			result.push_back(parent[currU]);
			// oldValue = __sync_fetch_and_add(&unfinishedTraversal[currU],1);
		
		}

		currU = parent[currU];
	}

	assignXY(x,y,currU,parent[currU]);
	
	omp_set_lock(&(lSetOfTreeEdges[currU]));
	setOfTreeEdges[currU].insert(currNonTreeEdge);
	omp_unset_lock(&(lSetOfTreeEdges[currU]));
	

	return currU;
}

int updateTreeEdgeValues(int lca, int u, int v, unweightedGraph &G, SpanningTree &T, Biconnected &B, LCAGraph &L)
{
	int currU = u;
	int x,y;
	int u1, v1;
	int oldValue;
	assignXY(u1, v1, u ,v);
	
	pair<int,int> currNonTreeEdge = make_pair(u1,v1); 

	while(T.parent[currU] != lca)
	{
		// assignXY(x,y,currU,T.parent[currU]);

		oldValue = __sync_fetch_and_add(&B.unfinishedTraversal[currU],1);

		omp_set_lock(&(B.lSetOfTreeEdges[currU]));
		B.setOfTreeEdges[currU].insert(currNonTreeEdge);
		omp_unset_lock(&(B.lSetOfTreeEdges[currU]));


		if(oldValue == 0 && B.isLCA[T.parent[currU]] > 0)
		{
			oldValue = __sync_fetch_and_add(&L.headPartiallyAffectedLCA,1);
			L.partiallyAffectedLCA[oldValue] = T.parent[currU];
			// oldValue = __sync_fetch_and_add(&L.headPartiallyAffectedBaseVertices[T.parent[currU]],1);
			// L.partiallyAffectedBaseVertices[T.parent[currU]][oldValue] = currU;
			// result.push_back(parent[currU]);
		}

		currU = T.parent[currU];
	}

	// assignXY(x,y,currU,T.parent[currU]);
	
	omp_set_lock(&(B.lSetOfTreeEdges[currU]));
	B.setOfTreeEdges[currU].insert(currNonTreeEdge);
	omp_unset_lock(&(B.lSetOfTreeEdges[currU]));
	

	return currU;
}


int updateTreeEdgeValuesDelete(unweightedGraph &G, int* parent, int lca, int u,int v, int* unfinishedTraversal,
						 		vector<set<pair<int,int>> > &setOfTreeEdges, int* isLCA, vector<int> &result,
								omp_lock_t* lSetOfTreeEdges)
{
	int currU = u;
	int x,y;
	set<pair<int,int>>::iterator it;

	int u1, v1;
	assignXY(u1, v1, u ,v);
	int newValue;
	
	pair<int,int> currNonTreeEdge = make_pair(u1,v1); 
	
	while(parent[currU] != lca)
	{
		assignXY(x,y,currU,parent[currU]);

		newValue = __sync_add_and_fetch(&unfinishedTraversal[currU],-1);

		omp_set_lock(&(lSetOfTreeEdges[currU]));
		
		it = setOfTreeEdges[currU].find(currNonTreeEdge);
		if(it!= setOfTreeEdges[currU].end())
			setOfTreeEdges[currU].erase(it);
		omp_unset_lock(&(lSetOfTreeEdges[currU]));
		
		if(newValue == 0 && isLCA[parent[currU]] > 0)
			result.push_back(parent[currU]);

		currU = parent[currU];
	}
	assignXY(x,y,currU,parent[currU]);
		
	omp_set_lock(&(lSetOfTreeEdges[currU]));	
	
	it = setOfTreeEdges[currU].find(currNonTreeEdge);
	if(it!= setOfTreeEdges[currU].end())
		setOfTreeEdges[currU].erase(it);
	
	omp_unset_lock(&(lSetOfTreeEdges[currU]));

	return currU;
}


int updateTreeEdgeValuesDelete(int lca, int u, int v, unweightedGraph &G, SpanningTree &T, Biconnected &B, 
								LCAGraph &L)
{
	int currU = u;
	int x,y;
	set<pair<int,int>>::iterator it;

	int u1, v1;
	assignXY(u1, v1, u ,v);
	int newValue;
	
	pair<int,int> currNonTreeEdge = make_pair(u1,v1); 
	
	while(T.parent[currU] != lca)
	{
		// assignXY(x,y,currU,parent[currU]);

		newValue = __sync_add_and_fetch(&B.unfinishedTraversal[currU],-1);

		omp_set_lock(&(B.lSetOfTreeEdges[currU]));
		
		it = B.setOfTreeEdges[currU].find(currNonTreeEdge);
		if(it!= B.setOfTreeEdges[currU].end())
			B.setOfTreeEdges[currU].erase(it);
		omp_unset_lock(&(B.lSetOfTreeEdges[currU]));
		
		if(newValue == 0 && B.isLCA[T.parent[currU]] > 0)
		{
			newValue = __sync_fetch_and_add(&L.headPartiallyAffectedLCA,1);
			L.partiallyAffectedLCA[newValue] = T.parent[currU];
			// newValue = __sync_fetch_and_add(&L.headPartiallyAffectedBaseVertices[T.parent[currU]],1);
			// L.partiallyAffectedBaseVertices[T.parent[currU]][newValue] = currU;
		}

		currU = T.parent[currU];
	}
	// assignXY(x,y,currU,parent[currU]);
		
	omp_set_lock(&(B.lSetOfTreeEdges[currU]));	
	
	it = B.setOfTreeEdges[currU].find(currNonTreeEdge);
	if(it!= B.setOfTreeEdges[currU].end())
		B.setOfTreeEdges[currU].erase(it);
	
	omp_unset_lock(&(B.lSetOfTreeEdges[currU]));

	return currU;
}



void constructLCAGraph(int lca, unweightedGraph &G, int* parent, int u, int v, int* unfinishedTraversal,
						vector<set<pair<int,int>> > &setOfTreeEdges, int* isLCA, vector<int> &partiallyAffectedLCA,
						omp_lock_t* lLcaGraph, omp_lock_t* lSetOfTreeEdges, dynamicGraph &globalGraph, int* isBase)
{
	int baseVertexU,baseVertexV;
	if(u==lca)
		baseVertexU = v;
	else
		baseVertexU = updateTreeEdgeValues(G, parent, lca, u, v, unfinishedTraversal,setOfTreeEdges,isLCA,
											partiallyAffectedLCA, lSetOfTreeEdges);
	
	if(v==lca)
		baseVertexV = u;
	else
		baseVertexV = updateTreeEdgeValues(G,parent,lca,v,u ,unfinishedTraversal,setOfTreeEdges,isLCA,
											partiallyAffectedLCA, lSetOfTreeEdges);
	

	
	if(isTreeEdge(parent,baseVertexU,lca))
		__sync_fetch_and_add(&isBase[baseVertexU],1);

	if(isTreeEdge(parent,baseVertexV,lca))
		__sync_fetch_and_add(&isBase[baseVertexV],1);

	if(!isTreeEdge(parent,baseVertexU,lca) || !isTreeEdge(parent,baseVertexV,lca))
		return;

	
	omp_set_lock(&(lLcaGraph[baseVertexU]));
	globalGraph.addEdge1(baseVertexU,baseVertexV);
	omp_unset_lock(&(lLcaGraph[baseVertexU]));

	omp_set_lock(&(lLcaGraph[baseVertexV]));
	globalGraph.addEdge1(baseVertexV,baseVertexU);
	omp_unset_lock(&(lLcaGraph[baseVertexV]));
	
}

bool constructLCAGraph(int lca, int u, int v, unweightedGraph &G, SpanningTree &T, Biconnected &B, LCAGraph &L)
{
	bool structureChanged = false;

	int baseVertexU,baseVertexV;
	if(u==lca)
		baseVertexU = v;
	else
		baseVertexU = updateTreeEdgeValues(lca, u, v, G, T, B, L);
	
	if(v==lca)
		baseVertexV = u;
	else
		baseVertexV = updateTreeEdgeValues(lca, v, u, G, T, B, L);
	

	bool isBaseVertexU = isTreeEdge(T.parent,baseVertexU,lca);
	bool isBaseVertexV = isTreeEdge(T.parent,baseVertexV,lca);
	
	int index;
	if(isBaseVertexU)
	{
		if(__sync_fetch_and_add(&B.isBase[baseVertexU],1) == 0)
		{
			structureChanged = true;
		}
		// index = __sync_fetch_and_add(&L.headAffectedBaseVertices[lca],1);
		// L.affectedBaseVertices[lca][index] = baseVertexU;
	}

	if(isBaseVertexV)
	{
		if(__sync_fetch_and_add(&B.isBase[baseVertexV],1) == 0)
		{
			structureChanged = true;
		}
		// index = __sync_fetch_and_add(&L.headAffectedBaseVertices[lca],1);
		// L.affectedBaseVertices[lca][index] = baseVertexV;
	}

	if(!isBaseVertexU || !isBaseVertexV)
		return structureChanged;

	
	omp_set_lock(&(L.lLcaGraph[baseVertexU]));
	if(L.globalGraph.addEdge1(baseVertexU,baseVertexV))
		structureChanged = true;
	omp_unset_lock(&(L.lLcaGraph[baseVertexU]));

	omp_set_lock(&(L.lLcaGraph[baseVertexV]));
	L.globalGraph.addEdge1(baseVertexV,baseVertexU);
	omp_unset_lock(&(L.lLcaGraph[baseVertexV]));

	return structureChanged;
}

void destructLCAGraph(int lca, unweightedGraph &G, int* parent, int u, int v, int* unfinishedTraversal,
					 	vector<set<pair<int,int>> > &setOfTreeEdges, int* isLCA, vector<int> &partiallyAffectedLCA,
					 	omp_lock_t* lLcaGraph, omp_lock_t* lSetOfTreeEdges, dynamicGraph &globalGraph,
						int* isBase)
{
	int baseVertexU,baseVertexV;
	if(u==lca)
		baseVertexU = v;
	else
		baseVertexU = updateTreeEdgeValuesDelete(G,parent,lca,u,v,unfinishedTraversal,setOfTreeEdges,isLCA,
													partiallyAffectedLCA, lSetOfTreeEdges);
	if(v==lca)
		baseVertexV = u;
	else
		baseVertexV = updateTreeEdgeValuesDelete(G,parent,lca,v,u,unfinishedTraversal,setOfTreeEdges,isLCA,
													partiallyAffectedLCA, lSetOfTreeEdges);
	

	if(isTreeEdge(parent,baseVertexU,lca))
	{
		__sync_fetch_and_add(&isBase[baseVertexU],-1);
		if(isBase[baseVertexU] < 0)
			printf("Problem in isBase value\n");
	}

	if(isTreeEdge(parent,baseVertexV,lca))
	{
		__sync_fetch_and_add(&isBase[baseVertexV],-1);
		if(isBase[baseVertexV] < 0)
			printf("Problem in isBase value\n");
	}

	if(!isTreeEdge(parent,baseVertexU,lca) || !isTreeEdge(parent,baseVertexV,lca))
		return;
	
	omp_set_lock(&(lLcaGraph[baseVertexU]));
	globalGraph.deleteEdge1(baseVertexU,baseVertexV);
	omp_unset_lock(&(lLcaGraph[baseVertexU]));

	omp_set_lock(&(lLcaGraph[baseVertexV]));
	globalGraph.deleteEdge1(baseVertexV,baseVertexU);
	omp_unset_lock(&(lLcaGraph[baseVertexV]));
}

bool destructLCAGraph(int lca, int u, int v, unweightedGraph &G, SpanningTree &T, Biconnected &B, LCAGraph &L)
{
	bool structureChanged = false;

	int baseVertexU,baseVertexV;
	if(u==lca)
		baseVertexU = v;
	else
		baseVertexU = updateTreeEdgeValuesDelete(lca, u, v, G, T, B, L);
	if(v==lca)
		baseVertexV = u;
	else
		baseVertexV = updateTreeEdgeValuesDelete(lca, v, u, G, T, B, L);
	

	bool isBaseVertexU = isTreeEdge(T.parent,baseVertexU,lca);
	bool isBaseVertexV = isTreeEdge(T.parent,baseVertexV,lca);

	int index;

	if(isBaseVertexU)
	{
		if(__sync_add_and_fetch(&B.isBase[baseVertexU],-1) == 0)
		{
			structureChanged = true;
		}
	}

	if(isBaseVertexV)
	{
		if(__sync_add_and_fetch(&B.isBase[baseVertexV],-1) == 0)
		{
			structureChanged = true;
		}
	}

	if(!isBaseVertexU || !isBaseVertexV)
		return structureChanged;
	
	omp_set_lock(&(L.lLcaGraph[baseVertexU]));
	if(L.globalGraph.deleteEdge1(baseVertexU,baseVertexV))
		structureChanged = true;
	omp_unset_lock(&(L.lLcaGraph[baseVertexU]));

	omp_set_lock(&(L.lLcaGraph[baseVertexV]));
	L.globalGraph.deleteEdge1(baseVertexV,baseVertexU);
	omp_unset_lock(&(L.lLcaGraph[baseVertexV]));

	return structureChanged;
}