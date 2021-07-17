#include"deleteEdge.h"
#include"addEdge.h"
#include"threads.h"
#include"lca.h"
#include<iostream>
#include<set>
#include<vector>
#include <sstream>
#include<omp.h>
#include<queue>
#include <chrono> 

using namespace std::chrono;


void insertEdge(int u, int v, map<pair<int,int>,pair<int,int>> &sGEdge, vector<pair<int,int>> &newTreeEdge)
{

	//u has to be first and v has to be second
	pair<int,int> curr;
	pair<int,int> next;
	if(u>v)
	{
		int temp = u;
		u = v;
		v = temp;
		curr.first = u;
		curr.second = v;
		next = sGEdge[curr];
		temp = next.first;
		next.first = next.second;
		next.second = temp;
	}
	else
	{	
		curr.first = u;
		curr.second = v;
		next = sGEdge[curr];
	}
	newTreeEdge.push_back(next);
	
}

int spanningTree(dynamicGraph &G, vector<pair<int,int>> &newTreeEdge, map<pair<int,int>,pair<int,int>> &sGEdge, 
	omp_lock_t &lMap)
{
	// G.print();
	map<int, map<int,unweightedEdge>>::iterator it;
	map<int,unweightedEdge>::iterator jt;
	it = G.adjacencyList.begin();
	
	queue<int> q;
	set<int> visited;
	int source = it->first;
	q.push(source);
	visited.insert(source);
	// std::cout<<"Starting BFS"<<endl;
	
	while(q.size()>0)
	{
		int curr = q.front(); q.pop();
		it = G.adjacencyList.find(curr);

		for(jt = it->second.begin(); jt!= it->second.end(); jt++)
		{
			int next = jt->first;
			if(visited.find(next) == visited.end())
			{
				visited.insert(next);
				insertEdge(curr,next,sGEdge,newTreeEdge);
				q.push(next);
			}
		}
	}
	return source;
}

int findRoot(int u, int* parent)
{
	while(parent[u] != u)
		u = parent[u];
	return u;
}

void deleteNonTreeEdge(int u,int v, unweightedGraph &G, int* parent, int* isLCA, int* unfinishedTraversal,
						vector<set<pair<int,int>> > &setOfTreeEdges, omp_lock_t* lLcaGraph,
						omp_lock_t* lSetOfTreeEdges, vector<int> &result,vector<int> &partiallyAffectedLCA, 
						vector<int> &isVisited, int &findLCACount, dynamicGraph &globalGraph, int* isBase)
{

	if(u>v)
	{
		int temp = u;
		u = v;
		v = temp;
	}

	int LCAVertex = findLCA(G,parent,u,v,isVisited,findLCACount);

	__sync_fetch_and_add(&isLCA[LCAVertex],-1);

	destructLCAGraph(LCAVertex, G, parent, u, v, unfinishedTraversal, setOfTreeEdges, isLCA, partiallyAffectedLCA,
					lLcaGraph, lSetOfTreeEdges, globalGraph, isBase);
	
	result.push_back(LCAVertex);
}

void deleteNonTreeEdge(int u,int v, unweightedGraph &G, SpanningTree &T, Biconnected &B, LCAGraph &L,int tid)
{

	if(u>v)
	{
		int temp = u;
		u = v;
		v = temp;
	}

	int LCAVertex = findLCA(u,v,G,T,B,tid);

	__sync_fetch_and_add(&B.isLCA[LCAVertex],-1);

	bool structureChanged = destructLCAGraph(LCAVertex, u, v, G, T, B, L);
	
	if(structureChanged)
	{
		int head = __sync_fetch_and_add(&(L.headAffectedLCA),1);
		L.affectedLCA[head] = LCAVertex;
	}
}



void deleteEdgeParallel(unweightedGraph &G, int count, int totalQueries, vector<query> &queries,
						  SpanningTree &T, Biconnected &B, LCAGraph &L)
{
	int* parent = T.parent;

	//Step 1: Remove Non-tree edges and exclude tree edges
	#pragma omp parallel for schedule(guided)
	for( int i= count ; i< (count+totalQueries) ; i++ )
	{		
		int u = queries[i].u;
		int v = queries[i].v;
		int tid = omp_get_thread_num();
		
		if(parent[u] == v || parent[v] == u)
		{
			continue;
		}

		deleteNonTreeEdge(u,v,G,T,B,L,tid);
						 //adds edge and updates the lca graph
	}
	
	
	

	
	// Step 2 - Collect nontree edges from every tree edge to be deleted.

	set<pair<int,int>> commonBucket;
	vector<pair<int,int>> vecCommonBucket;		
	// omp_lock_t lCommonBucket;
	// omp_init_lock (&(lCommonBucket));
	

	//Step 2.1 - Collect Support of each tree edge to be deleted in a set
	int totalTreeQueries = 0;
	// #pragma omp parallel for schedule(guided)
	for( int i= count ; i< (count+totalQueries) ; i++ )
	{		
		int u = queries[i].u;
		int v = queries[i].v;
		// int tid = omp_get_thread_num();		
		if(parent[u] != v && parent[v] != u) continue;
		__sync_fetch_and_add(&totalTreeQueries,1);
		if(parent[u] == v)
		{
			int temp = u;
			u = v;
			v = temp;
		}
		// u is parent of v
		
		// omp_set_lock(&(lCommonBucket));
		commonBucket.insert(B.setOfTreeEdges[v].begin(), B.setOfTreeEdges[v].end());
		// omp_unset_lock(&(lCommonBucket));
	}

	// printf("Time For collecting support of tree edges to be deleted =  %f ms\n",(time2-time1)*1000);


	if(totalTreeQueries == 0) return;

	// printf("Total Tree edges = %d and Total Support = %d \n",totalTreeQueries,int(commonBucket.size()));
	
	// Step 2.2 - Copy the support set into a list
	copyFromSetToVector(commonBucket, vecCommonBucket);
	


	// Step 2.3 - Delete All the supporting edges.

	#pragma omp parallel for schedule(guided)
	for(int i = 0; i < vecCommonBucket.size(); i++)
	{
		int u = vecCommonBucket[i].first;
		int v = vecCommonBucket[i].second;
		int tid = omp_get_thread_num();
		
		deleteNonTreeEdge(u,v,G,T,B,L,tid);
	}
	

	//Step 3 - Remove tree edges by updating the parent array



	#pragma omp parallel for
	for( int i= count ; i< (count+totalQueries) ; i++ )
	{		
		int u = queries[i].u;
		int v = queries[i].v;
		
		if(parent[u] != v && parent[v] != u) continue;

		if(parent[u] == v)
		{
			int temp = u;
			u = v;
			v = temp;
		}
		// u is parent of v
		parent[v] = v;
	}
	

	//Step 4 - Create supergraph by identifying the end points of each non-tree edge in parallel
	
	dynamicGraph sG;
	omp_lock_t lSG;
	omp_init_lock (&(lSG));

	map<pair<int,int>,pair<int,int>> sGEdge;
	omp_lock_t lMap;
	omp_init_lock (&(lMap));

	#pragma omp parallel for
	for(int i = 0; i<vecCommonBucket.size();i++)
	{
		int u = vecCommonBucket[i].first;
		int v = vecCommonBucket[i].second;
		int rootU = findRoot(u,parent);
		int rootV = findRoot(v,parent);

		omp_set_lock(&lSG);
		bool flag = sG.addEdge(rootU, rootV);
		omp_unset_lock(&lSG);
		
		if(flag)
		{
			pair<int,int> key,value;
			key.first = rootU;
			key.second = rootV;
			value.first = u;
			value.second = v;
			if(rootU > rootV)
			{
				key.first = rootV;
				key.second = rootU;
				value.first = v;
				value.second = u;
			}
			omp_set_lock(&lMap);
			sGEdge[key] = value;
			omp_unset_lock(&lMap);
		}
	}
	
	//Step 5 - Create Spanning tree on supergraph using BFS (serially for now)

	vector<pair<int,int>> newTreeEdge;
	int newRoot = spanningTree(sG,newTreeEdge,sGEdge,lMap);
	// printf("newRoot = %d\n",newRoot);
	G.root = newRoot;

	//Step 6 - Add edges in the spanning tree
	#pragma omp parallel for
	for(int i=0;i<newTreeEdge.size();i++)
	{
		int u = newTreeEdge[i].first;
		int v = newTreeEdge[i].second;
		int tid = omp_get_thread_num();
	
		addTreeEdge(u, v, B.countTreeEdge+i, G, T, B, L, tid);
	}

	

	B.countTreeEdge = B.countTreeEdge + newTreeEdge.size();
	

	// Step 7 - Add the remaining support edges 

	#pragma omp parallel for
	for(int i = 0; i < vecCommonBucket.size(); i++)
	{
		int u = vecCommonBucket[i].first;
		int v = vecCommonBucket[i].second;
		int tid = omp_get_thread_num();
			
		if(u == parent[v] || v == parent[u]) continue;

		addNonTreeEdge(u,v,G,T,B,L,tid);

	}

}	

void deleteEdgeParallel(unweightedGraph &G, int* parent,int* isLCA, int* unfinishedTraversal, 
						vector<set<pair<int,int>> > &setOfTreeEdges, int count, int totalQueries,
						vector<query> &queries, omp_lock_t* lLcaGraph, omp_lock_t* lSetOfTreeEdges, 
						vector<vector<int>> &affectedLCA, vector<vector<int>> &partiallyAffectedLCA, 
						int &countTreeEdge, vector<vector<int>> &isSpecial, vector<vector<int>> &rank,
						vector<vector<int>> &isVisited, int* findLCACount, dynamicGraph &globalGraph, int* isBase)
{


	#pragma omp parallel for
	for( int i= count ; i< (count+totalQueries) ; i++ )
	{		
		int u = queries[i].u;
		int v = queries[i].v;
		int tid = omp_get_thread_num();
		
		if(parent[u] == v || parent[v] == u)
		{
			continue;
		}

		deleteNonTreeEdge(u,v,G,parent,isLCA,unfinishedTraversal,setOfTreeEdges, lLcaGraph,lSetOfTreeEdges,
						affectedLCA[tid], partiallyAffectedLCA[tid],isVisited[tid], findLCACount[tid],globalGraph,
						isBase);
						 //adds edge and updates the lca graph
	}

	//deleting all the tree edges
	set<pair<int,int>> commonBucket;
	vector<pair<int,int>> vecCommonBucket;		
	omp_lock_t lCommonBucket;
	omp_init_lock (&(lCommonBucket));


	//Phase 1
	//collect nontree edges from every tree edge to be deleted
	

	int totalTreeQueries = 0;
	#pragma omp parallel for
	for( int i= count ; i< (count+totalQueries) ; i++ )
	{		
		int u = queries[i].u;
		int v = queries[i].v;
		int tid = omp_get_thread_num();		
		if(parent[u] != v && parent[v] != u) continue;
		__sync_fetch_and_add(&totalTreeQueries,1);
		if(parent[u] == v)
		{
			int temp = u;
			u = v;
			v = temp;
		}
		// u is parent of v
		
		omp_set_lock(&(lCommonBucket));
		commonBucket.insert(setOfTreeEdges[v].begin(), setOfTreeEdges[v].end());
		omp_unset_lock(&(lCommonBucket));
	}

	if(totalTreeQueries == 0) return;
	
	copyFromSetToVector(commonBucket, vecCommonBucket);

	#pragma omp parallel for
	for(int i = 0; i < vecCommonBucket.size(); i++)
	{
		int u = vecCommonBucket[i].first;
		int v = vecCommonBucket[i].second;
		int tid = omp_get_thread_num();
		
		deleteNonTreeEdge(u,v,G,parent,isLCA,unfinishedTraversal,setOfTreeEdges, lLcaGraph,lSetOfTreeEdges,
							affectedLCA[tid],partiallyAffectedLCA[tid], isVisited[tid], findLCACount[tid],
							globalGraph, isBase);
	}


	//Phase 2
	
	//remove tree edges by updating the parent array


	#pragma omp parallel for
	for( int i= count ; i< (count+totalQueries) ; i++ )
	{		
		int u = queries[i].u;
		int v = queries[i].v;
		
		if(parent[u] != v && parent[v] != u) continue;

		if(parent[u] == v)
		{
			int temp = u;
			u = v;
			v = temp;
		}
		// u is parent of v
		parent[v] = v;
	}

	//create supergraph by identifying the end points of each non-tree edge in parallel
	
	// start2 = high_resolution_clock::now();
	dynamicGraph sG;
	omp_lock_t lSG;
	omp_init_lock (&(lSG));

	map<pair<int,int>,pair<int,int>> sGEdge;
	omp_lock_t lMap;
	omp_init_lock (&(lMap));

	#pragma omp parallel for
	for(int i = 0; i<vecCommonBucket.size();i++)
	{
		int u = vecCommonBucket[i].first;
		int v = vecCommonBucket[i].second;
		int rootU = findRoot(u,parent);
		int rootV = findRoot(v,parent);

		omp_set_lock(&lSG);
		bool flag = sG.addEdge(rootU, rootV);
		omp_unset_lock(&lSG);
		
		if(flag)
		{
			pair<int,int> key,value;
			key.first = rootU;
			key.second = rootV;
			value.first = u;
			value.second = v;
			if(rootU > rootV)
			{
				key.first = rootV;
				key.second = rootU;
				value.first = v;
				value.second = u;
			}
			omp_set_lock(&lMap);
			sGEdge[key] = value;
			omp_unset_lock(&lMap);
		}
	}
	
	//spanning tree on supergraph using BFS (serially for now)

	vector<pair<int,int>> newTreeEdge;
	spanningTree(sG,newTreeEdge,sGEdge,lMap);


	#pragma omp parallel for
	for(int i=0;i<newTreeEdge.size();i++)
	{
		int u = newTreeEdge[i].first;
		int v = newTreeEdge[i].second;
		int tid = omp_get_thread_num();
	
		addTreeEdge(u,v,G,parent, isLCA, unfinishedTraversal,setOfTreeEdges, lLcaGraph,lSetOfTreeEdges,
					affectedLCA[tid], countTreeEdge+i, isSpecial[tid], rank[tid], tid,globalGraph, isBase);
	}

	countTreeEdge = countTreeEdge+newTreeEdge.size();
	
	// Phase 3


	#pragma omp parallel for
	for(int i = 0; i < vecCommonBucket.size(); i++)
	{
		int u = vecCommonBucket[i].first;
		int v = vecCommonBucket[i].second;
		int tid = omp_get_thread_num();
			
		if(u == parent[v] || v == parent[u]) continue;

		addNonTreeEdge(u,v,G,parent,isLCA,unfinishedTraversal,setOfTreeEdges, 
					lLcaGraph,lSetOfTreeEdges,affectedLCA[tid],partiallyAffectedLCA[tid],
					isVisited[tid], findLCACount[tid],globalGraph, isBase);

	}

}	