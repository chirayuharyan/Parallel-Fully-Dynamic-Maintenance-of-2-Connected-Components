#include"addEdge.h"
#include"lca.h"
#include<iostream>
#include<vector>
#include<set>
#include <sstream>
#include<omp.h>



void copyFromSetToVector(set<pair<int,int>> s, vector<pair<int,int>> &v)
{
	set<pair<int,int>>::iterator it;
	for(it = s.begin(); it!= s.end(); it++)
	{
		// string curr = *it;
		// stringstream currStream(curr);
		// pair<int,int> temp;
		// currStream>>temp.first>>temp.second;
		v.push_back(*it);
	}
}

void getSpecial(int &xPrev, int &xSpecial, int* parent, vector<int> &isSpecial, int currNumber, int tid)
{
	int old = xSpecial;
	while(isSpecial[xSpecial] != currNumber)
	{
		// if(tid==0 && xSpecial == 13417)
		// std::cout<<tid<<" Special Vertex = "<<xSpecial<<endl;
		xPrev = xSpecial;
		xSpecial = parent[xSpecial];
		// if(parent[xSpecial] == xSpecial && isSpecial[xSpecial] != currNumber)
		// 	std::cout<<"Problem "<<endl;
	}
}

void getSpecial(int &xPrev, int &xSpecial, int* parent, int* isSpecial, int currNumber, int tid)
{
	int old = xSpecial;
	while(isSpecial[xSpecial] != currNumber)
	{
		// if(tid==0 && xSpecial == 13417)
		// std::cout<<tid<<" Special Vertex = "<<xSpecial<<endl;
		xPrev = xSpecial;
		xSpecial = parent[xSpecial];
		// if(parent[xSpecial] == xSpecial && isSpecial[xSpecial] != currNumber)
		// 	std::cout<<"Problem "<<endl;
	}
}

void adjustLCAGraph(int oldLCA, int oldBase1, int oldBase2, int newLCA, int newBase1, int newBase2,
					 int* unfinishedTraversal, int* isLCA, int* parent, omp_lock_t* lLcaGraph,
					 dynamicGraph &globalGraph,int x,int y, int* isBase)
{
	bool isOldBase1,isOldBase2,isNewBase1,isNewBase2;
	isOldBase1 = isTreeEdge(parent,oldBase1,oldLCA);
	isOldBase2 = isTreeEdge(parent,oldBase2,oldLCA);
	isNewBase1 = isTreeEdge(parent,newBase1,newLCA);
	isNewBase2 = isTreeEdge(parent,newBase2,newLCA);


	if(isOldBase1)
	{
		__sync_fetch_and_add(&isBase[oldBase1],-1);
	}

	if(isOldBase2)
	{
		__sync_fetch_and_add(&isBase[oldBase2],-1);
	}

	if(isOldBase1 && isOldBase2)
	{
		
		omp_set_lock(&lLcaGraph[oldBase1]);
		globalGraph.deleteEdge1(oldBase1,oldBase2);
		omp_unset_lock(&lLcaGraph[oldBase1]);

		omp_set_lock(&lLcaGraph[oldBase2]);
		globalGraph.deleteEdge1(oldBase2,oldBase1);
		omp_unset_lock(&lLcaGraph[oldBase2]);

	}


	if(isNewBase1)
	{
		__sync_fetch_and_add(&isBase[newBase1],1);
	}

	if(isNewBase2)
	{
		__sync_fetch_and_add(&isBase[newBase2],1);
	}


	if(isNewBase1 && isNewBase2)
	{
	
		omp_set_lock(&lLcaGraph[newBase1]);
		globalGraph.addEdge1(newBase1,newBase2);
		omp_unset_lock(&lLcaGraph[newBase1]);

		omp_set_lock(&lLcaGraph[newBase2]);
		globalGraph.addEdge1(newBase2,newBase1);
		omp_unset_lock(&lLcaGraph[newBase2]);	
	}



	__sync_fetch_and_add(&isLCA[oldLCA],-1);
	
	__sync_fetch_and_add(&isLCA[newLCA],1);

	
	if(isOldBase1)
	{
		__sync_fetch_and_add(&unfinishedTraversal[oldBase1],1);
	}

	if(isOldBase2)
	{
		__sync_fetch_and_add(&unfinishedTraversal[oldBase2],1);
	}

	if(newLCA == parent[newBase1])
	{
		__sync_fetch_and_add(&unfinishedTraversal[newBase1],-1);
	}
	else if(newBase1 == parent[newLCA])
	{
		__sync_fetch_and_add(&unfinishedTraversal[newLCA],-1);
	}

	if(newLCA == parent[newBase2])
	{
		__sync_fetch_and_add(&unfinishedTraversal[newBase2],-1);
	}
	else if(newBase2 == parent[newLCA])
	{
		__sync_fetch_and_add(&unfinishedTraversal[newLCA],-1);
	}

}

void adjustLCAGraph(int oldLCA, int oldBase1, int oldBase2, int newLCA, int newBase1, int newBase2,
					int x,int y, SpanningTree &T, Biconnected &B, LCAGraph &L)
{
	bool isOldBase1,isOldBase2,isNewBase1,isNewBase2;
	isOldBase1 = isTreeEdge(T.parent,oldBase1,oldLCA);
	isOldBase2 = isTreeEdge(T.parent,oldBase2,oldLCA);
	isNewBase1 = isTreeEdge(T.parent,newBase1,newLCA);
	isNewBase2 = isTreeEdge(T.parent,newBase2,newLCA);


	if(isOldBase1)
	{
		__sync_fetch_and_add(&B.isBase[oldBase1],-1);
	}

	if(isOldBase2)
	{
		__sync_fetch_and_add(&B.isBase[oldBase2],-1);
	}

	if(isOldBase1 && isOldBase2)
	{
		
		omp_set_lock(&L.lLcaGraph[oldBase1]);
		L.globalGraph.deleteEdge1(oldBase1,oldBase2);
		omp_unset_lock(&L.lLcaGraph[oldBase1]);

		omp_set_lock(&L.lLcaGraph[oldBase2]);
		L.globalGraph.deleteEdge1(oldBase2,oldBase1);
		omp_unset_lock(&L.lLcaGraph[oldBase2]);

	}


	if(isNewBase1)
	{
		__sync_fetch_and_add(&B.isBase[newBase1],1);
	}

	if(isNewBase2)
	{
		__sync_fetch_and_add(&B.isBase[newBase2],1);
	}


	if(isNewBase1 && isNewBase2)
	{
	
		omp_set_lock(&L.lLcaGraph[newBase1]);
		L.globalGraph.addEdge1(newBase1,newBase2);
		omp_unset_lock(&L.lLcaGraph[newBase1]);

		omp_set_lock(&L.lLcaGraph[newBase2]);
		L.globalGraph.addEdge1(newBase2,newBase1);
		omp_unset_lock(&L.lLcaGraph[newBase2]);	
	}



	__sync_fetch_and_add(&B.isLCA[oldLCA],-1);
	
	__sync_fetch_and_add(&B.isLCA[newLCA],1);

	
	if(isOldBase1)
	{
		__sync_fetch_and_add(&B.unfinishedTraversal[oldBase1],1);
	}

	if(isOldBase2)
	{
		__sync_fetch_and_add(&B.unfinishedTraversal[oldBase2],1);
	}

	if(newLCA == T.parent[newBase1])
	{
		__sync_fetch_and_add(&B.unfinishedTraversal[newBase1],-1);
	}
	else if(newBase1 == T.parent[newLCA])
	{
		__sync_fetch_and_add(&B.unfinishedTraversal[newLCA],-1);
	}

	if(newLCA == T.parent[newBase2])
	{
		__sync_fetch_and_add(&B.unfinishedTraversal[newBase2],-1);
	}
	else if(newBase2 == T.parent[newLCA])
	{
		__sync_fetch_and_add(&B.unfinishedTraversal[newLCA],-1);
	}

}

void updateNonTreeEdge(int x, int y, int* parent, vector<int> &isSpecial, vector<int> &rank, int currNumber,
					 vector<int> &P, int* unfinishedTraversal, int* isLCA, omp_lock_t* lLcaGraph,
					 int tid, dynamicGraph &globalGraph, int* isBase)
{
	
	int xPrev = y, xSpecial = x, yPrev = x, ySpecial = y;
	getSpecial(xPrev, xSpecial, parent, isSpecial, currNumber,tid);
	
	getSpecial(yPrev, ySpecial, parent, isSpecial, currNumber,tid);
	
	
	if(rank[xSpecial] > rank[ySpecial])
	{
		adjustLCAGraph(xSpecial, xPrev, P[rank[xSpecial]-1], ySpecial, yPrev, P[rank[ySpecial]+1],
						unfinishedTraversal, isLCA, parent, lLcaGraph, globalGraph, x, y, isBase);
	}
	else
	{
		adjustLCAGraph(ySpecial, yPrev, P[rank[ySpecial]-1], xSpecial, xPrev, P[rank[xSpecial]+1],
			 			unfinishedTraversal, isLCA, parent, lLcaGraph, globalGraph, x, y, isBase);
	}


	// affectedLCA[xSpecial] = true;
	// affectedLCA[ySpecial] = true;
}

void updateNonTreeEdge(int x, int y, vector<int> &P, SpanningTree &T, Biconnected &B, LCAGraph &L, int tid,
						int currNumber)
{
	int* isSpecial = B.isSpecial[tid];
	int* rank = B.rank[tid];

	
	int xPrev = y, xSpecial = x, yPrev = x, ySpecial = y;
	
	getSpecial(xPrev, xSpecial, T.parent, isSpecial, currNumber,tid);
	
	getSpecial(yPrev, ySpecial, T.parent, isSpecial, currNumber,tid);
	
	
	if(rank[xSpecial] > rank[ySpecial])
	{
		// adjustLCAGraph(xSpecial, xPrev, P[rank[xSpecial]-1], ySpecial, yPrev, P[rank[ySpecial]+1],
		// 				unfinishedTraversal, isLCA, parent, lLcaGraph, globalGraph, x, y, isBase);

		adjustLCAGraph(xSpecial, xPrev, P[rank[xSpecial]-1], ySpecial, yPrev, P[rank[ySpecial]+1],
						x,y,T,B,L);
	}
	else
	{
		// adjustLCAGraph(ySpecial, yPrev, P[rank[ySpecial]-1], xSpecial, xPrev, P[rank[xSpecial]+1],
		// 	 			unfinishedTraversal, isLCA, parent, lLcaGraph, globalGraph, x, y, isBase);

		adjustLCAGraph(ySpecial, yPrev, P[rank[ySpecial]-1], xSpecial, xPrev, P[rank[xSpecial]+1],
						x,y,T,B,L);

	}


	// affectedLCA[xSpecial] = true;
	// affectedLCA[ySpecial] = true;
}


void addNonTreeEdge(int u,int v, unweightedGraph &G, int* parent, int* isLCA, int* unfinishedTraversal, 
					vector<set<pair<int,int>> > &setOfTreeEdges, omp_lock_t* lLcaGraph,
					omp_lock_t* lSetOfTreeEdges, vector<int> &result, vector<int> &partiallyAffectedLCA, 
					vector<int> &isVisited, int &findLCACount, dynamicGraph &globalGraph, int* isBase)
{
	
	if(u>v)
	{
		int temp = u;
		u = v;
		v = temp;
	}

	int LCAVertex = findLCA(G,parent,u,v,isVisited,findLCACount);

	__sync_fetch_and_add(&isLCA[LCAVertex],1);


	constructLCAGraph(LCAVertex, G, parent, u, v, unfinishedTraversal, setOfTreeEdges, isLCA,
						partiallyAffectedLCA, lLcaGraph, lSetOfTreeEdges, globalGraph, isBase);
	
	result.push_back(LCAVertex);
}

void addNonTreeEdge(int u,int v, unweightedGraph &G, SpanningTree &T, Biconnected &B, LCAGraph &L,int tid)
{
	
	if(u>v)
	{
		int temp = u;
		u = v;
		v = temp;
	}

	int LCAVertex = findLCA(u,v,G,T,B,tid);

	__sync_fetch_and_add(&(B.isLCA[LCAVertex]),1);

	bool structureChanged = constructLCAGraph(LCAVertex, u, v, G, T, B, L);
	
	if(structureChanged)
	{
		int head = __sync_fetch_and_add(&(L.headAffectedLCA),1);
		L.affectedLCA[head] = LCAVertex;
	}
	
}

void addTreeEdge(int u,int v, int count,unweightedGraph &G, SpanningTree &T, Biconnected &B, LCAGraph &L,int tid)
{
	double time1,time2;

	// printf("***1 - Add tree edge Started\n");
	int *isSpecial = B.isSpecial[tid];
	int* rank = B.rank[tid];
	int* parent = T.parent;

	time1=omp_get_wtime();
	vector<int> P; //special Path
	vector<pair<int,int>> affectedNonTreeEdges;

	int curr = v;
	P.push_back(v);
	
	isSpecial[curr] = count;
	set<pair<int,int>> affectedNonTreeEdgesSet;

	if(parent[curr]!=curr)
		affectedNonTreeEdgesSet = B.setOfTreeEdges[curr];
	
	int currRank = 0;
	rank[curr] = currRank;
	// affectedLCA[curr] = false; //probably this is not required

	set<pair<int,int>>::iterator it;
	while(parent[curr] != curr)
	{
		curr = parent[curr];
		P.push_back(curr);
		isSpecial[curr] = count;
		for(it = B.setOfTreeEdges[curr].begin(); it!= B.setOfTreeEdges[curr].end();it++)
		{
			affectedNonTreeEdgesSet.insert(*it);
		}
		
		currRank++;
		rank[curr] = currRank;
		// affectedLCA[curr] = false;
	}

	time2 = omp_get_wtime();
	// printf("Time For Step 1 =  %f ms with size of P = %d amd affectedNonTreeEdges size = %d\n",(time2-time1)*1000,int(P.size()),int(affectedNonTreeEdgesSet.size()));


	//step 1 : get the path from v to root and create following data structure
		//datastructure : list of integer in order v .... root, number rank to each spl vertex.
		//datastructure 2: list of non-tree edges which are affected

	time1=omp_get_wtime();
	copyFromSetToVector(affectedNonTreeEdgesSet,affectedNonTreeEdges);
	// #pragma omp barrier	//as add tree edge is happining in parallel no need to make this parallel 
	for(int i=0; i<affectedNonTreeEdges.size(); i++)
	{
		int x = affectedNonTreeEdges[i].first;
		int y = affectedNonTreeEdges[i].second;


		//update from here
		updateNonTreeEdge(x,y,P,T,B,L,tid,count);

	}
	time2 = omp_get_wtime();
	// printf("Time For Step 2 =  %f ms\n",(time2-time1)*1000);
	//step 2: for each non tree edge (u,v) identify the first ancestor for
			// each endpoint of the edge which is a spl vertice


	//step 3: update the lca graph by taking the old base-vertices and new base vertices

	//repeat step 2 and 3 for all non-tree edges
	
	time1=omp_get_wtime();
	parent[P[0]] = P[0];
	for(int i=1; i<P.size(); i++)
	{
		parent[P[i]] = P[i-1];
	}

	for(int i= P.size()-1; i>0;i--)
	{
		B.setOfTreeEdges[P[i]] = B.setOfTreeEdges[P[i-1]];
		B.unfinishedTraversal[P[i]] = B.unfinishedTraversal[P[i-1]];
	}

	B.unfinishedTraversal[P[0]] = 0;
	B.setOfTreeEdges[P[0]].clear();

	// std::cout<<"Step 4 completed "<<tid<<" "<<count<<endl;
	//step 4: reverse the parent - child pointers and also transfer set of tree edges

	parent[v] = u;

	// std::cout<<"Step 5 completed "<<tid<<" "<<count<<endl;
	// step 5: make parent of v as u 

	// count++; //count is updated after the function is called

	int indexHead = __sync_fetch_and_add(&(L.headAffectedLCA),P.size());
	for(int i = 0;i<P.size();i++)
		L.affectedLCA[indexHead+i] = P[i];

	time2 = omp_get_wtime();
	// printf("Time For Step 2 =  %f ms\n",(time2-time1)*1000);

}



void addTreeEdge(int u,int v, unweightedGraph &G, int* parent,int* isLCA, int* unfinishedTraversal,
					vector<set<pair<int,int>> > &setOfTreeEdges, omp_lock_t* lLcaGraph,
					omp_lock_t* lSetOfTreeEdges, vector<int> &result,
					int count, vector<int> &isSpecial, vector<int> &rank, int tid, dynamicGraph &globalGraph, 
					int* isBase)
{
	
	vector<int> P; //special Path
	vector<pair<int,int>> affectedNonTreeEdges;

	int curr = v;
	P.push_back(v);
	
	isSpecial[curr] = count;
	set<pair<int,int>> affectedNonTreeEdgesSet;

	if(parent[curr]!=curr)
		affectedNonTreeEdgesSet = setOfTreeEdges[curr];
	
	int currRank = 0;
	rank[curr] = currRank;
	// affectedLCA[curr] = false; //probably this is not required

	set<pair<int,int>>::iterator it;
	while(parent[curr] != curr)
	{
		curr = parent[curr];
		P.push_back(curr);
		isSpecial[curr] = count;
		for(it = setOfTreeEdges[curr].begin(); it!= setOfTreeEdges[curr].end();it++)
		{
			affectedNonTreeEdgesSet.insert(*it);
		}
		
		currRank++;
		rank[curr] = currRank;
		// affectedLCA[curr] = false;
	}


	//step 1 : get the path from v to root and create following data structure
		//datastructure : list of integer in order v .... root, number rank to each spl vertex.
		//datastructure 2: list of non-tree edges which are affected


	copyFromSetToVector(affectedNonTreeEdgesSet,affectedNonTreeEdges);
	// #pragma omp barrier	//as add tree edge is happining in parallel no need to make this parallel 
	for(int i=0; i<affectedNonTreeEdges.size(); i++)
	{
		int x = affectedNonTreeEdges[i].first;
		int y = affectedNonTreeEdges[i].second;
		updateNonTreeEdge(x,y,parent,isSpecial,rank,count,P,unfinishedTraversal,isLCA,
			lLcaGraph,tid,globalGraph,isBase);

	}

	//step 2: for each non tree edge (u,v) identify the first ancestor for
			// each endpoint of the edge which is a spl vertice


	//step 3: update the lca graph by taking the old base-vertices and new base vertices

	//repeat step 2 and 3 for all non-tree edges
	
	parent[P[0]] = P[0];
	for(int i=1; i<P.size(); i++)
	{
		parent[P[i]] = P[i-1];
	}

	for(int i= P.size()-1; i>0;i--)
	{
		setOfTreeEdges[P[i]] = setOfTreeEdges[P[i-1]];
		unfinishedTraversal[P[i]] = unfinishedTraversal[P[i-1]];
	}

	unfinishedTraversal[P[0]] = 0;
	setOfTreeEdges[P[0]].clear();

	// std::cout<<"Step 4 completed "<<tid<<" "<<count<<endl;
	//step 4: reverse the parent - child pointers and also transfer set of tree edges

	parent[v] = u;

	// std::cout<<"Step 5 completed "<<tid<<" "<<count<<endl;
	// step 5: make parent of v as u 

	// count++; //count is updated after the function is called

	for(int i = 0;i<P.size();i++)
		result.push_back(P[i]);

}



