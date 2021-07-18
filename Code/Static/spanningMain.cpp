#include<omp.h>
#include<iostream>
#include"unweightedGraph.h"
#include"spanningTree.h"

using namespace std;
#define NUM_THREADS 8

int main(int argc, char* argv[])
{
	int totalThreads = NUM_THREADS;
	omp_set_num_threads(totalThreads);

	ifstream inputGraph;
	inputGraph.open(argv[1]);
	unweightedGraph G(inputGraph);

	vector<int> root;
	
	int *parent = new int[G.totalVertices];
	omp_lock_t* lParent = new omp_lock_t[G.totalVertices];

	bool* spanningTreeVisited = new bool[G.totalVertices];
	omp_lock_t* lSpanningTreeVisited = new omp_lock_t[G.totalVertices];

	for(int i=0;i<G.totalVertices;i++)
	{
		spanningTreeVisited[i] = false;
		omp_init_lock(&(lParent[i]));
		omp_init_lock(&(lSpanningTreeVisited[i]));
	}

	double time1,time2;
	
	time1=omp_get_wtime();
	spanningTree(G,parent,root,spanningTreeVisited,lSpanningTreeVisited);
	time2=omp_get_wtime();
	cout<< "Time taken for parallel frontier spanning tree: "<< ((time2-time1)*1000) << " milliseconds" << endl;

	for(int i=0;i<G.totalVertices;i++)
	{
		spanningTreeVisited[i] = false;
		parent[i] = -1;
	}
	
	
}