#include <omp.h>
#include"unweightedGraph.h"
#include"dynamicGraph.h"
#include"bfs.h"
#include"spanningTree.h"
#include"lca.h"
#include"isCutVertex.h"
#include"addEdge.h"
#include"deleteEdge.h"
#include"threads.h"
#include<iostream>
#include<set>
#include<algorithm>
#include <chrono> 
#include<cfloat>

// #define NUM_THREADS 16

using namespace std;
using namespace std::chrono;

void postProcessing(unweightedGraph &G, int* parent, vector<set<pair<int,int>>> &setOfTreeEdges, int* isLCA,
 					bool* cutVertex, bool* cutVertexAsBridge, ofstream &output)
{
	set<int> cutVertices;
	vector<string> bridges; 
	output<<"Bridges"<<endl;
	int totalBridges = 0;
	
	for(int i=0;i<G.totalVertices;i++)
	{
		int u = i;
		int v = parent[u];

		if(u==v) continue;

		if(setOfTreeEdges[u].size() > 0)
			continue;
		
		totalBridges++;
			/*if(G.offset[u+1] - G.offset[u] > 1)
				cutVertices.insert(u);
			if(G.offset[v+1] - G.offset[v] > 1)
				cutVertices.insert(v);*/

		// if(G.degree[u] > 1)
		// 	cutVertices.insert(u);
		// if(G.degree[v] > 1)
		// 	cutVertices.insert(v);
			

			/*if(G.adjacencyList[u].size() > 1)
			{
				cutVertices.insert(u);
			}
			if(G.adjacencyList[v].size() > 1)
			{
				cutVertices.insert(v);
			}*/
			
		if(u>v)
		{
			int temp = u;
			u = v;
			v = temp;
		}

		bridges.push_back(to_string(u) +" "+ to_string(v));
			// cout<<u<<" "<<v<<endl;
	}
		
	


	sort(bridges.begin(), bridges.end());
	for(int i=0;i<bridges.size();i++)
	{
		output<<bridges[i]<<endl;
	}

	int totalCutVertices = 0;
	output<<"Cut Vertices"<<endl;
	for(int i=0;i<G.totalVertices;i++)
	{
		if(isLCA[i] == 0 && cutVertex[i])
		{
			printf("%d is not lca but cutVertex is true\n",i);
		}
		if((isLCA[i] > 0 && cutVertex[i]) || cutVertexAsBridge[i])
		{
			totalCutVertices++;
			// cout<<i<<endl;
			cutVertices.insert(i);
		}

	}

	for(set<int>::iterator it = cutVertices.begin(); it!= cutVertices.end();it++)
	{
		output<<*it<<endl;
	}
	// cout<<endl;
	// cout<<totalBridges<<" "<<cutVertices.size()<<endl;
}


void biconnected(unweightedGraph &G, SpanningTree &T, Biconnected &B, LCAGraph &L)
{
	
	spanningTree(G,T.parent);
	

	#pragma omp parallel for schedule(guided)
	for(int i=0;i<G.totalVertices;i++)
	{
		int u = i;
		for(int j = G.offset[i]; j < G.offset[i+1]; j++)
		{
			int v = G.neighbour[j];
			int tid = omp_get_thread_num();
			if(T.parent[u] == v || T.parent[v] == u) continue;
			if(u>v) continue;
			
			addNonTreeEdge(u,v,G,T,B,L,tid); 
								//adds edge and updates the lca graph

		}
		
	}

	

	isCutVertexAll(L.globalGraph, B.unfinishedTraversal, T.parent, B.isLCA, B.cutVertex, G.totalVertices, G.root,
					B.isBase, L.localParent, L.isSafe, L.level);
	

	#pragma omp parallel for schedule(guided)
	for(int i=0;i<G.totalVertices;i++)
	{
		int u = i;
		int v = T.parent[u];

		if(u==v) continue;

		if(B.setOfTreeEdges[u].size()>0)
			continue;
		
		if(G.degree[u] > 1)
			B.cutVertexAsBridge[u] = true;
		if(G.degree[v] > 1)
			B.cutVertexAsBridge[v] = true;
		
	}

}




void readQuery(ifstream &queryFile, vector<query> &queries, int &totalBatch)
{
	queryFile>>totalBatch;
	for(int i = 0;i<totalBatch;i++)
	{
		query data,curr;
		queryFile>>data.u>>data.v;
		queries.push_back(data);
		for(int j=0;j<data.v;j++)
		{
			queryFile>>curr.u>>curr.v;
			queries.push_back(curr);
		}
	}
}


string getBatchPercent(char *input)
{
	int i=0;
	while(input[i] != '\0' && input[i] != '_')
	{
		i++;
	}
	while(input[i] != '\0' && !(input[i]>='0' && input[i]<='9'))
	{
		i++;
	}
	string result = "";
	bool firstChar = true;
	while(input[i]!='\0' && input[i]>='0' && input[i]<='9')
	{
		if(firstChar)
		{
			firstChar = false;
			result += input[i];
			if(input[i] == '0')
				result += '.';
		}
		else
		{
			result += input[i];	
		}
		i++;
	}

	return result;

}

void initializeDatastructures(unweightedGraph &G,SpanningTree &T, Biconnected &B, LCAGraph &L)
{
	//tree
	T.parent = new int[G.totalVertices];
	T.levels = new int[G.totalVertices];

	//biconnected
	B.unfinishedTraversal = new int[G.totalVertices];
	B.isBase = new int[G.totalVertices];
	B.isLCA = new int[G.totalVertices];
	B.setOfTreeEdges.resize(G.totalVertices);
	B.lSetOfTreeEdges = new omp_lock_t[G.totalVertices];
	B.cutVertex = new bool[G.totalVertices];
	B.cutVertexAsBridge = new bool[G.totalVertices];

	B.isSpecial = new int*[NUM_THREADS];
	for(int i=0;i<NUM_THREADS;i++)
	{
		B.isSpecial[i] = new int[G.totalVertices];
	}

	B.countTreeEdge = 0;
	
	B.rank = new int*[NUM_THREADS];
	for(int i=0;i<NUM_THREADS;i++)
	{
		B.rank[i] = new int[G.totalVertices];
	}

	B.findLCAisVisited = new int*[NUM_THREADS];
	for(int i=0;i<NUM_THREADS;i++)
	{
		B.findLCAisVisited[i] = new int[G.totalVertices];
	}

	B.findLCACount = new int[NUM_THREADS];
	for(int i=0; i<NUM_THREADS;i++)
		B.findLCACount[i] = 0;

	//LCA Graph

	L.affectedLCA = new int[G.totalVertices]; // size could matter
	L.headAffectedLCA = 0;

	L.partiallyAffectedLCA = new int[G.totalVertices]; //size could give seg fault
	L.headPartiallyAffectedLCA = 0;

	L.partiallyAffectedBaseVertices = new int*[G.totalVertices];
	L.headPartiallyAffectedBaseVertices = new int[G.totalVertices];
	
	L.affectedBaseVertices = new int*[G.totalVertices];
	L.headAffectedBaseVertices = new int[G.totalVertices];
	

	L.globalGraph.adjacencyList1.resize(G.totalVertices);
	L.lLcaGraph = new omp_lock_t[G.totalVertices];

	L.myQueue = new int*[NUM_THREADS];



	//initialization

	for(int i=0; i<NUM_THREADS; i++)
	{
		L.myQueue[i] = new int[G.totalVertices];
	}


	for(int i=0;i<G.totalVertices;i++)
	{	
		T.parent[i] = -1;
		T.levels[i] = -1;
	
		B.unfinishedTraversal[i] = 0;
		B.isBase[i] = 0;
		B.isLCA[i] = 0;	
		omp_init_lock(&(B.lSetOfTreeEdges[i])); 
		B.cutVertex[i] = false;
		B.cutVertexAsBridge[i] = false;
		
		L.partiallyAffectedBaseVertices[i] = new int[1000];
		L.headPartiallyAffectedBaseVertices[i] = 0;

		L.affectedBaseVertices[i] = new int[1000];
		L.headAffectedBaseVertices[i] = 0;

		omp_init_lock(&(L.lLcaGraph[i]));
		
		for(int j = 0; j<NUM_THREADS;j++)
		{
			B.isSpecial[j][i] = -1;
			B.rank[j][i] = -1;
			B.findLCAisVisited[j][i] = -1;
		}
	}

		

}


int main(int argc, char* argv[])
{
	ifstream inputGraph;
	inputGraph.open(argv[1]);

	ifstream inputGraphOG;
	inputGraphOG.open(argv[2]);



	int totalBatchFiles = atoi(argv[3]);

	ofstream output;
	output.open(argv[4+totalBatchFiles]);
	
	unweightedGraph G(inputGraph);
	unweightedGraph OG(inputGraphOG);
	// printf("Read input Graph\n");
	//Read input graph
	
	omp_set_num_threads(NUM_THREADS);
	

	//Spanning Tree
	SpanningTree T;

	//Biconnected datastructure
	Biconnected B;

	//LCA Graph
	LCAGraph L;

	//PostProcessing
	int* postProcessed = new int[G.totalVertices];
	for(int i=0;i<G.totalVertices;i++)
	{
		postProcessed[i] = -1;
	}

	
	initializeDatastructures(G,T,B,L);
	// initialized datastructures
	// printf("initialized datastructures\n");
	

	biconnected(G, T, B, L);

	
	postProcessing(G, T.parent, B.setOfTreeEdges, B.isLCA, B.cutVertex, 
		B.cutVertexAsBridge, output);
	
	output<<endl;	

	// cout<<"Preprocessing done"<<endl;
	
	
	//Actual Program

	for(int currBatchFile = 0; currBatchFile < totalBatchFiles; currBatchFile++)
	{
		string batchSize = getBatchPercent(argv[4+currBatchFile]);
		ifstream queryFile;
		queryFile.open(argv[4+currBatchFile]);
		//Step 1: Reading Query File 
		
		vector<query> queries;
		int totalBatch;
		readQuery(queryFile,queries,totalBatch);
		int count = 0;

	

		L.headAffectedLCA = 0;
		L.headPartiallyAffectedLCA = 0;


		printf("Batch Size = %s\n",&batchSize[0]);
		output<<"Batch Size = "<< batchSize <<endl;

		double totalInsertTime = 0, totalDeleteTime = 0, totalTime = 0,minDelete = DBL_MAX; 
		
		for(int currBatchNumber = 0 ; currBatchNumber < totalBatch ; currBatchNumber++ )
		{
			// cout<<"Current Batch "<<currBatchNumber<<endl;
			// printf("\nBatch %d\n",currBatchNumber);
			int type = queries[count].u;
			int totalQueries = queries[count].v;
			count++;

			
			
			// Updating the degree of every vertex and not 

			for( int i= count ; i< (count+totalQueries) ; i++ )
			{
				int u = queries[i].u;
				int v = queries[i].v;

				if(type == 1)
				{
					G.addEdgeDegree(u,v);
				}
				else
				{
					G.deleteEdgeDegree(u,v);
				}
				
			}

			// auto start = high_resolution_clock::now(); 

			double time1,time2,finalTime = 0;
			time1=omp_get_wtime();

			if(type == 1)
			{
				#pragma omp parallel for schedule(guided)
				for( int i = count ; i< (count+totalQueries) ; i++ )
				{		
					int u = queries[i].u;
					int v = queries[i].v;
					int tid = omp_get_thread_num();
					

					addNonTreeEdge(u,v,G,T,B,L,tid);
									// adds edge and updates the lca graph
				}
				
			}
			else
			{
				deleteEdgeParallel(G, count, totalQueries, queries, T, B, L);
				
			}

			time2 = omp_get_wtime();
			// cout<<"Time For add/delete edges "<<(time2-time1)*1000<<" ms"<<endl;
			finalTime += (time2-time1);
			// cout<<"Batch complete"<<endl;


	        
	        time1 = omp_get_wtime();
				        
			
			isCutVertexAffected(OG, G, T, B, L, postProcessed, currBatchNumber);

			// isCutVertexAll(L.globalGraph, B.unfinishedTraversal, T.parent, B.isLCA, B.cutVertex, G.totalVertices, G.root,
			// 		B.isBase, L.localParent, L.isSafe, L.level);
		
			//Finding end point of bridges are cutvertices or not

			#pragma omp parallel for schedule(guided)
			for(int i=0;i<G.totalVertices;i++)
			{
				B.cutVertexAsBridge[i] = false;
			}

			#pragma omp parallel for schedule(guided)
			for(int i=0;i<G.totalVertices;i++)
			{
				int u = i;
				int v = T.parent[u];

				if(u==v) continue;

				if(B.setOfTreeEdges[u].size()>0)
					continue;
				
				if(G.degree[u] > 1)
					B.cutVertexAsBridge[u] = true;
				if(G.degree[v] > 1)
					B.cutVertexAsBridge[v] = true;
				
			}




			time2=omp_get_wtime();
			// cout<<"Time For connected components "<<(time2-time1)*1000<<" ms"<<endl;
			
			finalTime += (time2-time1);

			postProcessing(G, T.parent, B.setOfTreeEdges, B.isLCA, B.cutVertex,
				B.cutVertexAsBridge,output);
		
			// cout << "Time taken: "<< ((finalTime)*1000) << " milliseconds" << endl;
			
			output << "Time taken: "<< ((finalTime)*1000) << " milliseconds" << endl;
			output<<endl;
			
			totalTime += ((finalTime)*1000);
			count += totalQueries;
			if(type==1)
				totalInsertTime+= ((finalTime)*1000);
			else
			{
				minDelete = min(minDelete,finalTime*1000);
				totalDeleteTime+= ((finalTime)*1000);
			}
		}	
		output<<"Average Insert Time = "<<totalInsertTime/(totalBatch/2)<<endl;
		output<<"Average Delete Time = "<<totalDeleteTime/(totalBatch/2)<<endl;
		output<<"Average Time = "<<totalTime/totalBatch<<endl;

		// cout<<"Average Insert Time = "<<totalInsertTime/(totalBatch/2)<<endl;
		// cout<<"Average Delete Time = "<<totalDeleteTime/(totalBatch/2)<<endl;
		// cout<<"Minimum Delete Time = "<<minDelete<<endl;
		cout<<"Average Time = "<<totalTime/totalBatch<<endl;



	}



	delete[] B.unfinishedTraversal;
	for(int i=0;i<G.totalVertices;i++)
	{
		omp_destroy_lock(&B.lSetOfTreeEdges[i]);
		omp_destroy_lock(&(L.lLcaGraph[i]));
	}
	
}

//   199040