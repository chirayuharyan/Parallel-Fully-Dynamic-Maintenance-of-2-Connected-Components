#include <omp.h>
#include"unweightedGraph.h"
#include"isCutVertex.h"
#include"bfs.h"
#include"threads.h"
#include<iostream>
#include<set>
#include<algorithm>

using namespace std;

typedef struct SpanningTree
{
	int *parent, *level,*offset,*children;
}SpanningTree;

typedef struct LCAGraph
{
	int **baseVertex1, **baseVertex2, *head, *componentNumber;
	bool *isSafe;
}LCAGraph;

typedef struct Biconnected
{
	bool *partOfFundamental, *unfinishedTraversal, *isLCA, *isBase, *cutVertex;
	int **affectedLCA, *head;
}Biconnected;

void postProcessing(unweightedGraph &G, int* parent, bool* partOfFundamental, bool* isLCA, bool* cutVertex, ofstream &output)
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

		if(partOfFundamental[u])
			continue;

		// if((u == 174535 && v == 174537)|| (v == 174535 && u == 174537))
		// 	printf("%d,%d is a bridge \n",u,v);
		
		totalBridges++;
		
		// if(G.degree[u] > 1)
		// {
		// 	// if(u == 174537)
		// 	// 		printf("%d is a cutvertex because edge %d,%d is a bridge\n",u,u,v);
		// 	cutVertices.insert(u);
		// }
		// if(G.degree[v] > 1)
		// {
		// 	// if(v == 174537)
		// 	// 		printf("%d is a cutvertex because edge %d,%d is a bridge\n",v,u,v);
		// 	cutVertices.insert(v);
		// }

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
	for(unsigned int i=0;i<bridges.size();i++) 
	{
		output<<bridges[i]<<endl;
	}

	int totalCutVertices = 0;
	output<<"Cut Vertices"<<endl;
	for(int i=0;i<G.totalVertices;i++)
	{
		// if(isLCA[i] && cutVertex[i])
		if( cutVertex[i])
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
	queryFile.close();
}


void queryToArray(vector<query> &queries, int*& queriesU, int*& queriesV, int &totalQueries)
{
	totalQueries = queries.size();
	queriesU = new int[totalQueries];
	queriesV = new int[totalQueries];
	for(int i=0;i<totalQueries;i++)
	{
		queriesU[i] = queries[i].u;
		queriesV[i] = queries[i].v;
	}
}

double biconnected(unweightedGraph &G, ofstream &output, Biconnected *B, SpanningTree *T, LCAGraph *H)
{
	bool* partOfFundamental = B->partOfFundamental;
	bool* unfinishedTraversal = B->unfinishedTraversal;
	bool* isLCA = B->isLCA;
	bool* isBase = B->isBase;
	bool* cutVertex = B->cutVertex ;

	int* parent = T->parent;
	int* level = T->level;
	int* offset = T->offset;
	int* children = T->children;

	int* componentNumber = H->componentNumber;
	int** baseVertex1 = H->baseVertex1;
	int** baseVertex2 = H->baseVertex2;
	int* head = H->head;
	bool* isSafe = H->isSafe;

	//-----------------------Algorithm------------------------------//	
	
	//---------------------Step 1 - Create Spanning Tree----------------------------//	
	

	double time1,time2,totalTime = 0;
	
	time2 = bfs(G, G.root, parent, level);

	totalTime += time2;

	// cout<< "Time taken for spanning tree: "<< ((time2)*1000) << " milliseconds" << endl;
	// cout<<endl;
	
	//------------------Step 1.2 - Creating children array for SpanningTree----------------------------//	

	for(int i=0;i<G.totalVertices;i++)
	{
		if(i == G.root) continue;
		offset[parent[i]]++;
	}

	int curr = 0, next = 0;
	int *currPointer = new int[G.totalVertices]; 
	for(int i=0;i<G.totalVertices;i++)
	{
		next += offset[i];
		offset[i] = curr;
		currPointer[i] = curr;
		curr = next;
	}

	for(int i=0;i<G.totalVertices;i++)
	{
		int p_i = parent[i];
		int index = currPointer[p_i];
		currPointer[p_i]++;
		children[index] = i;
	}


	//---------------------Step 2 - Add Non-tree Edge--------------------------------//	
	
	time1=omp_get_wtime();

	#pragma omp parallel 
	{
		int tid = omp_get_thread_num();
		int *b1 = baseVertex1[tid];
		int *b2 = baseVertex2[tid];
		int *h = &head[tid];
		#pragma omp for schedule(guided)
		for(int i=0;i<G.totalEdges/2;i++)
		{
			int u = G.U[i];
			int v = G.V[i];

			if(parent[u] == v || parent[v] == u) continue; 	//tree edge

			
			if(level[u] < level[v])
			{
				partOfFundamental[v] += true;
				unfinishedTraversal[v]+= true;
				v = parent[v];
			}
			if(level[v] < level[u])
			{
				partOfFundamental[u] += true;
				unfinishedTraversal[u]+= true;
				u = parent[u];
			}

			while(parent[u] != parent[v])
			{
				partOfFundamental[u] += true;
				partOfFundamental[v] += true;
				unfinishedTraversal[u]+= true;
				unfinishedTraversal[v]+= true;
				v = parent[v];
				u = parent[u];
			}
			partOfFundamental[u] += true;
			partOfFundamental[v] += true;
			isLCA[parent[u]] += true;

			isBase[u] += true;
			isBase[v] += true;

			// baseVertex1[tid][head[tid]] = u;
			// baseVertex2[tid][head[tid]] = v;
			// head[tid]++;
			
			b1[*h] = u;
			b2[*h] = v;
			(*h)++;

		}
	}


	time2=omp_get_wtime();
	
	// cout<< "Time taken for adding non-tree edges: "<< ((time2-time1)*1000) << " milliseconds" << endl;
	// cout<<endl;
	
	totalTime += (time2-time1);
	//-------------Step 3 - Calculate whether LCA is a cutvertex--------------------//	

	
	#pragma omp parallel for
	for(int i=0;i<G.totalVertices;i++)
	{
		if(isBase[i])
			componentNumber[i] = i;
	}

	time2 = markCutVertices(baseVertex1,baseVertex2, head, unfinishedTraversal, isLCA, isBase, parent, cutVertex, componentNumber, G.totalVertices, G.root, isSafe);
	
	totalTime += time2;

	time1=omp_get_wtime();	
	#pragma omp parallel for schedule(guided)
	for(int i=0;i<G.totalVertices;i++)
	{
		int u = i;
		int v = parent[u];

		if(u==v) continue;

		if(partOfFundamental[u])
			continue;
		
		if(G.offset[u+1] - G.offset[u] > 1)
		{
			cutVertex[u] += true;
		}
		if(G.offset[v+1] - G.offset[v] > 1)
		{
			cutVertex[v] += true;
		}
	}
	
	time2=omp_get_wtime();
	totalTime += (time2-time1);
	
	// cout<< "Time taken for calculating whether it is a cutvertex: "<< ((time2)*1000) << " milliseconds" << endl;
	// cout<<endl;
	
	//--------------------------------Completed--------------------------------//	
	
	
	
	
	// postProcessing(G, parent, partOfFundamental, isLCA, cutVertex,output);
	output << "Time taken: "<< ((totalTime)*1000) << " milliseconds" << endl;
	output<<endl;
	// cout<< "Time taken: "<< ((totalTime)*1000) << " milliseconds" << endl;
	// cout<<endl;
	
	// printf("unfinishedTraversal Values\n");
	// for(int i=0;i<G.totalVertices;i++)
	// {
	// 	printf("%d ",unfinishedTraversal[i]);
	// }
	// printf("\n");
	
	return ((totalTime)*1000);
}

double addBatch(unweightedGraph &G, int* U,int* V, int totalQueries, ofstream &output, Biconnected *B, SpanningTree *T, LCAGraph *H)
{
	bool* partOfFundamental = B->partOfFundamental;
	bool* unfinishedTraversal = B->unfinishedTraversal;
	bool* isLCA = B->isLCA;
	bool* isBase = B->isBase;
	bool* cutVertex = B->cutVertex;
	int** affectedLCA = B->affectedLCA;
	int* headAffectedLCA = B->head;

	int* parent = T->parent;
	int* level = T->level;
	int* offset = T->offset;
	int* children = T->children;

	int* componentNumber = H->componentNumber;
	int** baseVertex1 = H->baseVertex1;
	int** baseVertex2 = H->baseVertex2;
	int* head = H->head;
	bool* isSafe = H->isSafe;



	//-----------------------Algorithm------------------------------//	
	
	double time1,time2,totalTime = 0;
	//---------------------Step 1 - Add Non-tree Edge--------------------------------//	
	
	time1=omp_get_wtime();

	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		int *b1 = baseVertex1[tid];
		int *b2 = baseVertex2[tid];
		int *h = &head[tid];
		// int *affected = affectedLCA[tid];
		// int *hAffectedLCA = &headAffectedLCA[tid];
		// *hAffectedLCA = 0;

		#pragma omp for schedule(guided)
		for(int i=0;i<totalQueries;i++)
		{
			int u = U[i];
			int v = V[i];



			if(level[u] < level[v])
			{
				while(level[u] < level[parent[v]])
				{

					// if(isLCA[v])
					// {
					// 	affected[*hAffectedLCA] = v;
					// 	(*hAffectedLCA)++;
					// }
					partOfFundamental[v] += true;
					unfinishedTraversal[v]+= true;
					v = parent[v];
				}
				if(parent[v] == u) // back edges
				{
					isLCA[u] += true;
					partOfFundamental[v] += true;
					isBase[v] += true;
					continue;
				}
				else
				{
					partOfFundamental[v] += true;
					unfinishedTraversal[v]+= true;
					v = parent[v];
				}
			}

			if(level[v] < level[u])
			{
				while(level[v] < level[parent[u]])
				{
					// if(isLCA[u])
					// {
					// 	affected[*hAffectedLCA] = u;
					// 	(*hAffectedLCA)++;
					
					partOfFundamental[u] += true;
					unfinishedTraversal[u]+= true;
					u = parent[u];
				}

				if(parent[u] == v) // back edges
				{
					isLCA[v] += true;
					partOfFundamental[u] += true;
					isBase[u] += true;
					continue;
				}
				else
				{
					partOfFundamental[u] += true;
					unfinishedTraversal[u]+= true;
					u = parent[u];
				}	
			}

			
			while(parent[u] != parent[v])
			{
				// if(isLCA[u])
				// {
				// 	affected[*hAffectedLCA] = u;
				// 	(*hAffectedLCA)++;
				// }
				// if(isLCA[v])
				// {
				// 	affected[*hAffectedLCA] = v;
				// 	(*hAffectedLCA)++;
				// }

				partOfFundamental[u] += true;
				partOfFundamental[v] += true;
				unfinishedTraversal[u]+= true;
				unfinishedTraversal[v]+= true;
				v = parent[v];
				u = parent[u];
			}

			partOfFundamental[u] += true;
			partOfFundamental[v] += true;
			isLCA[parent[u]] += true;
			// affected[*hAffectedLCA] = parent[u];
			// (*hAffectedLCA)++;

			isBase[u] += true;
			isBase[v] += true;
			
			b1[*h] = u;
			b2[*h] = v;
			(*h)++;


		}
	}



	time2=omp_get_wtime();
	
	// cout<< "Time taken for adding non-tree edges: "<< ((time2-time1)*1000) << " milliseconds" << endl;
	// cout<<endl;
	


	totalTime += (time2-time1);
	//-------------Step 3 - Calculate whether LCA is a cutvertex--------------------//	


	time2 = markCutVerticesAffected(baseVertex1,baseVertex2, head, affectedLCA, headAffectedLCA, unfinishedTraversal, isLCA,
	isBase, parent, offset, children , cutVertex, componentNumber, G.totalVertices, G.root, isSafe);
	
	


	totalTime += time2;
	

	time1=omp_get_wtime();	
	#pragma omp parallel for schedule(guided)
	for(int i=0;i<G.totalVertices;i++)
	{
		int u = i;
		int v = parent[u];

		if(u==v) continue;

		if(partOfFundamental[u])
			continue;
		
		if(G.offset[u+1] - G.offset[u] > 1)
		{
			cutVertex[u] += true;
		}
		if(G.offset[v+1] - G.offset[v] > 1)
		{
			cutVertex[v] += true;
		}
	}
	
	time2=omp_get_wtime();
	totalTime += (time2-time1);
	// cout<< "Time taken for calculating whether it is a cutvertex: "<< ((time2)*1000) << " milliseconds" << endl;
	// cout<<endl;


	
	//--------------------------------Completed--------------------------------//	
	
	
	
	
	// postProcessing(G, parent, partOfFundamental, isLCA, cutVertex,output);
	output << "Time taken: "<< ((totalTime)*1000) << " milliseconds" << endl;
	output<<endl;
	// cout<< "Time taken: "<< ((totalTime)*1000) << " milliseconds" << endl;
	// cout<<endl;
	
	
	return ((totalTime)*1000);
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

	if(input[i] == 'k')
		result += 'k';

	return result;

}


int main(int argc, char* argv[])
{
	ifstream inputGraph;
	inputGraph.open(argv[1]);

	int totalBatchFiles = atoi(argv[2]);

	ofstream output;
	output.open(argv[3+totalBatchFiles]);


	// cout<<"Started"<<endl;
	unweightedGraph G(inputGraph);
	cout<<"Created Graph"<<endl;
	//Read input graph

	int totalThreads = NUM_THREADS;
	omp_set_num_threads(totalThreads);
	// printf("Total Threads = %d\n",totalThreads);



//-------------------------------------------DATASTRUCTES -------------------------------------------------//

	//---------------------General-----------------------------------//
	
	//---------------------BiCC Algo Data----------------------------//	
	Biconnected B;
	B.partOfFundamental = new bool[G.totalVertices];
	B.unfinishedTraversal = new bool[G.totalVertices];
	B.isLCA = new bool[G.totalVertices];
	B.isBase = new bool[G.totalVertices];
	B.cutVertex = new bool[G.totalVertices];
	B.affectedLCA = new int*[NUM_THREADS];
	B.head = new int[NUM_THREADS];
	
	//---------------------Spanning Tree Data -----------------------//	
	
	SpanningTree T;
	T.parent = new int[G.totalVertices];
	T.level = new int[G.totalVertices];
	T.offset = new int[G.totalVertices+1];
	T.children = new int[G.totalVertices];

	
	//---------------------isCutvertex Data--------------------------//	

	LCAGraph H;
	
	H.componentNumber = new int[G.totalVertices];
	H.baseVertex1 = new int*[NUM_THREADS];
	H.baseVertex2 = new int*[NUM_THREADS];
	H.head = new int[NUM_THREADS];
	H.isSafe = new bool[G.totalVertices];
	
	//---------------------Initialization --------------------------//	
	

	for(int i=0;i<G.totalVertices;i++)
	{ 
		B.partOfFundamental[i] = false;
		B.unfinishedTraversal[i] = false;
		B.isLCA[i] = false;
		B.isBase[i] = false;
		B.cutVertex[i] = false;

		T.parent[i] = -1;
		T.level[i] = -1;
		T.offset[i] = 0;
		H.componentNumber[i] = -1;
	}
	T.offset[G.totalVertices] = G.totalVertices-1;

	for(int i=0;i<NUM_THREADS;i++)
	{
		B.head[i] = 0;
		B.affectedLCA[i] = new int[G.totalVertices];

		H.head[i] = 0;
		H.baseVertex1[i] = new int[G.totalEdges];
		H.baseVertex2[i] = new int[G.totalEdges];
	}



//--------------------------------------------Preprocessing-------------------------------------------------//

	
	biconnected(G,output,&B,&T,&H);
	
	// cout<<"Computed on original graph"<<endl;
	
	
	
//--------------------------------------------------Program-------------------------------------------------//

	//Step 1: Reading Query File 
	
	// string queryPrefix(argv[3]);
	for(int currBatchFile = 0; currBatchFile < totalBatchFiles; currBatchFile++)
	{
		string batchSize = getBatchPercent(argv[2+currBatchFile+1]);
		ifstream queryFile;
		queryFile.open(argv[2+currBatchFile+1]);
		vector<query> queries;
		int *queriesU, *queriesV, totalQueries;
		int totalBatch;
		readQuery(queryFile,queries,totalBatch);
		queryToArray(queries,queriesU,queriesV,totalQueries);
		int count = 0;

		printf("Batch Size = %s\n",&batchSize[0]);
		output<<"Batch Size = "<< batchSize <<endl;
		
		double totalTime = 0; 
		for(int x=0 ; x < totalBatch ; x++ )
		{
			// cout<<"Current Batch "<<x<<endl;
			int totalQueries = queries[count].v;
			count++;


			for( int i= count ; i< (count+totalQueries) ; i++ )
			{
				int u = queriesU[i];
				int v = queriesV[i];

				G.degree[u]++;
				G.degree[v]++;

			}

			// cout<<"G updated "<<endl;

			double time;
			time = addBatch(G,queriesU+count,queriesV+count,totalQueries,output,&B,&T,&H);

			totalTime += time;
			count += totalQueries;
		}	
		output<<"Average Time = "<<totalTime/totalBatch<<endl;
		cout<<"Average Time = "<<totalTime/totalBatch<<endl;
	}
	
}