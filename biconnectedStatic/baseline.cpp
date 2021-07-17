#include <omp.h>
#include"unweightedGraph.h"
#include"isCutVertex.h"
#include"bfs.h"
#include"threads.h"
#include<iostream>
#include<set>
#include<algorithm>



using namespace std;



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
		

		totalBridges++;
		
		// if(G.offset[u+1] - G.offset[u] > 1)
		// {
		// 	cutVertices.insert(u);
		// }
		// if(G.offset[v+1] - G.offset[v] > 1)
		// {
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
	for(int i=0;i<bridges.size();i++) 
	{
		output<<bridges[i]<<endl;
	}

	int totalCutVertices = 0;
	output<<"Cut Vertices"<<endl;
	for(int i=0;i<G.totalVertices;i++)
	{
		// if(isLCA[i] && cutVertex[i])
		if(cutVertex[i])
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

double biconnected(unweightedGraph &G, ofstream &output)
{
	//---------------------General-----------------------------------//
	
	//---------------------BiCC Algo Data----------------------------//	
	bool* partOfFundamental = new bool[G.totalVertices];
	bool* unfinishedTraversal = new bool[G.totalVertices];
	
	//---------------------Spanning Tree Data -----------------------//	
	
	int *parent = new int[G.totalVertices];
	int* level = new int[G.totalVertices];

	
	//---------------------isCutvertex Data--------------------------//	

	bool* cutVertex = new bool[G.totalVertices];
	bool* cutVertex2 = new bool[G.totalVertices];
	bool* isLCA = new bool[G.totalVertices];
	bool* isSafe = new bool[G.totalVertices];
	int* componentNumber = new int[G.totalVertices];
	bool* isBase = new bool[G.totalVertices];

	int** baseVertex1 = new int*[NUM_THREADS];
	int** baseVertex2 = new int*[NUM_THREADS];
	int* head = new int[NUM_THREADS];
	
	//---------------------Initialization --------------------------//	
	

	for(int i=0;i<G.totalVertices;i++)
	{ 

		partOfFundamental[i] = false;
		unfinishedTraversal[i] = false;
		isLCA[i] = false;
		level[i] = -1;
		cutVertex[i] = false;
		cutVertex2[i] = false;
		isSafe[i] = false;
		componentNumber[i] = -1;
		isBase[i] = false;
		parent[i] = -1;
	}

	for(int i=0;i<NUM_THREADS;i++)
	{
		head[i] = 0;
		baseVertex1[i] = new int[G.totalEdges];
		baseVertex2[i] = new int[G.totalEdges];
	}

	//-----------------------Algorithm------------------------------//	
	
	//---------------------Step 1 - Create Spanning Tree----------------------------//	
	
	// printf("The root of the tree is %d\n",G.root );
	double time1,time2,totalTime = 0;
	
	time2 = bfs(G, G.root, parent, level);

	totalTime += time2;

	// cout<< "Time taken for spanning tree: "<< ((time2)*1000) << " milliseconds" << endl;
	// cout<<endl;

	// printf("**parent[%d] = %d\n",174535,parent[174535]);
	// printf("**parent[%d] = %d\n",174536,parent[174536]);
	// 
	// cout<<"Spanning Tree Constructed"<<endl;
	

	/*for(int i=0;i<G.totalVertices;i++)
	{
		if(parent[i]==-1)
		{
			cout<<"Problem"<<endl;
		}
	}*/

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

			// if((G.U[i] == 174535 && G.V[i]== 174536) || (G.V[i] == 174535 && G.U[i]== 174536) )
			// {
			// 	printf("************NonTree %d,%d edge is getting added\n",u,v);
			// }
			
			if(parent[u] == v || parent[v] == u) continue; 	//tree edge

			// if(v == 122 || u == 122)
			// {
			// 	printf("Non tree edge %d,%d\n",G.U[i],G.V[i]);
			// }
			if(level[u] < level[v])
			{
				// if(v == 122)
				// {
				// 	printf("Non tree edge %d,%d\n",G.U[i],G.V[i]);
				// }

				partOfFundamental[v] += true;
				unfinishedTraversal[v]+= true;
				v = parent[v];
			}
			if(level[v] < level[u])
			{
				// if(u == 122)
				// {
				// 	printf("Non tree edge %d,%d\n",G.U[i],G.V[i]);
				// }

				partOfFundamental[u] += true;
				unfinishedTraversal[u]+= true;
				u = parent[u];
			}

			while(parent[u] != parent[v])
			{
				// if(v == 122 || u == 122)
				// {
				// 	printf("Non tree edge %d,%d\n",G.U[i],G.V[i]);
				// }
				partOfFundamental[u] += true;
				partOfFundamental[v] += true;
				unfinishedTraversal[u]+= true;
				unfinishedTraversal[v]+= true;
				v = parent[v];
				u = parent[u];
			}
			// if(v == 122 || u == 122)
			// {
			// 	printf("Non tree edge %d,%d\n",G.U[i],G.V[i]);
			// }
			partOfFundamental[u] += true;
			partOfFundamental[v] += true;
			isLCA[parent[u]] += true;

			isBase[u] += true;
			isBase[v] += true;


			// if(parent[u] == 174537)
			// {
			// 	printf("*** base edge %d,%d added for non tree edge %d,%d\n",u,v,G.U[i],G.V[i]);
			// }

			// baseVertex1[tid][head[tid]] = u;
			// baseVertex2[tid][head[tid]] = v;
			// head[tid]++;
			
			b1[*h] = u;
			b2[*h] = v;
			(*h)++;

			// if((G.U[i] == 174535 && G.V[i]== 174536) || (G.U[i] == 174535 && G.V[i]== 174536))
			// {
			// 	printf("**********LCA of NonTree Edge %d,%d is %d \n",G.U[i],G.V[i],parent[u]);
			// }

		}
	}


	// printf("Parent of [%d] is %d\n",122,parent[122]);
	// printf("partOfFundamental[122] = %d\n",partOfFundamental[122]);
		

	// for(int i=190429;i<=190435;i++)
	// {
	// 	printf("unfinishedTraversal[%d] = %d\n",i,unfinishedTraversal[i]);
	// }

	time2=omp_get_wtime();
	
	
	totalTime += (time2-time1);
	//-------------Step 3 - Calculate whether LCA is a cutvertex--------------------//	

	// printf("isCutVertex[%d] == %d\n",174537,cutVertex[174537]);
	
	time2 = markCutVertices(baseVertex1,baseVertex2,head, unfinishedTraversal, isLCA, isBase,parent, cutVertex, isSafe, componentNumber, G.totalVertices, G.root);
	
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

	// printf("isCutVertex[%d] == %d\n",174537,cutVertex[174537]);
	
	//--------------------------------Completed--------------------------------//	
	
	
	
	
	postProcessing(G, parent, partOfFundamental, isLCA, cutVertex,output);
	output << "Time taken: "<< ((totalTime)*1000) << " milliseconds" << endl;
	output<<endl;
	

	delete[] partOfFundamental;
	delete[] unfinishedTraversal;
	
	
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

	return result;

}

int main(int argc, char* argv[])
{
	ifstream inputGraph;
	inputGraph.open(argv[1]);

	int totalBatchFiles = atoi(argv[2]);

	ofstream output;
	output.open(argv[3+totalBatchFiles]);

	unweightedGraph G(inputGraph);
	// cout<<"Created Graph"<<endl;
	// //Read input graph

	int totalThreads = NUM_THREADS;
	omp_set_num_threads(totalThreads);



	biconnected(G,output);
	
	// cout<<"Computed on original graph"<<endl;
	
	//Actual Program

	
	//Step 1: Reading Query File 
	for(int currBatchFile = 0; currBatchFile < totalBatchFiles; currBatchFile++)
	{
		string batchSize = getBatchPercent(argv[3+currBatchFile]);
		ifstream queryFile;
		queryFile.open(argv[3+currBatchFile]);
	
		vector<query> queries;
		int totalBatch;
		readQuery(queryFile,queries,totalBatch);
		int count = 0;

		printf("Batch Size = %s\n",&batchSize[0]);
		output<<"Batch Size = "<< batchSize <<endl;
	
		
		// cout<<"Starting Batch"<<endl;

		double totalInsertTime = 0, totalDeleteTime = 0, totalTime = 0; 
		for(int x=0 ; x < totalBatch ; x++ )
		{
			// cout<<"Current Batch "<<x<<endl;
			int type = queries[count].u;
			int totalQueries = queries[count].v;
			count++;

			
			// G.print();

			for( int i= count ; i< (count+totalQueries) ; i++ )
			{
				int u = queries[i].u;
				int v = queries[i].v;

				if(type == 1)
				{
					G.addEdgeSet(u,v);
				}
				else
				{
					G.deleteEdgeSet(u,v);
				}

			}
			G.recreate();

			// cout<<"G updated "<<endl;

			double time;
			time = biconnected(G,output);
			totalTime += time;
			if(type==1)
				totalInsertTime += time;
			else
				totalDeleteTime += time;

			count += totalQueries;
		}	
		
		output<<"Average Time = "<<totalTime/totalBatch<<endl;
		cout<<"Average Time = "<<totalTime/totalBatch<<endl;
	}
}