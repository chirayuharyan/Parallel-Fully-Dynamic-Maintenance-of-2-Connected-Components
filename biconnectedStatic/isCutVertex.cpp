#include"isCutVertex.h"
#include<queue>
#include<iostream>
#include<unordered_map>
#include<omp.h>
#include<math.h>
using namespace std;




void mergeVector(vector<vector<pair<int,int>>> &vecIn, vector<pair<int,int>> &vecOut)
{
	int totalSize = 0;
	int *starting = new int[vecIn.size()+1];
	starting[0] = 0;
	for(int i=0; i<vecIn.size();i++)
	{
		totalSize += vecIn[i].size();
		starting[i+1] = totalSize;
	}
	vecOut.resize(totalSize);
	
	#pragma omp parallel
	{
		int i = omp_get_thread_num();
		for(int j=0;j<vecIn[i].size();j++)
		{
			vecOut[starting[i]+j] = vecIn[i][j];
		}
	}

}



double markCutVerticesOld(int** baseVertex1, int** baseVertex2, int* head, bool* isSafe, bool *isLCA, bool* isBase, int* parent, bool *isCutVertex, bool* isSafeOld, int* componentNumber, int n, int root)
{
	


	//-----------------------------Initialization--------------------------------------//
	
	// vector<pair<int,int>> edgeListCombined;
	// mergeVector(edgeList,edgeListCombined);
	int *totalChild = new int[NUM_THREADS];
	double time1 = 0, time2 = 0;
	
	for(int i=0; i<NUM_THREADS;i++)
		totalChild[i] = 0;

	int start, end, block1,tid, i,block2;
	block1 = ceil(double(n)/NUM_THREADS);
	// block2 = ceil(double(edgeListCombined.size())/NUM_THREADS);
	int iteration = 0;
	bool flag = true;


	time1=omp_get_wtime();

	
	

	//-----------------------------Algorithm--------------------------------------//
	
	
	
	
	// printf("Starting \n");
	
	#pragma omp parallel private(tid,i)
	{
		tid = omp_get_thread_num();
		int start1 = tid * block1;
		int end1 = min(start1 + block1 -1, n-1);

		int start2 = 0;
		int end2 = head[tid]-1;

		int u,v,comU,comV,parentCom;
		bool temp;

		for(i=start1;i<=end1;i++)
		{
			if(isBase[i])
			{
				componentNumber[i] = i;
				// isSafe[i] = unfinishedTraversal[i];
			}
			// isBase[i]?componentNumber[i] = i,isSafe[i] = unfinishedTraversal[i]:NULL;
		}
		
		
		while(flag)
		{
			#pragma omp barrier //let all the threads come inside while
			#pragma omp single 
			{
				// printf("Step 1 completed %d\n",tid);
				flag = false;
				iteration++;
			}

			for(i = start2; i<=end2 ;i++)
			{	
				u = baseVertex1[tid][i];
				v = baseVertex2[tid][i];	
				comU = componentNumber[u];
				comV = componentNumber[v];
				if(comU == comV) continue;
				flag += true;
				
				if(iteration%2)
				{
					if(comU < comV) //min
						componentNumber[comV] = comU;
					else
						componentNumber[comU] = comV;
				}
				else
				{
					if(comU > comV) //max
						componentNumber[comV] = comU;
					else
						componentNumber[comU] = comV;
				}
				
				// temp = isSafe[comU]  + isSafe[comV]; 		
				// if(temp)
				// 	isSafe[parentCom] += temp; // check for bugs, componentNumber[comU] might get changed by another thread

			}

			#pragma omp barrier
			
			for(i=start1; i<=end1; i++)
			{
				if(componentNumber[i] == -1) continue;
				while( componentNumber[i] != componentNumber[componentNumber[i]])
					componentNumber[i] = componentNumber[componentNumber[i]];
			}

			/*#pragma omp barrier
			for(i = start2; i<=end2 ;i++)
			{
				edgeListCombined[i].first = componentNumber[edgeListCombined[i].first];
				edgeListCombined[i].second = componentNumber[edgeListCombined[i].first];
			}*/
		}

		#pragma omp barrier

		for(i=start1;i<=end1;i++)
		{
			if(isBase[i] && isSafe[i])
			{
				isSafe[componentNumber[i]] += true;
			}
		}

		#pragma omp barrier

		for(i=start1;i<=end1;i++)
		{
			if(componentNumber[i]==i && !isSafe[i])
			{		

				isCutVertex[parent[i]] += true; 

				if(parent[i]==root)
					totalChild[tid]++;
			}
		}

	}
	
	// printf("Total Iteration = %d\n", iteration);

	

	int total = 0;
	for(int i=0;i<NUM_THREADS;i++)
		total+=totalChild[i];

	isCutVertex[root] = (total > 1);

	time2=omp_get_wtime();


	return (time2-time1);

}	

double markCutVertices(int** baseVertex1, int** baseVertex2, int* head, bool* isSafe, bool *isLCA, bool* isBase, int* parent, bool *isCutVertex, bool* isSafeOld, int* componentNumber, int n, int root)
{
	


	//-----------------------------Initialization--------------------------------------//
	
	// vector<pair<int,int>> edgeListCombined;
	// mergeVector(edgeList,edgeListCombined);
	int *totalChild = new int[NUM_THREADS];
	double time1 = 0, time2 = 0;
	
	for(int i=0; i<NUM_THREADS;i++)
		totalChild[i] = 0;

	int start, end, block1,tid, i,block2;
	block1 = ceil(double(n)/NUM_THREADS);
	// block2 = ceil(double(edgeListCombined.size())/NUM_THREADS);
	int iteration = 0;
	bool flag = true;


	time1=omp_get_wtime();

	
	

	//-----------------------------Algorithm--------------------------------------//
	
	
	
	
	// printf("Starting \n");
	
	#pragma omp parallel private(tid,i)
	{
		tid = omp_get_thread_num();
		int start1 = tid * block1;
		int end1 = min(start1 + block1 -1, n-1);

		int start2 = 0;
		int end2 = head[tid]-1;

		int u,v,comU,comV,parentCom;
		bool temp;

		for(i=start1;i<=end1;i++)
		{
			if(isBase[i])
			{
				componentNumber[i] = i;
				// isSafe[i] = unfinishedTraversal[i];
			}
			// isBase[i]?componentNumber[i] = i,isSafe[i] = unfinishedTraversal[i]:NULL;
		}
		
		
		while(flag)
		{
			#pragma omp barrier //let all the threads come inside while
			#pragma omp single 
			{
				// printf("Step 1 completed %d\n",tid);
				flag = false;
				iteration++;
			}

			for(i = start2; i<=end2 ;i++)
			{	
				u = baseVertex1[tid][i];
				v = baseVertex2[tid][i];	
				comU = componentNumber[u];
				comV = componentNumber[v];
				if(comU == comV) continue;
				flag += true;
				
				if(iteration%2)
				{
					if(comU < comV) //min
						componentNumber[comV] = comU;
					else
						componentNumber[comU] = comV;
				}
				else
				{
					if(comU > comV) //max
						componentNumber[comV] = comU;
					else
						componentNumber[comU] = comV;
				}
				
				// temp = isSafe[comU]  + isSafe[comV]; 		
				// if(temp)
				// 	isSafe[parentCom] += temp; // check for bugs, componentNumber[comU] might get changed by another thread

			}

			#pragma omp barrier
			
			for(i=start1; i<=end1; i++)
			{
				if(componentNumber[i] == -1) continue;
				while( componentNumber[i] != componentNumber[componentNumber[i]])
					componentNumber[i] = componentNumber[componentNumber[i]];
			}

			int newHead = 0;
			#pragma omp barrier
			for(i = start2; i<=end2 ;i++)
			{
				u = baseVertex1[tid][i];
				v = baseVertex2[tid][i];
				if(componentNumber[u] != componentNumber[v])
				{
					baseVertex1[tid][newHead] = componentNumber[u];
					baseVertex2[tid][newHead] = componentNumber[v];
					newHead++;
				}
			}

			end2 = newHead - 1;
		}

		#pragma omp barrier

		for(i=start1;i<=end1;i++)
		{
			// if(componentNumber[i] == 174535)
			// 	printf("%d is in component of 174535 with isSafe[%d] = %d\n",i,i,isSafe[i]);

			// if(isBase[i] && parent[i] == 174537)
			// {
			// 	printf("%d is base vertex with unfinishedTraversal[%d] = %d with componentNumber[%d] = %d\n",i,i,isSafe[i],i,componentNumber[i]);
			// }
			if(isBase[i] && isSafe[i])
			{
				isSafe[componentNumber[i]] += true;
			}
		}



		#pragma omp barrier

		for(i=start1;i<=end1;i++)
		{
			if(componentNumber[i]==i && !isSafe[i])
			{		
				// if(parent[i] == 174537)
				// 	printf("Marked as cutvertex by edge %d,%d \n", i,parent[i]);
				isCutVertex[parent[i]] += true; 

				if(parent[i]==root)
					totalChild[tid]++;
			}
		}

	}
	
	// printf("Total Iteration = %d\n", iteration);

	

	int total = 0;
	for(int i=0;i<NUM_THREADS;i++)
		total+=totalChild[i];

	isCutVertex[root] = (total > 1);

	time2=omp_get_wtime();

	int totalComponents = 0;
	for(int i=0;i<n; i++)
	{
		if(componentNumber[i]==i)
			totalComponents++;

	}
	// std::cout<<"Total components = "<<totalComponents<<endl;

	return (time2-time1);

}	




