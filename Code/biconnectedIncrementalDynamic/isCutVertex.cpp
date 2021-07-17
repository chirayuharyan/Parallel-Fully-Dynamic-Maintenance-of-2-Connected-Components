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
	for(unsigned int i=0; i<vecIn.size();i++)
	{
		totalSize += vecIn[i].size();
		starting[i+1] = totalSize;
	}
	vecOut.resize(totalSize);
	
	#pragma omp parallel
	{
		int i = omp_get_thread_num();
		for(unsigned int j=0;j<vecIn[i].size();j++)
		{
			vecOut[starting[i]+j] = vecIn[i][j];
		}
	}

}


double markCutVertices(int** baseVertex1, int** baseVertex2, int* head, bool* unfinishedTraversal, bool *isLCA, bool* isBase, int* parent, bool *isCutVertex, int* componentNumber, int n, int root, bool* isSafe)
{
	


	//-----------------------------Initialization--------------------------------------//
	
	// vector<pair<int,int>> edgeListCombined;
	// mergeVector(edgeList,edgeListCombined);
	int *totalChild = new int[NUM_THREADS];
	// bool* isSafe = new bool[n];

	
	double time1 = 0, time2 = 0;
	
	for(int i=0; i<NUM_THREADS;i++)
		totalChild[i] = 0;

	int  block1,tid, i;
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

		int u,v,comU,comV;
		

		for(i=start1;i<=end1;i++)
		{
			if(isBase[i])
			{
				componentNumber[i] = i;
				isSafe[i] = unfinishedTraversal[i];
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
				// printf("Iteration = %d\n",iteration);
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
				{
					componentNumber[i] = componentNumber[componentNumber[i]];
					// if(i == 239243)
					// 	printf("componentNumber[%d] = %d\n",i,componentNumber[i]);
				}
			}

			// #pragma omp single 
			// {
			// 	printf("shortcutting completed\n");
			// }

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
			head[tid] = newHead;

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

	int totalComponents = 0;
	for(int i=0;i<n; i++)
	{
		if(componentNumber[i]==i)
			totalComponents++;

	}
	// std::cout<<"Total components = "<<totalComponents<<endl;

	return (time2-time1);

}	

double markCutVerticesAffected(int** baseVertex1, int** baseVertex2, int* head, int** affectedLCA, 
	int* headAffectedLCA, bool* unfinishedTraversal, bool *isLCA, bool* isBase, int* parent,int* offset, 
	int* children, bool *isCutVertex, int* componentNumber, int n, int root, bool* isSafe)
{
	

	//-----------------------------Initialization--------------------------------------//
	
	int *totalChild = new int[NUM_THREADS];

	// bool *affectedLCAFlag = new bool[n];
	// for(int i=0;i<n;i++)
	// 	affectedLCAFlag[i] = false; 
	// bool* isSafe = new bool[n];
	
	double time1 = 0, time2 = 0;
	
	for(int i=0; i<NUM_THREADS;i++)
		totalChild[i] = 0;

	int block1,tid, i;
	block1 = ceil(double(n)/NUM_THREADS);
	// block2 = ceil(double(edgeListCombined.size())/NUM_THREADS);
	int iteration = 0;
	bool flag = true;


	time1=omp_get_wtime();

	
	// #pragma omp parallel for
	// for(int i=0;i<n;i++)
	// {
	// 	// isSafe[i] = false;
	// 	isSafe[i] = unfinishedTraversal[i];
	// }

	// cout<<"baseVertex1[0] = "<<baseVertex1[0]<<endl;
	// cout<<"baseVertex1[0][0] = "<<baseVertex1[0][0]<<endl;
	// cout<<"affectedLCA[0] = "<<affectedLCA[0]<<endl;
	// cout<<"affectedLCA[0][0] = "<<affectedLCA[0][0]<<endl;
	// printf("headAffectedLCA[0] = %d\n",headAffectedLCA[0]);
	// for(int i=0;i<headAffectedLCA[0];i++)
	// {
	// 	printf("%d ",affectedLCA[0][i]);
	// }	
	// printf("\n");
	

	//-----------------------------Algorithm--------------------------------------//
	
	
	
	
	// printf("Starting \n");
	
	#pragma omp parallel private(tid,i)
	{
		tid = omp_get_thread_num();
		int start1 = tid * block1;
		int end1 = min(start1 + block1 -1, n-1);

		int start2 = 0;
		int end2 = head[tid]-1;

		int u,v,comU,comV;
		

		// int *affected = affectedLCA[tid];
		// int *hAffectedLCA = &headAffectedLCA[tid];

		// int totalAffectedLCA = headAffectedLCA[tid];


		for(i = start1 ; i<=end1; i++)
		{
			isCutVertex[i] = false;
		}
		for(i = start1 ; i<=end1; i++)
		{
			isSafe[i] = unfinishedTraversal[i];
		}
		for(i = start1 ; i<=end1; i++)
		{
			if(isBase[i] && componentNumber[i] == -1)
			{
				componentNumber[i] = i;
			}
		}

		// for(i=0; i<totalAffectedLCA;i++)
		// {
		// 	isCutVertex[affected[i]] = false; 
		// 	if(affected[i]<0 ||affected[i]>=n)
		// 		printf("LCA greater than totalVertices %d\n",affected[i]);
		// 	affectedLCAFlag[affected[i]] += true;
		// }

		// for(i=0; i<totalAffectedLCA;i++)
		// {
		// 	lca = affected[i];
		// 	// printf("Affected LCA = %d\n",lca);
		// 	endG = offset[lca+1];
		// 	for(u = offset[lca]; u < endG; u++)
		// 	{
		// 		currChild = children[u];
		// 		// printf("Parent = %d and currChild = %d\n",lca,currChild);
		// 		// if(!isBase[currChild]) continue;
		// 		isSafe[currChild] = unfinishedTraversal[currChild];
		// 		if(componentNumber[i] == -1)
		// 			componentNumber[i] = i;
		// 	}
		// }

		// for(i = start2; i<=end2 ;i++)
		// {	
		// 	u = baseVertex1[tid][i];
		// 	v = baseVertex2[tid][i];	
		// 	isSafe[u] += unfinishedTraversal[u];
		// 	isSafe[v] += unfinishedTraversal[v];
		// 	if(componentNumber[u] == -1)
		// 		componentNumber[u] = u;
		// 	if(componentNumber[v] == -1)
		// 		componentNumber[v] = v;
		// }



		// for(i=start1;i<=end1;i++)
		// {
		// 	if(isBase[i])
		// 	{
		// 		componentNumber[i] = i;
		// 		// isSafe[i] = unfinishedTraversal[i];
		// 	}
		// 	// isBase[i]?componentNumber[i] = i,isSafe[i] = unfinishedTraversal[i]:NULL;
		// }
		
		
		while(flag)
		{
			#pragma omp barrier //let all the threads come inside while
			#pragma omp single 
			{
				// printf("Step 1 completed %d\n",tid);
				flag = false;
				iteration++;
				// printf("Iteration = %d\n",iteration);
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
				{
					componentNumber[i] = componentNumber[componentNumber[i]];
					// if(i == 239243)
					// 	printf("componentNumber[%d] = %d\n",i,componentNumber[i]);
				}
			}

			// #pragma omp single 
			// {
			// 	printf("shortcutting completed\n");
			// }

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
			head[tid] = newHead;

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
				// if(affectedLCAFlag[parent[i]]==false) continue;
				// if(parent[i] == 174537)
				// 	printf("Marked as cutvertex by edge %d,%d \n", i,parent[i]);
				isCutVertex[parent[i]] += true; 
				// if(affectedLCAFlag[parent[i]] == false)
				// 	printf("LCA vertex %d is not marked as affected \n",parent[i]);
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




