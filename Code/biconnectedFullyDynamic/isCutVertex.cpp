#include"isCutVertex.h"
#include<queue>
#include<iostream>
#include<unordered_map>
#include<omp.h>
#include"threads.h"
using namespace std;
 
#define THREAD_QUEUE_SIZE 1024
#define ALPHA 15.0
#define BETA 24.0

bool bfs(dynamicGraph &G, unordered_map<int,int> &myParent, int root, int* unfinishedTraversal, int lca, int* parent)
{
	unordered_map<int,bool> visited;
	int totalUnfinished = 0;
	queue<int> q;
	q.push(root);
	myParent[root] = root;

    // bool isThis = false;
    // if(lca == 86976)
    // {
    //     isThis = true;
    // }


	while(q.size()>0)
	{
		int curr = q.front(); q.pop();
       
		while(q.size()>0 && visited.find(curr) != visited.end())
		{
			curr = q.front(); q.pop();
		}
		if(visited.find(curr) != visited.end()) break;

		int x = lca;
		int y = curr;
		if(x>y)
		{
			x = curr;
			y = lca;
		}
		// cout<<"Checking unfinishedTraversal of tree edge "<<x<<","<<y<<endl;
		if(parent[curr]==lca)
		totalUnfinished += unfinishedTraversal[curr];

		visited[curr] = true;
		map<int,unweightedEdge>::iterator it;
		for(it = G.adjacencyList[curr].begin(); it != G.adjacencyList[curr].end(); it++)
		{

			int next = G.adjacencyList[curr][it->first].v;
			if(visited.find(next) != visited.end()) continue;
			myParent[next] = curr;
			q.push(next);
		}

	}

	// if(isThis)
 //    {
 //        cout<<"TOtal unfinished of 86976 is"<<totalUnfinished<<endl;
 //    }
	if(totalUnfinished == 0)
		return true;
	return false;
}

bool isCutVertex(dynamicGraph &G, int* unfinishedTraversal, int lca, int* parent, int &totalComponentsGlobal)
{
	bool result = false; //if any one component is not safe it is a cutVertex
	map<int,map<int,unweightedEdge>>::iterator it;
	unordered_map<int,int> myParent;

	bool isRoot = (parent[lca] == lca);
	int totalComponents = 0;
	
	for(it = G.adjacencyList.begin(); it!=G.adjacencyList.end(); it++)
	{
		if(myParent.find(it->first) != myParent.end())	continue;
		result += bfs(G,myParent,it->first,unfinishedTraversal,lca,parent); //returns true if this component is not safe
		totalComponents++;
		// if(isRoot && totalComponents > 1) return true;
	}
    totalComponentsGlobal += totalComponents;

	// cout<<"For LCA vertex "<<lca<<" total components are "<<totalComponents<<endl;
	if(isRoot == false)
		return result;
	else
		return (totalComponents>1);
}



/* ------------------------------------- 2 imp -----------------------*/

inline void empty_queue(int* thread_queue, int& thread_queue_size, 
                        int* queue_next, int& queue_size_next)
{
	int start_offset;

	#pragma omp atomic capture
		start_offset = queue_size_next += thread_queue_size;

	start_offset -= thread_queue_size;
	for (int i = 0; i < thread_queue_size; ++i)
		queue_next[start_offset + i] = thread_queue[i];
	thread_queue_size = 0;
}

inline void add_to_queue(int* thread_queue, int& thread_queue_size, 
                         int* queue_next, int& queue_size_next, int vert)
{
  thread_queue[thread_queue_size++] = vert;

  if (thread_queue_size == THREAD_QUEUE_SIZE)
    empty_queue(thread_queue, thread_queue_size, 
                queue_next, queue_size_next);
}



double bfs(dynamicGraph &G, int root, int* parents, int* levels, int num_verts, bool* isSafe, int* unfinishedTraversal, int* spanningTree)
{
    // double avg_out_degree = G.totalEdges/(double)G.totalVertices;
	double avg_out_degree = 5; //hardcoded

    int* queue = new int[num_verts];
    int* queue_next = new int[num_verts];
    int queue_size = 0;  
    int queue_size_next = 0;

    queue[0] = root;
    queue_size = 1;
    parents[root] = root;
    levels[root] = 0;

    bool safeComponent = (unfinishedTraversal[root]>0);

   

    // level_queues[0] = new int[1];
    // level_queues[0][0] = root;
    // level_counts[0] = 1;

    int level = 1;
    int num_descs = 0;
    int local_num_descs = 0;
    bool use_hybrid = false;
    bool already_switched = false;

    double time = omp_get_wtime();

    #pragma omp parallel
    {
        int thread_queue[ THREAD_QUEUE_SIZE ];
        int thread_queue_size = 0;

        while (queue_size)
        {
            if (!use_hybrid)
            {
                #pragma omp for schedule(guided) reduction(+:local_num_descs) nowait
                for (int i = 0; i < queue_size; ++i)
                {
                    int vert = queue[i];

                    // unsigned out_degree = out_degree(G, vert);
                    // unsigned out_degree = G.adjacencyList1[vert].size();
                    // int* outs = out_vertices(G, vert);
                    // for (unsigned j = 0; j < out_degree; ++j)
                    for (auto it  = G.adjacencyList1[vert].begin(); it != G.adjacencyList1[vert].end(); it++)
                    {      
                        int out = it->first;
                        if (levels[out] < 0)
                        {
                            levels[out] = level;
                            parents[out] = root;
                            ++local_num_descs;
                            if(unfinishedTraversal[out] > 0)
                                safeComponent += true;
                            if(spanningTree[root] != spanningTree[out])
                            {
                                cout<<"Problem parent dont match root = "<<root<<" out = "<<out<<" parent[root]= "<<spanningTree[root]<<" parent[out]= "<<spanningTree[out]<<endl;
                            }
                            add_to_queue(thread_queue, thread_queue_size, 
                            queue_next, queue_size_next, out);
                        }
                    }
                }
            }
            else
            {
                int prev_level = level - 1;

                #pragma omp for schedule(guided) reduction(+:local_num_descs) nowait
                for (int vert = 0; vert < num_verts; ++vert)
                {
                    if (levels[vert] < 0)
                    {
                        // unsigned out_degree = out_degree(G, vert);
                        // int* outs = out_vertices(G, vert);
                        for (auto it  = G.adjacencyList1[vert].begin(); it != G.adjacencyList1[vert].end(); it++)
                    {
                            int out = it->first;
                            if (levels[out] == prev_level)
                            {
                                levels[vert] = level;
                                parents[vert] = root;
                                ++local_num_descs;
                                if(unfinishedTraversal[vert] > 0)
                                    safeComponent += true;
                                if(spanningTree[root] != spanningTree[vert])
                                {
                                    cout<<"Problem parent dont match root = "<<root<<" out = "<<out<<" parent[root]= "<<spanningTree[root]<<" parent[out]= "<<spanningTree[out]<<endl;
                                }
                                add_to_queue(thread_queue, thread_queue_size, 
                                queue_next, queue_size_next, vert);
                                break;
                            }
                        }
                    }
                }
            }

            empty_queue(thread_queue, thread_queue_size, queue_next, queue_size_next);
            #pragma omp barrier

            #pragma omp single
            { 
                num_descs += local_num_descs;

                if (!use_hybrid)
                {  
                    double edges_frontier = (double)local_num_descs * avg_out_degree;
                    double edges_remainder = (double)(num_verts - num_descs) * avg_out_degree;
                    if ((edges_remainder / ALPHA) < edges_frontier && edges_remainder > 0 && !already_switched)
                        use_hybrid = true;
                }
                else
                {
                    if ( ((double)num_verts / BETA) > local_num_descs  && !already_switched)
                    {
                        use_hybrid = false;
                        already_switched = true;
                    }
                }
                
                local_num_descs = 0;
                // ++num_levels;

                queue_size = queue_size_next;
                queue_size_next = 0;
                int* temp = queue;
                queue = queue_next;
                queue_next = temp;
                ++level;

            } // end single

        }
    } // end parallel
    double time2 = omp_get_wtime();

    delete [] queue;
    delete [] queue_next;

    isSafe[root] += safeComponent;
    return time2 - time;
}

void bfsSerial(dynamicGraph &G, int root, int* parents, int* levels, int num_verts, bool* isSafe, 
                int* unfinishedTraversal, int* spanningTree, int* q)
{
    // queue<int> q;
    int front = 0, tail = 0;
    // q.push(root);
    q[tail] = root;
    tail++;
    parents[root] = root;
    int level = 1;
    bool safeComponent = (unfinishedTraversal[root]>0);

    // while(q.size()>0)
    while(front<tail)    
    {
        int totalQueue = tail - front;
        for(int i=front; i < totalQueue+front; i++)
        {

            int curr = q[i];
    
            // cout<<"Checking unfinishedTraversal of tree edge "<<x<<","<<y<<endl;
            if(unfinishedTraversal[curr] > 0)
                safeComponent += true;

            if(levels[curr] != -1) continue;

            levels[curr] = level;
            
            for(auto it  = G.adjacencyList1[curr].begin(); it != G.adjacencyList1[curr].end(); it++)
            {
                int next = it->first;
                if(levels[next] != -1) continue;
                parents[next] = root;
                q[tail] = next;
                tail++;
                // q.push(next);
            }

        }
        front = front + totalQueue;
        level++;
    }

    isSafe[root] += safeComponent;



}


void isCutVertexAll(dynamicGraph &G, int* unfinishedTraversal, int* parent, int* isLCA,
                    bool *isCutVertex, int totalVertices, int root, int* isBase, int* &localParent, bool* &isSafe,
                    int* &level)
{
	localParent = new int[totalVertices];
	isSafe = new bool[totalVertices];
	level = new int[totalVertices];
	#pragma omp parallel 
	{
		#pragma omp  for nowait
		for(int i=0;i<totalVertices;i++)
		{
			localParent[i] = -1;
		}
		#pragma omp for nowait
		for(int i=0;i<totalVertices;i++)
		{
			isSafe[i] = false;
		}
        #pragma omp for nowait
        for(int i=0;i<totalVertices;i++)
        {
            if(parent[i] ==i)
                root = i;
        }

        #pragma omp for nowait
        for(int i=0;i<totalVertices;i++)
        {
            isCutVertex[i] = false;
        }

		#pragma omp for
		for(int i=0;i<totalVertices;i++)
		{
			level[i] = -1;
		}
	}

    int totalComponentsRoot = 0;
	for(int i=0;i<totalVertices;i++)
	{
		if(isBase[i] == 0 || localParent[i] != -1)
			continue;
        if(parent[i]==root)
        {
            totalComponentsRoot++;
        }
		bfs(G, i, localParent, level,totalVertices,isSafe,unfinishedTraversal,parent);
        // if(!isSafe[i])
        // {
        //     isCutVertex[parent[i]] += true;
        //     break;
        // }
	}

	
	#pragma omp parallel for
	for(int i=0;i<totalVertices;i++)
	{
    
		if(localParent[i] == i && !(isSafe[i]) )
		{
			isCutVertex[parent[i]] += true;
		}
	}
    isCutVertex[root] = (totalComponentsRoot>1);

    // std::cout<<"Total Components = "<<totalComponents<<endl;

}


// void isCutVertexAffectedOld(dynamicGraph &G, int* unfinishedTraversal, int* parent, int* isLCA, bool* isCutVertex,
//                         int totalVertices, int root, int* isBase, int* &localParent, bool* &isSafe, int* &level,
//                         vector<vector<int>> &children, vector<vector<int>> &affectedLCA, int* postProcessed,
//                         int batchNumber)
// {
//     int totalComponentsRoot = 0;
//     for(int i=0;i<NUM_THREADS;i++)
//     {
//         int currChild;
//         for(int j = 0;j<affectedLCA[i].size();j++)
//         {
//             int currAffectedLCA = affectedLCA[i][j];
//             isCutVertex[currAffectedLCA] = false;
//         }
//     }
//     // #pragma omp parallel for
//     for(int i=0;i<NUM_THREADS;i++)
//     {
//         int currChild;
//         for(int j = 0;j<affectedLCA[i].size();j++)
//         {
//             int currAffectedLCA = affectedLCA[i][j];
//             if(postProcessed[currAffectedLCA] == batchNumber) continue;
//             postProcessed[currAffectedLCA] = batchNumber;
//             // isCutVertex[currAffectedLCA] = false;

//             for(int k=0;k<children[currAffectedLCA].size();k++)
//             {
//                 currChild = children[currAffectedLCA][k]; 
//                 localParent[currChild] = -1;
//                 isSafe[currChild] = false;
//                 level[currChild] = -1;
//             }
            
//             for(int k=0;k<children[currAffectedLCA].size();k++)
//             {
//                 currChild = children[currAffectedLCA][k]; 
//                 if(isBase[currChild] == 0 || localParent[currChild] != -1)
//                     continue;
//                 if(parent[currChild]==root)
//                 {
//                     __sync_fetch_and_add(&totalComponentsRoot,1);
//                 }
//                 bfsSerial(G, currChild, localParent, level,totalVertices,isSafe,unfinishedTraversal,parent,q[tid]);
//             }

//             for(int k=0;k<children[currAffectedLCA].size();k++)
//             {
//                 currChild = children[currAffectedLCA][k];
                
//                 if(localParent[currChild] == currChild && !(isSafe[currChild]))
//                 {
//                     isCutVertex[parent[currChild]] += true;
//                 }
//             }

//         }
//     }
    
//     if(postProcessed[root] == batchNumber)
//         isCutVertex[root] = (totalComponentsRoot>1);
// }


// void isCutVertexAffected(dynamicGraph &G, int* unfinishedTraversal, int* parent, int* isLCA, bool* isCutVertex,
//                         int totalVertices, int root, int* isBase, int* &localParent, bool* &isSafe, int* &level,
//                         vector<vector<int>> &children, int* affectedLCA, int* postProcessed,
//                         int batchNumber, int totalAffectedLCA, int** q, unweightedGraph &OG)
void isCutVertexAffected(unweightedGraph &OG, unweightedGraph &G, SpanningTree &T, Biconnected &B, LCAGraph &L,
                        int* postProcessed, int batchNumber)
{
    // printf("Affected LCA list head = %d\n ",L.headAffectedLCA);
    set<int> affectedLCASet;
    
    for(int i = 0; i < L.headAffectedLCA;i++)
    {
        affectedLCASet.insert(L.affectedLCA[i]);
    }
    
    L.headAffectedLCA = 0;
    
    //Step 2: Convert affectedLCA set into a list

    int *affectedLCAList = new int[affectedLCASet.size()];
    int pointer = 0;    
    for(auto it = affectedLCASet.begin(); it!=  affectedLCASet.end(); it++)
    {
        affectedLCAList[pointer] = *it;
        pointer++;
    }


    //Step 3: Collect Partially Affected LCAs
    set<int> partiallyAffectedLCASet;
    for(int i = 0; i < L.headPartiallyAffectedLCA; i++)
    {
        partiallyAffectedLCASet.insert(L.partiallyAffectedLCA[i]);
    }
    L.headPartiallyAffectedLCA = 0;

    //Step 4: Convert PartiallyAffected LCA set into list

    int *partiallyAffectedLCAList = new int[partiallyAffectedLCASet.size()];
    pointer = 0;
    for(auto it = partiallyAffectedLCASet.begin(); it!=  partiallyAffectedLCASet.end(); it++)
    {
        partiallyAffectedLCAList[pointer] = *it;
        pointer++;
    }

    // time2=omp_get_wtime();
    // cout<<"Time For constructing list of affected and partially affected LCA "<<(time2-time1)*1000<<" ms"<<endl;
    
    // time1 = omp_get_wtime();
    
    // Step 5: Work for partially affected LCAs
    int totalPartiallyAffectedLCA = partiallyAffectedLCASet.size();
    
    #pragma omp parallel for schedule(guided)
    for(int i = 0; i < totalPartiallyAffectedLCA; i++)
    {
        int currLCA,currChild;
        currLCA = partiallyAffectedLCAList[i];
        
        B.cutVertex[currLCA] = false;
        
        int startNeighbour = OG.offset[currLCA];
        int endNeighbour = OG.offset[currLCA+1]; 

        for(int j = startNeighbour; j < endNeighbour;j++)
        {
            currChild = OG.neighbour[j];
            if(T.parent[currChild]!=currLCA)
                continue;
            L.isSafe[currChild] = B.unfinishedTraversal[currChild] > 0;
        }


        for(int j = startNeighbour; j < endNeighbour;j++)
        {
            // currChild = children[currLCA][j];
            currChild = OG.neighbour[j];
            if(T.parent[currChild]!=currLCA)
                continue;
            L.isSafe[L.localParent[currChild]] += L.isSafe[currChild];
        }

        
        for(int j = startNeighbour; j < endNeighbour;j++)
        {
            // currChild = children[currLCA][j];
            currChild = OG.neighbour[j];
            if(T.parent[currChild]!=currLCA)
                continue;
            if(L.localParent[currChild] == currChild && (!L.isSafe[currChild]))
                B.cutVertex[currLCA] = true;
        }
        
    }



    int totalComponentsRoot = 0, totalAffectedLCA = affectedLCASet.size();
    bool isRootAffected = false;
    // printf("totalAffectedLCA = %d\n",totalAffectedLCA);

    #pragma omp parallel for schedule(guided)
    for(int i=0;i<totalAffectedLCA;i++)
    {
        int currAffectedLCA = affectedLCAList[i];
        if(currAffectedLCA == G.root)
            isRootAffected = true;
        B.cutVertex[currAffectedLCA] = false;
        
    }
    
    #pragma omp parallel for schedule(guided)
    for(int i=0;i<totalAffectedLCA;i++)
    {
        int tid = omp_get_thread_num();
        int currChild;
        int currAffectedLCA = affectedLCAList[i];

        int startNeighbour = OG.offset[currAffectedLCA];
        int endNeighbour = OG.offset[currAffectedLCA+1]; 
        
        for(int k = startNeighbour; k < endNeighbour; k++)
        {
            currChild = OG.neighbour[k];
            if(T.parent[currChild]!=currAffectedLCA)
                continue;
            L.localParent[currChild] = -1;
            L.isSafe[currChild] = false;
            L.level[currChild] = -1;
        }
            
        for(int k = startNeighbour; k < endNeighbour; k++)
        {
            currChild = OG.neighbour[k];
            
            if(T.parent[currChild]!=currAffectedLCA)
                continue;
            if(B.isBase[currChild] == 0 || L.localParent[currChild] != -1)
                continue;
            
            
            if(T.parent[currChild]==G.root)
            {
                __sync_fetch_and_add(&totalComponentsRoot,1);
            }

            bfsSerial(L.globalGraph, currChild, L.localParent, L.level, G.totalVertices, L.isSafe, B.unfinishedTraversal, T.parent, L.myQueue[tid]);
            
        }

        for(int k = startNeighbour; k < endNeighbour; k++)
        {
            currChild = OG.neighbour[k];
            if(T.parent[currChild]!=currAffectedLCA)
                continue;
            
            if(L.localParent[currChild] == currChild && !(L.isSafe[currChild]))
            {
                B.cutVertex[T.parent[currChild]] += true;
            }
        }

       
    }
    
    if(isRootAffected)
        B.cutVertex[G.root] = (totalComponentsRoot>1);
}