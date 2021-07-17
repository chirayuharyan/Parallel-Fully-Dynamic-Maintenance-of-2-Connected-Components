#include "spanningTree.h"
#include "unweightedGraph.h"
#include <vector>
#include <queue>
#include <iostream>




void bfs2(unweightedGraph &G, int root, int* parent, bool* visited, omp_lock_t* lVisited)
{
	vector<int>* frontier = new vector<int>;
	vector<int>* nextFrontier = new vector<int>;
	frontier->push_back(root);
	visited[root] = true;
	parent[root] = root;
			
	while(frontier->size()!=0)
	{
		#pragma omp parallel for
		for(int i=0;i<frontier->size();i++)
		{
			int curr,next;
			
			curr = (*frontier)[i];
			for(int j = G.offset[curr]; j < G.offset[curr+1]; j++)
			{
				bool flag = false;
				next = G.neighbour[j];
				// omp_set_lock(&(lVisited[next]));
				if(visited[next] == false)
				{
					visited[next] = true;
					parent[next] = curr;
					flag = true;
				}
				// omp_unset_lock(&(lVisited[next]));
				if(!flag) continue;
				#pragma omp critical
				{
					nextFrontier->push_back(next);				
				}	
			}
		}
		frontier = nextFrontier;
		nextFrontier = new vector<int>;
	}
}

void bfs4(unweightedGraph &G, int root, int* parent, bool* visited, omp_lock_t* lVisited)
{
	vector<bool>* frontier = new vector<bool>(G.totalVertices, false);
	vector<bool>* nextFrontier = new vector<bool>(G.totalVertices, false);
	(*frontier)[root] = true;
	visited[root] = true;
	parent[root] = root;
	bool continueFlag = true;
	while(continueFlag)
	{
		continueFlag = false;
		#pragma omp parallel for
		for(int i=0;i<G.totalVertices;i++)
		{
			if((*frontier)[i] == false) continue;
			int curr,next;
			
			curr = i;
			for(int j = G.offset[curr]; j < G.offset[curr+1]; j++)
			{
				bool flag = false;
				next = G.neighbour[j];
				// omp_set_lock(&(lVisited[next]));
				if(visited[next] == false)
				{
					continueFlag += true;
					visited[next] = true;
					parent[next] = curr;
					flag = true;
				}
				// omp_unset_lock(&(lVisited[next]));
				if(!flag) continue;
				(*nextFrontier)[next] = true;					
			}
		}
		frontier = nextFrontier;
		nextFrontier = new vector<bool>(G.totalVertices,false);
	}
}

void bfs3(unweightedGraph &G, int root, int* parent, bool* visited, omp_lock_t* lVisited)
{
	queue<int> q;
	parent[root] = root;
	q.push(root);
	visited[root] = true;
	while(q.size()>0)
	{
		int total = q.size();
		#pragma omp parallel for
		for(int i=0;i<total;i++)
		{
			int curr;
			#pragma omp critical
			{
				curr = q.front(); q.pop();	
			}
		

			for(int i=G.offset[curr]; i< G.offset[curr+1];i++)
			{
				int next = G.neighbour[i];
				if(visited[next] == false)
				{
					visited[next] = true;
					parent[next] = curr;
					 #pragma omp critical
						q.push(next);
				}
			}
		}
	}
}

void bfs(unweightedGraph &G, int root, int* parent, bool* visited, omp_lock_t* lVisited, int *level)
{
	queue<int> q;
	parent[root] = root;
	q.push(root);
	level[root] = 0;
	// visited[root] = true;
	while(q.size()>0)
	{
		int curr = q.front(); q.pop();
			
		for(int i=G.offset[curr]; i< G.offset[curr+1];i++)
		{
			int next = G.neighbour[i];
			
			if(level[next] == -1)
			{
				// visited[next] = true;
				parent[next] = curr;
				level[next] = level[curr]+1;
				q.push(next);
			}
		}
	}
}




void spanningTree(unweightedGraph &G, int* parent,int *level, bool* visited, omp_lock_t* lVisited, vertex_set* frontier, vertex_set* new_frontier)
{
	bfs_hybrid(G, G.root, parent, level, frontier, new_frontier);
	// bfs(G,G.root,parent,visited,lVisited,level);
	// bfs_top_down(G, G.root , parent, level, frontier, new_frontier);
	// bfs_bottom_up(G, G.root , parent, level, frontier, new_frontier);
}
