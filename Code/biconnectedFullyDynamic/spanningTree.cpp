#include "spanningTree.h"
#include "unweightedGraph.h"
#include <vector>
#include <queue>
#include <iostream>

void spanningTree(unweightedGraph &G, int* parent)
{
	vector<bool> visited(G.totalVertices,false);
	// std::cout<<G.totalVertices<<" "<<G.totalEdges<<endl;
	// for(int i=0;i<G.totalVertices;i++)
	// {
	// 	if(visited[i]) continue;
	// 	root.push_back(i); 
	// 	bfs(G,i,parent,visited);
	// }
	bfs(G,G.root,parent,visited);
}

void bfs(unweightedGraph &G, int root, int* parent, vector<bool> &visited)
{
	queue<int> q;
	parent[root] = root;
	q.push(root);
	visited[root] = true;
	while(q.size()>0)
	{
		int curr = q.front(); q.pop();
		
		for(int i=G.offset[curr]; i< G.offset[curr+1];i++)
		{
			int next = G.neighbour[i];
			if(visited[next] == false)
			{
				visited[next] = true;
				parent[next] = curr;
				q.push(next);
			}
		}
	}
}