#include<iostream>
#include <omp.h>
#include"unweightedGraph.h"
#include"bfs.h"
#include <sstream>
#include <unordered_map>
#include <set>
#include <utility>
using namespace std;

void readGraph(ifstream &inputGraph, unweightedGraph &G)
{
	string line;
	int totalVertices,totalEdges;
	for(int i = 0;i<2;i++)
		getline(inputGraph,line);
	printf("Starting now\n");
	
	getline(inputGraph,line);
	cout<<line<<endl;
	stringstream ss;
	ss<<line;
	string temp;
	ss>>temp>>temp>>totalVertices>>temp>>totalEdges;
	cout<<totalVertices<<" "<<totalEdges<<endl;
	getline(inputGraph,line);
	

	int u,v,i=0;
	int *U = new int[totalEdges];
	int *V = new int[totalEdges];
	bool flag = false;
	while(!inputGraph.eof())
	{
		// getline(inputGraph,line);
		// cout<<line<<endl;
		// ss<<line;
		inputGraph>>u>>v;
		// cout<<u<<" "<<v<<endl;
		U[i] = u;
		V[i] = v;
		i++;
		if(!flag && (u >= totalVertices || v >= totalVertices))
		{
			flag = true;
		}
	}


	if(flag)
	{
		printf("Renaming is required\n");
		unordered_map<int,int> m;
		int count = 0;
		for(int i=0;i<totalEdges;i++)
		{
			int u = U[i];
			int v = V[i];
			if(m.find(u) == m.end())
			{
				m[u] = count;
				count++;
			}
			if(m.find(v) == m.end())
			{
				m[v] = count;
				count++;
			}
		}

		for(int i=0;i<totalEdges;i++)
		{
			U[i] = m[U[i]];
			V[i] = m[V[i]];
		}
		printf("Renaming complete with count = %d and totalVertices = %d\n",count,totalVertices);
	}


	/* ----------------------Creating set --------------------------------*/

	set<pair<int,int>> edgeSet;
	for(int i=0;i<totalEdges;i++)
	{
		edgeSet.insert(make_pair(U[i],V[i]));
		edgeSet.insert(make_pair(V[i],U[i]));
	}

	delete [] U;
	delete [] V;

	G.edgeSet = edgeSet;
	G.totalVertices = totalVertices;
	G.recreate();
	printf("Old Total Edges = %d\n",G.totalEdges);


	vector<int> root;
	int *parent = new int[G.totalVertices];
	int *levels = new int[G.totalVertices];

	#pragma omp parallel
	{
		#pragma omp for nowait
		for(int i=0;i<G.totalVertices;i++)
		{
			parent[i] = -1;
		}
		#pragma omp for
		for(int i=0;i<G.totalVertices;i++)
		{
			levels[i] = -1;
		}
	}
	spanningTree(G, root, parent,levels);

	if(root.size()>1)
	{
		printf("The graph is broken\n");
		for(unsigned int i=1;i<root.size();i++)
		{
			G.addEdgeSet(root[0],root[i]);
		}
		printf("Total edges added = %d\n",int(root.size()-1));
		G.recreate();
	}
	
	root.clear();
	#pragma omp parallel
	{
		#pragma omp for nowait
		for(int i=0;i<G.totalVertices;i++)
		{
			parent[i] = -1;
		}
		#pragma omp for
		for(int i=0;i<G.totalVertices;i++)
		{
			levels[i] = -1;
		}
	}
	spanningTree(G, root, parent,levels);

	if(root.size()==1)
	{
		printf("The graph is now connected with totalEdges = %d\n",int(G.edgeSet.size()));
	}
	else
	{
		printf("The graph is disconnected with number of roots = %d\n",int(root.size()));
	}

}


int main(int argc, char* argv[])
{
	int totalThreads = 64;
	omp_set_num_threads(totalThreads);


	ifstream inputGraph;
	inputGraph.open(argv[1]);

	ofstream output;
	output.open(argv[2]);

	cout<<"Started"<<endl;
	unweightedGraph G;
	cout<<"Created Graph"<<endl;

	readGraph(inputGraph,G);

	output<<G.totalVertices<<" "<<G.totalEdges<<endl;
	for(auto it = G.edgeSet.begin(); it!=G.edgeSet.end(); it++)
	{
		output<<it->first<<" "<<it->second<<endl;
	}

}