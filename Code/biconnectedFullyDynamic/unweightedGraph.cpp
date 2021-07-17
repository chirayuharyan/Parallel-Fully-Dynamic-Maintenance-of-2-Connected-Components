#include"unweightedGraph.h"
#include<iostream>
#include<vector>
using namespace std; 

unweightedGraph::unweightedGraph(ifstream &edgeList)
{
	edgeList>>totalVertices>>totalEdges;
	// cout<<"Total Vertices :"<<totalVertices<<" TotalEdges :"<<totalEdges<<endl;
	offset = new int[totalVertices+1];
	neighbour = new int[totalEdges];
	root = 0;
	for(int i=0;i<totalVertices;i++)
	{
		vector<edge> temp;
		adjacencyList.push_back(temp);
	}

	vector<int> totalNeighbours(totalVertices, 0);
	// cout<<"Total neighbour size = "<<totalNeighbours.size()<<endl; 
	int *U, *V; //To store edges as we do not want to read from file again
	U = new int[totalEdges];
	V = new int[totalEdges];
	for(int i=0;i<totalEdges;i++)
	{
		int u,v;
		edgeList>>u>>v;
		// cout<<"Edge : "<<u<<" "<<v<<endl;
		totalNeighbours[u]++;
		// cout<<1<<endl;
		U[i]=u;
		V[i]=v;
		edge curr;
		curr.u = u;
		curr.v = v;
		adjacencyList[u].push_back(curr);
		// cout<<2<<endl;
	}

	vector<int> edgeCounter(totalVertices);
	offset[0]=0;
	maxDegree = 0;
	for(int i=0;i<totalVertices;i++)
	{
		offset[i+1]=totalNeighbours[i]+offset[i];
		if(totalNeighbours[i] > maxDegree)
		{
			root = i;
			maxDegree = totalNeighbours[i];
		}
		edgeCounter[i] = offset[i];
	}

	

	for(int i=0;i<totalEdges;i++)
	{
		int u,v;
		u = U[i];
		v = V[i];

		int currIndex = edgeCounter[u];
		edgeCounter[u]++;
		neighbour[currIndex] = v;
	}

	degree = new int[totalVertices];
	for(int i=0;i<totalVertices;i++)
	{
		degree[i] = totalNeighbours[i];
	}


	delete [] U;
	delete [] V;
}

void unweightedGraph::addEdge(int u, int v)
{

	//adding u,v
	totalEdges++;
	int *newNeighbour = new int[totalEdges];
	for(int i=0;i<=u;i++)
	{
		for(int j = offset[i]; j<offset[i+1]; j++)
		{
			newNeighbour[j] = neighbour[j];
		}
	}
	newNeighbour[offset[u+1]] = v;
	for(int i = u+1; i<totalVertices; i++)
	{
		for(int j = offset[i]; j<offset[i+1]; j++)
		{
			newNeighbour[j+1] = neighbour[j];
		}
		offset[i] = offset[i]+1;
	}
	offset[totalVertices]++;
	neighbour = newNeighbour;
	
	edge curr;
	curr.u = u;
	curr.v = v;
	adjacencyList[u].push_back(curr);

	//adding v,u

	totalEdges++;
	int *newNeighbour1 = new int[totalEdges];
	for(int i=0;i<=v;i++)
	{
		for(int j = offset[i]; j<offset[i+1]; j++)
		{
			newNeighbour1[j] = neighbour[j];
		}
	}
	newNeighbour1[offset[v+1]] = u;
	for(int i = v+1; i<totalVertices; i++)
	{
		for(int j = offset[i]; j<offset[i+1]; j++)
		{
			newNeighbour1[j+1] = neighbour[j];
		}
		offset[i] = offset[i]+1;
	}
	offset[totalVertices]++;
	neighbour = newNeighbour1;

	// edge curr;
	curr.u = v;
	curr.v = u;
	adjacencyList[v].push_back(curr);	
}

void unweightedGraph::deleteEdge(int u, int v)
{
	//delete u,v
	// cout<<"Started Delete function"<<endl;
	totalEdges--;
	// cout<<"Total Edges = "<<totalEdges<<endl;

	// cout<<"Creating New array"<<endl;
	int *newNeighbour = new int[totalEdges];
	// cout<<"Created New array "<<endl;

	int point = 0;
	for(int i=0;i<totalVertices;i++)
	{
		// cout<<"Curr Vertex "<<i<<endl;
		for(int j = offset[i]; j<offset[i+1]; j++)
		{
			// cout<<"Current neighbour "<<neighbour[j]<<endl;
			if(i == u && neighbour[j] == v)
			{
				// cout<<"Found the edge u,v to be deleted"<<endl;
				continue;
			}
			newNeighbour[point] = neighbour[j];
			point++;
		}
		if(i>u)
		{
			offset[i]--;
		}
	}
	offset[totalVertices]--;
	neighbour = newNeighbour;

	//delete v,u 
	totalEdges--;
	int *newNeighbour1 = new int[totalEdges];

	point = 0;
	for(int i=0;i<totalVertices;i++)
	{
		// cout<<"Curr Vertex "<<i<<endl;
		for(int j = offset[i]; j<offset[i+1]; j++)
		{
			// cout<<"Current neighbour "<<neighbour[j]<<endl;
			if(i == v && neighbour[j] == u)
			{
				// cout<<"Found the edge v,u to be deleted"<<endl;
				continue;
			}
			newNeighbour1[point] = neighbour[j];
			point++;
		}
		if(i>v)
		{
			offset[i]--;
		}
	}
	offset[totalVertices]--;

	neighbour = newNeighbour1;

	// cout<<"Deleting from adjacencyList u to v"<<endl;

	//delete u,v adjacencyList
	int position = -1;
	for(int i=0;i<adjacencyList[u].size();i++)
	{
		if(adjacencyList[u][i].v == v)
		{
			position = i;
			// cout<<"Found the edge u,v to be deleted in adjacencyList with position"<<position<<endl;
			break;
		}
	}
	if(position == -1)
		cout<<"Position invalid"<<endl;
	adjacencyList[u].erase(adjacencyList[u].begin() + position);


	// cout<<"Deleting from adjacencyList v to u "<<endl;
	//delete v,u adjacencyList
	position = -1;
	for(int i=0;i<adjacencyList[v].size();i++)
	{
		if(adjacencyList[v][i].v == u)
		{
			position = i;
			// cout<<"Found the edge v,u to be deleted in adjacencyList with position"<<position<<endl;
			break;
		}
	}
	if(position == -1)
		cout<<"Position invalid"<<endl;
	adjacencyList[v].erase(adjacencyList[v].begin() + position);


}

void unweightedGraph::addEdgeDegree(int u, int v)
{
	degree[u]++;
	degree[v]++;
}

void unweightedGraph::deleteEdgeDegree(int u, int v)
{
	degree[u]--;
	degree[v]--;
}

void unweightedGraph::print()
{
	cout<<"Total Edges = "<<totalEdges<<endl;
	cout<<"Total Vertices = "<<totalVertices<<endl;
	for(int i=0;i<adjacencyList.size();i++)
	{
		cout<<"Vertex "<<i<<" -> { ";
		for(int j=0;j<adjacencyList[i].size();j++)
		{
			printEdge(adjacencyList[i][j]);
			cout<<", ";
		}
		cout<<"}"<<endl;
	}
}

void unweightedGraph::printCSR()
{
	cout<<"Total Edges = "<<totalEdges<<endl;
	cout<<"Total Vertices = "<<totalVertices<<endl;
	for(int i=0;i<totalVertices;i++)
	{
		int u = i;
		cout<<"Vertex "<<u<<" -> { ";
		for(int j=offset[i]; j<offset[i+1];j++)
		{
			int v = neighbour[j];
			cout<<"( "<<u<<","<<v<<")";
			cout<<", ";
		}
		cout<<"}"<<endl;
	}
}

void printEdge(edge e)
{
	cout<<"( "<<e.u<<","<<e.v<<")";
}
