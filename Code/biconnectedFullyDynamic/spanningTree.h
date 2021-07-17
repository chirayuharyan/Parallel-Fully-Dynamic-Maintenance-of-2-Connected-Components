#ifndef SPANNINGTREE_H
#define SPANNINGTREE_H
#include "unweightedGraph.h"
#include <vector>

void spanningTree(unweightedGraph &G, int* parent);
void bfs(unweightedGraph &G, int root, int* parent,vector<bool> &visited);
#endif 