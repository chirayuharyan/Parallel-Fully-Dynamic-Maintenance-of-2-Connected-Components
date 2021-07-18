#ifndef SPANNINGTREE_H
#define SPANNINGTREE_H
#include "unweightedGraph.h"
#include <vector>
#include <omp.h>
#include "bfs.h"


void spanningTree(unweightedGraph &G, int* parent,int *level, bool* visited, omp_lock_t* lVisited, vertex_set* frontier, vertex_set* new_frontier);

void bfs(unweightedGraph &G, int root, int* parent, bool* visited, omp_lock_t* lVisited, vector<int> &level);
#endif