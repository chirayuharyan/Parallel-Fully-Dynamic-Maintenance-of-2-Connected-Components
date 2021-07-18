#ifndef __BFS_H__
#define __BFS_H__

//#define DEBUG
#define THREAD_QUEUE_SIZE 1024
#include "unweightedGraph.h"


struct vertex_set {
  // # of vertices in the set
  int count;
  // max size of buffer vertices 
  int max_vertices;
  // array of vertex ids in set
  int *vertices;
};


void bfs_top_down(unweightedGraph &G, int root ,int* parent, int* level, vertex_set* frontier, vertex_set* new_frontier);

void bfs_bottom_up(unweightedGraph &G, int root ,int* parent, int* level, vertex_set* frontier, vertex_set* new_frontier);

void bfs_hybrid(unweightedGraph &G, int root ,int* parent, int* level, vertex_set* frontier, vertex_set* new_frontier);

void vertex_set_init(vertex_set* list, int count);

double bfs(unweightedGraph &G, int root, int* parents, int* levels);



// void bfs_bottom_up(Graph graph, solution* sol);
// void bfs_hybrid(Graph graph, solution* sol);


inline void add_to_queue(int* thread_queue, int& thread_queue_size, 
                         int* queue_next, int& queue_size_next, int vert);
inline void empty_queue(int* thread_queue, int& thread_queue_size, 
                        int* queue_next, int& queue_size_next);

inline void add_to_queue(int* thread_queue, int& thread_queue_size, 
                         int* queue_next, int& queue_size_next, int vert)
{
  thread_queue[thread_queue_size++] = vert;

  if (thread_queue_size == THREAD_QUEUE_SIZE)
    empty_queue(thread_queue, thread_queue_size, 
                queue_next, queue_size_next);
}

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

#endif