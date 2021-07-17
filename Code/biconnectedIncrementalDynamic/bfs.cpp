#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>
#include "unweightedGraph.h"
#include <iostream>

#define NOT_VISITED_MARKER -1

#define ALPHA 15.0
#define BETA 24.0

using namespace std;

void vertex_set_clear(vertex_set* list) 
{
    list->count = 0;
}

void vertex_set_init(vertex_set* list, int count) 
{
    list->max_vertices = count;
    list->vertices = (int*)malloc(sizeof(int) * list->max_vertices);
    vertex_set_clear(list);
}

void prepareFrontier(vertex_set* frontier,int *level, int currLevel, int n)
{
    #pragma omp parallel for
    for(int i=0; i<n ; i++)
    {
        if(level[i] == currLevel)
        {
            int index = __sync_fetch_and_add(&frontier->count, 1);
            frontier->vertices[index] = i;
        }
    }
}

void top_down_step( unweightedGraph &G, int* parent, int* level, vertex_set* frontier, vertex_set* new_frontier)
{

    // double time1,time2;

    // time1=omp_get_wtime();

    #pragma omp parallel for
    for (int i=0; i<frontier->count; i++) 
    {

        int node = frontier->vertices[i];

        int start_edge = G.offset[node];
        int end_edge = G.offset[node+1];
        
        for (int neighbor=start_edge; neighbor<end_edge; neighbor++) 
        {
            int outgoing = G.neighbour[neighbor];
            // if(__sync_bool_compare_and_swap(&level[outgoing],NOT_VISITED_MARKER,level[node] + 1))
            if (level[outgoing] == NOT_VISITED_MARKER) 
            {
                level[outgoing] = level[node] + 1;
                parent[outgoing] = node;
                // int index;
                // #pragma omp critical
                //     index = new_frontier->count++;
                int index = __sync_fetch_and_add(&new_frontier->count, 1);
                new_frontier->vertices[index] = outgoing;
            }
        }
    }

    // time2=omp_get_wtime();
    // cout<< "Time taken for top down step: "<< ((time2-time1)*1000) << " milliseconds" << endl;
    // cout<<endl;
}

void bfs_top_down(unweightedGraph &G, int root ,int* parent, int* level, vertex_set* frontier, vertex_set* new_frontier)
{
    

    frontier->vertices[frontier->count++] = root;
    level[root] = 0;
    parent[root] = root;

    while (frontier->count != 0) 
    {

        vertex_set_clear(new_frontier);

        top_down_step(G, parent, level, frontier, new_frontier);

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }

   
}



void bottom_up_step(unweightedGraph &G, int* parent, int* level, vertex_set* frontier, vertex_set* new_frontier, bool &flag, int currLevel, int &frontierSize)
{
    // double time1,time2;

    // time1=omp_get_wtime();

    int n = G.totalVertices;
    #pragma omp parallel for reduction(+:frontierSize)
    for(int i=0;i<n;i++)
    {
        if(level[i] != NOT_VISITED_MARKER) continue;
        
        int start_edge = G.offset[i];   //as this graph is undirected all the edges incident are outgoing as well as incoming. 
        int end_edge = G.offset[i+1];   //hence using the same offset array for top down and bottom up

        for (int neighbor=start_edge; neighbor<end_edge; neighbor++) 
        {
            int incoming = G.neighbour[neighbor];
            if (level[incoming] == currLevel) 
            {
                frontierSize++;
                level[i] = level[incoming] + 1;
                parent[i] = incoming;    
                flag += true;
                break;
            }   
        }
    } 
    // time2=omp_get_wtime();
    // cout<< "Time taken for bottom up step: "<< ((time2-time1)*1000) << " milliseconds" << endl;
    // cout<<endl;  
}


void bfs_bottom_up(unweightedGraph &G, int root ,int* parent, int* level, vertex_set* frontier, vertex_set* new_frontier)
{
    level[root] = 0;
    parent[root] = root;
    bool flag = true;
    int currLevel = 0;

    while (flag) 
    {
        flag = false;
        int x = 0;
        bottom_up_step(G, parent, level, frontier, new_frontier, flag, currLevel,x);

        currLevel++;

    }
}

void bfs_hybrid(unweightedGraph &G, int root ,int* parent, int* level, vertex_set* frontier, vertex_set* new_frontier)
{
    frontier->vertices[frontier->count++] = root;
    level[root] = 0;
    parent[root] = root;
    bool flag = true;
    int currLevel = 0;
    int mode = 0; // mode 0 means top down and 1 means bottom up

    while (flag) 
    {
        flag = false;
        
        if(mode == 0)
        {
            vertex_set_clear(new_frontier);

            top_down_step(G, parent, level, frontier, new_frontier);

            // swap pointers
            vertex_set* tmp = frontier;
            frontier = new_frontier;
            new_frontier = tmp;

            if(frontier->count > 0) flag = true;
            if(frontier->count >= 0.02*G.totalVertices)
                mode = 1;
        }
        else
        {
            int frontierSize = 0;
            bottom_up_step(G, parent, level, frontier, new_frontier, flag, currLevel, frontierSize);
            if(frontierSize < 0.02*G.totalVertices)
            {
                mode = 0;
                prepareFrontier(frontier,level,currLevel+1,G.totalVertices);
            }
        }        
        
        currLevel++;

    }
}


double bfs(unweightedGraph &G, int root, int* parents, int* levels)
{
    int num_verts = G.totalVertices;
    double avg_out_degree = G.totalEdges/(double)G.totalVertices;

    int* queue = new int[num_verts];
    int* queue_next = new int[num_verts];
    int queue_size = 0;  
    int queue_size_next = 0;

    queue[0] = root;
    queue_size = 1;
    parents[root] = root;
    levels[root] = 0;
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

                    unsigned out_degree = out_degree(G, vert);
                    int* outs = out_vertices(G, vert);
                    for (unsigned j = 0; j < out_degree; ++j)
                    {      
                        int out = outs[j];
                        if (levels[out] < 0)
                        {
                            levels[out] = level;
                            parents[out] = vert;
                            ++local_num_descs;
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
                        unsigned out_degree = out_degree(G, vert);
                        int* outs = out_vertices(G, vert);
                        for (unsigned j = 0; j < out_degree; ++j)
                        {
                            int out = outs[j];
                            if (levels[out] == prev_level)
                            {
                                levels[vert] = level;
                                parents[vert] = out;
                                ++local_num_descs;
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

    return time2 - time;
}