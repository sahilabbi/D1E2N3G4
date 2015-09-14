#ifndef __GRAPH__

#define __GRAPH__

#include <stdbool.h>

//#define NUM_NODES 8
//#define NODE_DEGREE 3

typedef struct
{
	int adj[NUM_NODES][NUM_NODES];
} graph_t;

void copy_graph(graph_t * dst, graph_t * src);
void print_graph(graph_t * g);
//gets the next graph in sequence for the brute force search
//returns true if the graph has reached its final state
//and therefore cannot be incremented
//returns false otherwise
bool increment_graph(graph_t * graph);
bool is_regular(int * degree);

graph_t * initial_graph();
graph_t * first_graph();


#endif

