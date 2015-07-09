#ifndef __GRAPH__

#define __GRAPH__

#include <stdbool.h>

#define N 32
#define K 5

typedef struct
{
	bool adj[N][N];
} graph_t;

void copy_graph(graph_t * dst, graph_t * src);
//gets the next graph in sequence for the brute force search
//returns true if the graph has reached its final state
//and therefore cannot be incremented
//returns false otherwise
bool increment_graph(graph_t * graph);
bool is_regular(graph_t * graph);


#endif

