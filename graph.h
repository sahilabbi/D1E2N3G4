#ifndef __GRAPH__

#define __GRAPH__

#include <stdbool.h>

struct graph_t
{
	bool adj[N][N];
};

void copy_graph(graph_t * dst, graph_t * src);
//gets the next graph in sequence for the brute force search
void increment_graph(graph_t * graph);


#endif

