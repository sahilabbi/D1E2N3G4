#ifndef __GRAPH_SET_H__

#define __GRAPH_SET_H__

#include "graph.h"
#include "isomorphism.h"
#include <gsl/gsl_vector.h>

//Graph sets are sorted by eigenvalues
typedef struct
{
    gsl_vector ** eigenvalues;
    graph_t *** graphs;
    size_t * graph_set_sizes;
    size_t size;
} graph_set;

graph_set * graph_set_alloc();
void delete_graph_set(graph_set * gs);

//returns the index of the eigenvalues if a matching vector is found
//otherwise returns -1 * the index the vector should be before otherwise
int find_eigenvalues(graph_set * gs, gsl_vector * eigenval);

void insert_graph(graph_set * gs, graph_t * g);
bool check_isomorphism(graph_set * gs, graph_t * g);


#endif
