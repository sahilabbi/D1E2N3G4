#ifndef __GRAPH_SET_H__

#define __GRAPH_SET_H__

#include "graph.h"
#include "isomorphism.h"
#include <gsl/gsl_vector.h>

typedef struct
{
    gsl_vector * eigenvalues;
    graph_t ** graphs;
    size_t num_graphs;
} isospectral_group;

//Graph sets are sorted by eigenvalues
typedef struct
{
    isospectral_group * iso_groups;
    size_t size;
} graph_set;

graph_set * graph_set_alloc();
void delete_graph_set(graph_set * gs);
void print_graph_set(graph_set * gs);

//returns a pointer to the isospectral group containing eigenval
//returns NULL if gs contains no such isospectral group
isospectral_group * find_eigenvalues(graph_set * gs, gsl_vector * eigenval);

void insert_graph(graph_set * gs, graph_t * g);
bool check_isomorphism(graph_set * gs, graph_t * g);


#endif
