#ifndef __GRAPH_SET_H__

#define __GRAPH_SET_H__

#include "graph.h"
#include "isomorphism.h"
#include "comp_graph.h"
#include <gsl/gsl_vector.h>

typedef struct
{
    gsl_vector * eigenvalues;
    comp_graph ** graphs;
    size_t num_graphs;
} isospectral_group;

//Graph sets are sorted by eigenvalues
typedef struct
{
    isospectral_group * iso_groups;
    size_t size;
    size_t amount_alloced;
} graph_set;

graph_set * graph_set_alloc();
void delete_graph_set(graph_set * gs);
void print_graph_set(graph_set * gs);
void print_graph_set_compressed(graph_set * gs);

//returns a pointer to the isospectral group containing eigenval
//returns NULL if gs contains no such isospectral group
isospectral_group * find_eigenvalues(graph_set * gs, gsl_vector * eigenval);
//comparision functions for find_eigenvalues
int vector_compare(const void * a, const void * b);
int isospectral_group_compare(const void * a, const void * b);

void insert_graph(graph_set * gs, graph_t * g);
bool check_isomorphism(graph_set * gs, graph_t * g);


#endif
