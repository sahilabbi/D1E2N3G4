#ifndef __ISOMORPHISM__

#define __ISOMORPHISM__
#include <stdbool.h>

#include "graph.h"

int * eigenvalues(graph_t * g);
int * matrix_vector_mult(graph_t * g, int * v);

bool accurate_isomorphism(graph_t * g1, graph_t * g2);


#endif
